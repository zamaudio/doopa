/*
    doopa - A fast bam file deduplicator

    Copyright (c) 2018 Damien Zammit <damien@zamaudio.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <inttypes.h>
#include <unordered_map>
#include <map>
#include <tuple>
#include <utility>
#include <functional>
#include <math.h>

#include <openssl/sha.h>

#include "htslib/thread_pool.h"
#include "htslib/hfile.h"
#include "htslib/sam.h"

// Pack into 64 bits:
// chr  start   len
// ff8 7fffffff ffffff
#define PACK_CHRPOSLEN(chr, pos, len) \
    (uint64_t)( (((uint64_t)(chr) & 0x1ff) << 55) | \
                (((uint64_t)(pos) & 0x7fffffff) << 24) | \
                 ((uint64_t)(len) & 0xffffff) \
              )

#define PACK_STARTPOS(pos1, pos2) \
    (uint64_t)( (((uint64_t)(pos1)) << 32) | \
                 ((uint64_t)(pos2)) \
              )

#define EXTRACT_STARTPOS(chrposlen) \
    (uint64_t)( ((uint64_t)(chrposlen >> 24) & 0x7fffffff) )

#define ABS(x)  ((x < 0) ? (-x) : (x))

#define MAX_THREADS	8
#define FRAGMENT_BIN_SIZE 5
#define MAX_FRAGMENT_SIZE 2000

typedef struct {
  uint64_t lo;
  uint64_t hi;
} chrposlen_t;

size_t key_hash(const chrposlen_t& k) {
    uint64_t result;
    const uint8_t *ptr = (const uint8_t *)&k;
    uint8_t sha[SHA256_DIGEST_LENGTH] = {0};

    SHA256(ptr, 16, sha);

    result = *((uint64_t *)sha);
    return result;
}

void error(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    fflush(stderr);
    fprintf(stderr, "doopa: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    va_end(args);
}

/* Calculate the current read's start based on the stored cigar string. */
static int32_t unclipped_start(bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int32_t clipped = 0;
    uint32_t i;

    for (i = 0; i < b->core.n_cigar; i++) {
        char c = bam_cigar_opchr(cigar[i]);

        if (c == 'S' || c == 'H') { // clips
            clipped += bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }

    return b->core.pos - clipped + 1;
}

/* Calculate the current read's end based on the stored cigar string. */
static int32_t unclipped_end(bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int32_t end_pos, clipped = 0;
    int32_t i;

    end_pos = bam_endpos(b);

    // now get the clipped end bases (if any)
    // if we get to the beginning of the cigar string
    // without hitting a non-clip then the results are meaningless
    for (i = b->core.n_cigar - 1; i >= 0; i--) {
        char c = bam_cigar_opchr(cigar[i]);

        if (c == 'S' || c == 'H') { // clips
            clipped += bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }

    return end_pos + clipped;
}

/* Calculate the mate's unclipped start based on position and cigar string from MC tag. */
static int32_t unclipped_other_start(int32_t op, char *cigar) {
    char *c = cigar;
    int32_t clipped = 0;

    while (*c && *c != '*') {
        long num = 0;

        if (isdigit((int)*c)) {
            num = strtol(c, &c, 10);
        } else {
            num = 1;
        }

        if (*c == 'S' || *c == 'H') { // clips
            clipped += num;
        } else {
            break;
        }

        c++;
    }
    return op - clipped + 1;
}

/* Calculate the mate's unclipped end based on start position and cigar string from MC tag.*/
static int32_t unclipped_other_end(int32_t op, char *cigar) {
    char *c = cigar;
    int32_t refpos = 0;
    int skip = 1;

    while (*c && *c != '*') {
        long num = 0;

        if (isdigit((int)*c)) {
            num = strtol(c, &c, 10);
        } else {
            num = 1;
        }

        switch (*c) {
            case 'M':
            case 'D':
            case 'N':
            case '=':
            case 'X':
                refpos += num;
                skip = 0; // ignore initial clips
            break;

            case 'S':
            case 'H':
                if (!skip)
                    refpos += num;
            break;
        }

        c++;
    }
    return op + refpos;
}

/* Create a signature hash of the current read and its pair.
   Uses the unclipped start and end of read and its pair. */
static void make_key(chrposlen_t *key, bam1_t *bam) {
    int32_t chr1, chr2, start1, start2, len1, len2, tmp;
    uint8_t *data;
    char *cig;

    chr1   = bam->core.tid;
    start1 = unclipped_start(bam);
    len1   = unclipped_end(bam) - start1;
    chr2   = 0;
    start2 = 0;
    len2   = 0;

    if (bam->core.mtid >= 0) {
        chr2 = bam->core.mtid;
        if ((data = bam_aux_get(bam, "MC"))) {
            cig = bam_aux2Z(data);
            start2 = unclipped_other_start(bam->core.mpos, cig);
            tmp = unclipped_other_end(bam->core.mpos, cig) - start2;
            len2 = ABS(tmp);
        } else {
            start1 = bam->core.pos;
            len1 = bam_endpos(bam) - start1;
            start2 = bam->core.mpos;
            len2 = len1;
        }
    }
    key->lo = PACK_CHRPOSLEN(chr1, start1, len1);
    key->hi = PACK_CHRPOSLEN(chr2, start2, len2);
}

bool key_equal_to(const chrposlen_t& k1, const chrposlen_t& k2) {
    return key_hash(k1) == key_hash(k2);
}

typedef std::unordered_map<chrposlen_t, std::pair<uint64_t, uint64_t>,
        std::function<size_t(const chrposlen_t&)>,
        std::function<bool(const chrposlen_t&, const chrposlen_t&)> > doopa_t;

typedef std::map<uint64_t, uint64_t> fragment_t;

static inline uint64_t get_qualsum(const bam1_t *b, uint64_t *total, uint64_t *q30)
{
    int i;
    uint64_t sum, q, len;
    uint8_t *qual = bam_get_qual(b);

    len = b->core.l_qseq;
    for (i = sum = 0; i < len; ++i) {
        q = qual[i];
        sum += q;
        if (q30 && q >= 30) {
            (*q30)++;
        }
    }
    if (total)
        *total += len;
    return sum;
}

void print_frag_stats(fragment_t *frag_hist, uint64_t total_fragments)
{
    float half_dist = total_fragments * 0.5f;
    float csum = 0.f;
    float dist = 0.f;
    float mean = 0.f;
    float stdev = 0.f;
    float cstd = 0.f;
    float median = 0.f;
    int found = 0;

    for (fragment_t::iterator frag = frag_hist->begin(); frag != frag_hist->end(); frag++) {
        float mi = (float)frag->first * (float)FRAGMENT_BIN_SIZE + (float)FRAGMENT_BIN_SIZE * 0.5f;
        float curr = mi * (float)frag->second;
        dist += (float)frag->second;
        csum += curr;
        if (!found && (dist >= half_dist)) {
            /* L + ( (n/2 – F) / f ) * w */
            median = ((float)frag->first * (float)FRAGMENT_BIN_SIZE) + ((half_dist - dist) / (float)frag->second) * (float)FRAGMENT_BIN_SIZE;
            found = 1;
        }
    }
    mean = csum / (float)total_fragments;

    for (fragment_t::iterator frag = frag_hist->begin(); frag != frag_hist->end(); frag++) {
        float mi = (float)frag->first * (float)FRAGMENT_BIN_SIZE + (float)FRAGMENT_BIN_SIZE * 0.5f;
        float cmu = mi - mean;
        cstd += (float)frag->second * cmu * cmu;
    }
    stdev = sqrt(cstd / ((float)total_fragments - 1.f));

    error("Mean fragment size: %.4f", mean);
    error("Median fragment size: %.0f", median);
    error("Stdev fragment size: %.4f", stdev);
}

static void dedup_bam(const char *filename, bool stats_only)
{
    htsThreadPool p = {NULL, 0};
    uint64_t total_reads = 0;
    uint64_t paired_reads = 0;
    uint64_t mapped_reads = 0;
    uint64_t bases_above_q30 = 0;
    uint64_t total_bases = 0;
    uint64_t duplicate_reads = 0;
    uint64_t qualsum, existing_qual, fragment_bin;
    chrposlen_t key;
    hts_itr_t *iter;
    bam1_t *b;
    bam_hdr_t *hdr = NULL;
    samFile *out = NULL;
    samFile *in = NULL;
    hts_idx_t *idx = NULL;
    htsFormat _bam;
    hts_parse_format(&_bam, "bam");

    htsFile *fp = hts_open(filename,"r");
    if (!fp) {
        if (errno == ENOEXEC) {
            error("Couldn't understand format of \"%s\"", filename);
            exit(1);
        } else {
            error("Couldn't open \"%s\"", filename);
            exit(1);
        }
    }

    enum htsExactFormat format = hts_get_format(fp)->format;
    if (format != bam) {
        error("File \"%s\" is not a bam file", filename);
        exit(1);
    }

    hts_close(fp);

    in = sam_open(filename, "r");

    if ((idx = sam_index_load(in, filename)) == 0) {
        error("cannot open bam index");
        exit(1);
    }

    doopa_t mp((doopa_t::size_type)1000000, key_hash, key_equal_to);
    fragment_t fragment_histogram;

    out = sam_open("/dev/stdout", "w");

    if (out == NULL) { error("reopening standard output failed"); goto clean; }

    if (!(p.pool = hts_tpool_init(MAX_THREADS))) {
        error("error creating thread pool");
        goto clean;
    }
    hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
    hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);

    hdr = sam_hdr_read(in);
    if (hdr == NULL) {
        errno = 0; error("reading headers from \"%s\" failed", filename);
        goto clean;
    }

    if (!stats_only) {
        if (sam_hdr_write(out, hdr) != 0) {
            error("writing headers to standard output failed");
            goto clean;
        }
    }

    error("Start deduping...");

    iter = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
    b = bam_init1();
    if (b == NULL) { error("can't create record"); exit(1); }
    for(; sam_itr_next(in, iter, b) >= 0; total_reads++) {
        const bam1_core_t *c = &b->core;
        if (c->tid < 0) {
            continue;
        }
        // read must not be secondary, supplementary, unmapped or failed QC
        if (c->flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FQCFAIL)) {
            continue;
        }
        mapped_reads++;
        if (c->mtid >= 0) {
            if ((c->flag & BAM_FPROPER_PAIR) &&
                    c->isize > 0 && c->qual > 30) {
                if (c->isize > MAX_FRAGMENT_SIZE) {
                    fragment_bin = MAX_FRAGMENT_SIZE / FRAGMENT_BIN_SIZE;
                } else {
                    fragment_bin = c->isize / FRAGMENT_BIN_SIZE;
                }
                fragment_histogram[fragment_bin]++;
                paired_reads += 2;
            }
        }
        make_key(&key, b);
        qualsum = get_qualsum(b, &total_bases, &bases_above_q30);
        if (mp.count(key)) {
            // Key exists
            duplicate_reads++;
            existing_qual = std::get<1>(mp[key]);
            if (qualsum > existing_qual) {
                mp[key] = std::make_pair(total_reads, qualsum);
            }
        } else {
            mp[key] = std::make_pair(total_reads, qualsum);
        }
    }
    error("Total bases:\t%lld", total_bases);
    error("Bases above Q30:\t%lld", bases_above_q30);
    error("Total reads:\t%lld", total_reads);
    error("Paired reads:\t%lld", paired_reads);
    error("Mapped reads:\t%lld", mapped_reads);
    error("Duplicate reads:\t%lld", duplicate_reads);
    print_frag_stats(&fragment_histogram, paired_reads / 2);
    error("");

    error("Fragment Histogram:");
    error("Lower\tUpper\tFrequency");
    for (fragment_t::iterator frag = fragment_histogram.begin(); frag != fragment_histogram.end(); frag++) {
        error("%u\t%u\t%u", frag->first * FRAGMENT_BIN_SIZE, frag->first * FRAGMENT_BIN_SIZE + (FRAGMENT_BIN_SIZE - 1), frag->second);
    }

    hts_itr_destroy(iter);

    iter = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
    if (!stats_only) {
        for(total_reads = 0; sam_itr_next(in, iter, b) >= 0; total_reads++) {
            const bam1_core_t *c = &b->core;
            if (c->tid < 0) {
                /* Write unmapped reads as is */
                if (sam_write1(out, hdr, b) < 0) {
                    error("writing to standard output failed");
                    exit(1);
                }
                continue;
            }
            make_key(&key, b);
            qualsum = get_qualsum(b, NULL, NULL);
            if ((qualsum == std::get<1>(mp[key])) && (std::get<0>( mp[key] ) == total_reads)) {
                if (sam_write1(out, hdr, b) < 0) {
                    error("writing to standard output failed");
                    exit(1);
                }
            }
        }
    }
    error("Done");

    bam_destroy1(b);
    hts_itr_destroy(iter);

clean:
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(in);
    if (sam_close(out) < 0) {
        error("could not close output file");
    }
    if (p.pool) hts_tpool_destroy(p.pool);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        error("needs indexed bam file as input");
        return 1;
    }
    if (!strcmp(argv[1], "--statsonly")) {
        if (argc < 3) {
            error("needs indexed bam file as input");
            return 1;
        }
        dedup_bam(argv[2], true);
    } else {
        dedup_bam(argv[1], false);
    }
    return 0;
}

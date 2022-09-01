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
    uint8_t resbuf[32] = {0};

    SHA256(ptr, 16, resbuf);

    result = *((uint64_t *)resbuf);
    return result;
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
        if (q >= 30) {
            (*q30)++;
        }
    }
    *total += len;
    return sum;
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
    int32_t chr, start, stop, len, chr2, start2;
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
        mapped_reads++;
        chr = c->tid;
        start = c->pos;
        stop = bam_endpos(b);
        len = stop - start;
        chr2 = 0;
        start2 = 0;
        if (c->mtid >= 0) {
            chr2 = c->mtid;
            start2 = c->mpos;
            if ((c->flag & BAM_FPROPER_PAIR) == BAM_FPROPER_PAIR &&
                    (c->flag & BAM_FSECONDARY) != BAM_FSECONDARY &&
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
        qualsum = get_qualsum(b, &total_bases, &bases_above_q30);
        if (mp.count({PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, 0)})) {
            // Key exists
            duplicate_reads++;
            existing_qual = std::get<1>(mp[{PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, 0)}]);
            if (qualsum > existing_qual) {
                mp[{PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, 0)}] = std::make_pair(total_reads, qualsum);
            }
        } else {
            mp[{PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, 0)}] = std::make_pair(total_reads, qualsum);
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
            chr = c->tid;
            start = c->pos;
            stop = bam_endpos(b);
            len = stop - start;
            chr2 = 0;
            start2 = 0;
            if (c->mtid >= 0) {
                chr2 = c->mtid;
                start2 = c->mpos;
            }

            if (std::get<0>( mp[{PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, 0)}] ) == total_reads) {
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

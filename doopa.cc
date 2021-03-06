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
#include <tuple>
#include <utility>
#include <functional>

#include "htslib/thread_pool.h"
#include "htslib/hfile.h"
#include "htslib/sam.h"

typedef struct {
  uint64_t lo;
  uint64_t hi;
} chrposlen_t;

size_t key_hash(const chrposlen_t& k) {
    return (k.lo & 0x5a5a5a5a5a5a5a5aull) | (k.hi & 0xa5a5a5a5a5a5a5a5ull);
}

bool key_equal_to(const chrposlen_t& k1, const chrposlen_t& k2) {
    return key_hash(k1) == key_hash(k2);
}

typedef std::unordered_map<chrposlen_t, std::pair<uint64_t, uint64_t>,
        std::function<size_t(const chrposlen_t&)>,
        std::function<bool(const chrposlen_t&, const chrposlen_t&)> > doopa_t;

// Pack into 64 bits:
// chr  start   len
// ff8 7fffffff ffffff
#define PACK_CHRPOSLEN(chr, pos, len) \
    (uint64_t)( (((uint64_t)(chr) & 0x1ff) << 55) | \
                (((uint64_t)(pos) & 0x7fffffff) << 24) | \
                 ((uint64_t)(len) & 0xffffff) \
              )

#define MAX_THREADS	8

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

static void dedup_bam(const char *filename, bool stats_only)
{
    htsThreadPool p = {NULL, 0};
    uint64_t reads = 0;
    uint64_t mapped_reads = 0;
    uint64_t bases_above_q30 = 0;
    uint64_t total_bases = 0;
    uint64_t duplicate_reads = 0;
    uint64_t qualsum, existing_qual;
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
    for(; sam_itr_next(in, iter, b) >= 0; reads++) {
        if (b->core.tid < 0) {
            continue;
        }
        mapped_reads++;
        chr = b->core.tid;
        start = b->core.pos;
        stop = bam_endpos(b);
        len = stop - start;
        chr2 = 0;
        start2 = 0;
        if (b->core.mtid >= 0) {
            chr2 = b->core.mtid;
            start2 = b->core.mpos;
        }
        qualsum = get_qualsum(b, &total_bases, &bases_above_q30);
        if (mp.count({PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, len)})) {
            // Key exists
            duplicate_reads++;
            existing_qual = std::get<1>(mp[{PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, len)}]);
            if (qualsum > existing_qual) {
                mp[{PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, len)}] = std::make_pair(reads, qualsum);
            }
        } else {
            mp[{PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, len)}] = std::make_pair(reads, qualsum);
        }
    }
    error("Total bases:\t%lld", total_bases);
    error("Bases above Q30:\t%lld", bases_above_q30);
    error("Total reads:\t%lld", reads);
    error("Mapped reads:\t%lld", mapped_reads);
    error("Duplicate reads:\t%lld", duplicate_reads);
    hts_itr_destroy(iter);

    iter = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
    if (!stats_only) {
        for(reads = 0; sam_itr_next(in, iter, b) >= 0; reads++) {
            if (b->core.tid < 0) {
                /* Write unmapped reads as is */
                if (sam_write1(out, hdr, b) < 0) {
                    error("writing to standard output failed");
                    exit(1);
                }
                continue;
            }
            chr = b->core.tid;
            start = b->core.pos;
            stop = bam_endpos(b);
            len = stop - start;
            chr2 = 0;
            start2 = 0;
            if (b->core.mtid >= 0) {
                chr2 = b->core.mtid;
                start2 = b->core.mpos;
            }

            if (std::get<0>( mp[{PACK_CHRPOSLEN(chr, start, len), PACK_CHRPOSLEN(chr2, start2, len)}] ) == reads) {
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

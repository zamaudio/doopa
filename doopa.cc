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
#include <omp.h>

#include "htslib/thread_pool.h"
#include "htslib/hfile.h"
#include "htslib/sam.h"

typedef uint64_t chrposlen_t;

uint64_t key_hash(const chrposlen_t& k) {
    return k;
}

bool key_equal(const chrposlen_t& v1, const chrposlen_t& v2) {
    return (v1 == v2);
}

typedef std::unordered_map<chrposlen_t, std::pair<uint64_t, uint64_t>,
        std::function<uint64_t(const chrposlen_t&)>,
        std::function<bool(const chrposlen_t&, const chrposlen_t&)> > doopa_t;

// Pack into 64 bits:
// chr  start   len
// ff8 7fffffff ffffff
#define PACK_CHRPOSLEN(chr, pos, len) \
    (uint64_t)( (((uint64_t)(chr) & 0x1ff) << 55) | \
                (((uint64_t)(pos) & 0x7fffffff) << 24) | \
                 ((uint64_t)(len) & 0xffffff) \
              )

static inline uint64_t get_qualsum(const bam1_t *b)
{
    int i;
    uint64_t q;
    uint8_t *qual = bam_get_qual(b);
    for (i = q = 0; i < b->core.l_qseq; ++i) q += qual[i];
    return q;
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

static void dedup_bam(const char *filename)
{
    htsThreadPool p = {NULL, 0};
    uint64_t reads = 0;
    uint64_t qualsum, existing_qual;
    int32_t chr, start, stop, len;
    hts_itr_t *iter;
    bam1_t *b;
    bam_hdr_t *hdr = NULL;
    samFile *out = NULL;
    samFile *in = NULL;
    hts_idx_t *idx = NULL;
    htsFormat _bam;
    hts_parse_format(&_bam, "bam");

    in = sam_open_format(filename, "r", (const htsFormat *)&_bam);
    if (!in) {
        error("cannot open bam file");
    }

    if ((idx = sam_index_load(in, filename)) == 0) {
        error("cannot open bam index");
        exit(1);
    }

    doopa_t mp((doopa_t::size_type)1000000, key_hash, key_equal);

    out = sam_open_format("/dev/stdout", "w", (const htsFormat *)&_bam);

    if (out == NULL) { error("reopening standard output failed"); goto clean; }

    if (!(p.pool = hts_tpool_init(8))) {
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

    if (sam_hdr_write(out, hdr) != 0) {
        error("writing headers to standard output failed");
        goto clean;
    }

    error("Start deduping...");

    iter = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
    b = bam_init1();
    if (b == NULL) { error("can't create record"); exit(1); }
    for(; sam_itr_next(in, iter, b) >= 0; reads++) {
        if (b->core.tid < 0) {
            continue;
        }
        chr = b->core.tid;
        start = b->core.pos;
        stop = bam_endpos(b);
        len = stop - start;
        qualsum = get_qualsum(b);
        if (mp.count(PACK_CHRPOSLEN(chr, start, len))) {
            // Key exists
            existing_qual = std::get<1>(mp[PACK_CHRPOSLEN(chr, start, len)]);
            if (qualsum > existing_qual) {
                mp[PACK_CHRPOSLEN(chr, start, len)] = std::make_pair(reads, qualsum);
            }
        } else {
            mp[PACK_CHRPOSLEN(chr, start, len)] = std::make_pair(reads, qualsum);
        }
    }
    error("Found %d mapped reads, deduped to %d reads, now writing output", reads, mp.size());
    hts_itr_destroy(iter);

    iter = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
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
        if (std::get<0>( mp[PACK_CHRPOSLEN(chr, start, len)] ) == reads) {
            if (sam_write1(out, hdr, b) < 0) {
                error("writing to standard output failed");
                exit(1);
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
    dedup_bam(argv[1]);
    return 0;
}

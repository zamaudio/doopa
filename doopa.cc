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
#include <unordered_map>
#include <tuple>
#include <utility>
#include <functional>

#include "htslib/hfile.h"
#include "htslib/sam.h"

struct cheater {
    int one, two, three, reads;
};

typedef std::tuple<int32_t, int32_t, int32_t> chrposlen_t;

// Pack into 64 bits:
// chr  start   len
// ff fffffffe 1ffffff
uint64_t key_hash(const chrposlen_t& k) {
    return ((((uint64_t)std::get<0>(k) & 0xff) << 56) |
            (((uint64_t)std::get<1>(k) & 0x7fffffff) << 25) |
             ((uint64_t)std::get<2>(k) & 0x1ffffff) );
}

bool key_equal(const chrposlen_t& v1, const chrposlen_t& v2) {
    /*
    return ( (std::get<0>(v1) == std::get<0>(v2)) &&
             (std::get<1>(v1) == std::get<1>(v2)) &&
             (std::get<2>(v1) == std::get<2>(v2)));
    */
    return (key_hash(v1) == key_hash(v2));
}

typedef std::unordered_map<chrposlen_t, std::string,
        std::function<uint64_t(const chrposlen_t&)>,
        std::function<bool(const chrposlen_t&, const chrposlen_t&)> > doopa_t;

void error(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    fflush(stdout);
    fprintf(stderr, "doopa: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    va_end(args);
}

static htsFile *dup_stdout(const char *mode)
{
    int fd = dup(STDOUT_FILENO);
    hFILE *hfp = (fd >= 0)? hdopen(fd, mode) : NULL;
    return hfp? hts_hopen(hfp, "-", mode) : NULL;
}

static void dedup_sam(samFile *in, const char *filename)
{
    int ret = 0;
    char region[128] = { 0 };
    int32_t chr, start, stop, len, reads = 0;
    struct cheater *getreads;

    bam1_t *b = NULL;
    bam_hdr_t *hdr = NULL;
    samFile *out = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;

    if ((idx = sam_index_load(in, filename)) == 0) {
        error("cannot open bam index");
        exit(1);
    }

    getreads = (struct cheater *)idx;
    doopa_t mp((doopa_t::size_type)getreads->reads, key_hash, key_equal);

    hdr = sam_hdr_read(in);
    if (hdr == NULL) {
        errno = 0; error("reading headers from \"%s\" failed", filename);
        goto clean;
    }

    out = dup_stdout("w");
    if (out == NULL) { error("reopening standard output failed"); goto clean; }

    if (sam_hdr_write(out, hdr) != 0) {
        error("writing headers to standard output failed");
        goto clean;
    }

    b = bam_init1();
    if (b == NULL) { error("can't create record"); goto clean; }

    if ((iter = sam_itr_queryi(idx, HTS_IDX_START, 0, 0)) == 0) {
        error("can't create iterator from start");
        goto clean;
    }

    error("Found %d reads, deduping...", getreads->reads);

    while (sam_itr_next(in, iter, b) >= 0) {
        if (b->core.tid < 0) {
            continue;
        }
        reads++;
        chr = b->core.tid;
        start = b->core.pos;
        stop = bam_endpos(b);
        len = stop - start;
        snprintf(region, 64, "%s:%d-%d", hdr->target_name[chr], start, stop);
        mp[{chr, start, len}] = std::string(region);
    }
    hts_itr_destroy(iter);

    error("Found %d mapped reads, deduped to %d reads, now writing output", reads, mp.size());

    for (doopa_t::iterator itr = mp.begin(); itr != mp.end(); ++itr) {
        iter = sam_itr_querys(idx, hdr, itr->second.c_str());
        ret = sam_itr_next(in, iter, b);
        if (sam_write1(out, hdr, b) < 0) {
            error("writing to standard output failed");
            goto clean;
        }
        hts_itr_destroy(iter);
    }

    error("Done");

 clean:
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    if (out) hts_close(out);
}

int main(int argc, char **argv)
{
    hFILE *fp;
    htsFile *hts;

    if (argc == 2) {
        fp = hopen(argv[1], "r");
        if (fp == NULL) {
            error("can't open \"%s\"", argv[1]);
            return 1;
        }
        hts = hts_hopen(fp, argv[1], "r");
    } else {
        error("need bam file path and bai to exist");
        return 1;
    }

    if (hts) {
        switch (hts_get_format(hts)->category) {
        case sequence_data:
            dedup_sam(hts, argv[1]);
            break;
        default:
            error("can't open bam file \"%s\": unknown format", argv[1]);
            break;
        }

        if (hts_close(hts) < 0) error("closing \"%s\" failed", argv[1]);
        fp = NULL;
    }
    
    if (fp && hclose(fp) < 0) error("closing \"%s\" failed", argv[1]);
    return 0;
}

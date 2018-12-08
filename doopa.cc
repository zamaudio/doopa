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

#include "htslib/hts.h"
#include "htslib/hfile.h"
#include "htslib/sam.h"

typedef std::tuple<int32_t, int32_t, int32_t> chrposlen_t;

struct key_hash : public std::unary_function<chrposlen_t, std::size_t>
{
   std::size_t operator()(const chrposlen_t& k) const
   {
      return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
   }
};

struct key_equal : public std::binary_function<chrposlen_t, chrposlen_t, bool>
{
   bool operator()(const chrposlen_t& v0, const chrposlen_t& v1) const
   {
      return (
               std::get<0>(v0) == std::get<0>(v1) &&
               std::get<1>(v0) == std::get<1>(v1) &&
               std::get<2>(v0) == std::get<2>(v1)
             );
   }
};

typedef std::unordered_map<chrposlen_t, int64_t, key_hash, key_equal> doopa_t;

void error(const char *format, ...)
{
    int err = errno;
    va_list args;
    va_start(args, format);
    fflush(stdout);
    fprintf(stderr, "doopa: ");
    vfprintf(stderr, format, args);
    if (err) fprintf(stderr, ": %s\n", strerror(err));
    else fprintf(stderr, "\n");
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
    doopa_t mp;
    
    bam1_t *b = NULL;
    bam_hdr_t *hdr = NULL;
    samFile *out = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;

    if ((idx = sam_index_load(in, filename)) == 0) {
            error("cannot open bam index");
        goto clean;
    }

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

    iter = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);

    while (1) {
        ret = iter? sam_itr_next(in, iter, b) : sam_read1(in, hdr, b);
        if (ret < 0) {
                break;
        }
        mp[std::make_tuple(b->core.tid, b->core.pos, bam_endpos(b))] = iter->tid;
    }

    for (doopa_t::iterator itr = mp.begin(); itr != mp.end(); ++itr) {
        iter = sam_itr_queryi(idx, itr->second, std::get<1>(itr->first), std::get<2>(itr->first));
        ret = iter? sam_itr_next(in, iter, b) : sam_read1(in, hdr, b);
        if (ret < 0) {
            error("can't find read in bam");
            goto clean;
        }
        if (sam_write1(out, hdr, b) < 0) {
            error("writing to standard output failed");
            goto clean;
        }
    }

 clean:
    hts_itr_destroy(iter);
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

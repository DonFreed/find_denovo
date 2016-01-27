#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <inttypes.h>
#include "htslib/vcf.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"

KHASH_MAP_INIT_STR(ped, int)

/* 
 * Read a .ped file and return familial relationships as a map.
 */
khash_t(ped) *read_ped_file(const char *fnped, kstring_t *fam_ids)
{
    int ret, res, ln = 0, n = 0, *offsets = 0, max = 0, idx, i;
    khiter_t k;
    FILE *fp;
    char *s;
    kstring_t line = {0, 0, 0 };
    khash_t(ped) *h = kh_init(ped);
    
    fp = fopen(fnped, "r");
    while ((res = kgetline(&line, (kgets_func *)fgets, fp)) != EOF) {
        ln++;
        n = ksplit_core(line.s, '\t', &max, &offsets);
        if (n != 6) {
            fprintf(stderr, "Error: ped file is poorly formatted at line %d. Incorrect number of columns\n", ln);
            return 0;
        }
        if (line.s[offsets[2]] != '0') {
            line.s[offsets[3] - 1] = '\t';
            idx = fam_ids->l;
            res = kputsn(line.s + offsets[2], offsets[4] - offsets[2], fam_ids);
            //fprintf(stderr, "Adding %s to the hash table with key %u\n", line.s + offsets[1], kh_str_hash_func(line.s + offsets[1]));
            k = kh_put(ped, h, strdup(line.s + offsets[1]), &ret); // FIXME: strdup is less efficent than the allocation of a memory block; no free
            //if (ret == -1) fprintf(stderr, "Put failed\n");
            kh_value(h, k) = idx;
        }
        line.l = 0;
    }

    free(line.s);
    free(offsets);
    return h;
}

/*
inline void get_pa_gts(bcf_hdr_t *hdr, bf1_t *line, const char *pa_ids, int *fa_gt, int *mo_gt)
{
    int *offsets, n;
    kstring_t s = {0, 0, 0};
    kputs(pa_ids, &s);
    offsets = ksplit(&s, ';', &n);
}
*/

inline int fmt_to_int(uint8_t **_fmt, int type)
{
    uint8_t *fmt = *_fmt;
    int i, ret=0;
    switch (type) {
        case BCF_BT_INT8: ret = (int)*(fmt++); break;
        case BCF_BT_INT16: for (i = 0; i < 2; ++i) ret = (ret << 8) | (int)*(fmt++); break;
        case BCF_BT_INT32: for (i = 0; i < 4; ++i) ret = (ret << 8) | (int)*(fmt++); break;
        default: fprintf(stderr, "Error: unknown fmt type %d\n", type); abort(); break;
    }
    *_fmt = fmt;
    return ret;
}

int find_denovo(bcf_hdr_t *hdr, bcf1_t *line, int *dnv_vals, khash_t(ped) *h, kstring_t *fam_ids, int min_dp, int min_alt, int par_pl_pen, int prb_pl_pen)
{
    int i, j, res, n, max = 0, *offsets = 0, int_id, n_denovo = 0;
    uint8_t *gt, *pa_gts[2], *pa_gt_end[2];
    khiter_t k;
    kstring_t s = {0, 0, 0};
    bcf_fmt_t *gt_ptr = bcf_get_fmt(hdr, line, "GT"), *ad_ptr = bcf_get_fmt(hdr, line, "AD");
//    fprintf(stderr, "3\n");
    for (i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        int last_allele = -1, cur_allele, cur_pa[2] = {-1, -1}, denovo = -1;
        k = kh_get(ped, h, bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
        fprintf(stderr, "individual %d or %s\n", i, bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
        if ((k == kh_end(h))) {
            continue;
        }
//        fprintf(stderr, "5\n");
        gt = gt_ptr->p + (i * gt_ptr->size);
        res = kh_value(h, k);
        kputs(fam_ids->s + res, &s);
//        fprintf(stderr, "6\n");
        n = ksplit_core(s.s, '\t', &max, &offsets);
//        fprintf(stderr, "10\n");
        assert(n == 2);
        for (j = 0; j < 2; ++j) {
//            fprintf(stderr, "11\n");
            int_id = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, s.s + offsets[j]);
            fprintf(stderr, "parent %d column %d\n", j, int_id);
            pa_gts[j] = gt_ptr->p + (int_id * gt_ptr->size);
            pa_gt_end[j] = pa_gts[j] + gt_ptr->size;
        }
//        fprintf(stderr, "7\n");
        while (gt < gt_ptr->p + ((i + 1) * gt_ptr->size)) {
            cur_allele = gt_to_int(&gt, gt_ptr->type);
            if (bcf_gt_is_missing(cur_allele)) break;
            cur_allele = bcf_gt_allele(cur_allele);
            fprintf(stderr, "Current allele %d\n", cur_allele);
            if (cur_allele == last_allele) continue;
            for (j = 0; j < 2; ++j) {
                while (cur_allele > cur_pa[j] && pa_gts[j] < pa_gt_end[j]) {
                    cur_pa[j] = bcf_gt_allele(gt_to_int(&pa_gts[j], gt_ptr->type));
                }
            }
            if (cur_allele != cur_pa[0] && cur_allele != cur_pa[1]) {
                denovo = cur_allele;
            }
            last_allele = cur_allele;
        }
//        fprintf(stderr, "8\n");
        if (denovo >= 0) {
            fprintf(stderr, "Found de novo allele %d\n", denovo);
            ++n_denovo;
            dnv_vals[i] = denovo + 1;
        }
//        fprintf(stderr, "9\n");
        s.l = 0;
        s.s[0] = '\0';
    }
    free(s.s);
    free(offsets);
    return n_denovo;
}

int main(int argc, char *argv[])
{
    int c, min_dp = 20, min_alt = 3, par_pl_pen = 20, prb_pl_pen = 20, help = 0, compression = 7, i, res, *dnv_vals = 0, found_dnv, n_samples;
    char *fnin = 0, *fnout = 0, *fndenovo = 0, *fnped = 0;
    char out_mode = 'v';
    kstring_t fam_ids = {0, 0, 0};
    khash_t(ped) *h;
    khiter_t k;
    htsFile *fin;
    bcf_hdr_t *hdr;
    bcf1_t *line;
    while ((c = getopt(argc, argv, "c:a:s:l:i:p:o:d:O:h")) >= 0) {
        if (c == 'c') min_dp = atoi(optarg);
        else if (c == 'a') min_alt = atoi(optarg);
        else if (c == 's') denovo_penalty = atoi(optarg);
        else if (c == 'l') compression = atoi(optarg);
        else if (c == 'i') fnin = strdup(optarg);
        else if (c == 'p') fnped = strdup(optarg);
        else if (c == 'o') fnout = strdup(optarg);
        else if (c == 'd') fndenovo = strdup(optarg);
        else if (c == 'O') out_mode = optarg[0];
        else if (c == 'h') help = 1;
        else break;
    }
    if (argc > optind || help) {
        fprintf(stderr, "\n");
        fprintf(stderr, "About:   Identify de novo variants in a VCF/BCF file using Phred-scaled genotype likelihoods.\n");
        fprintf(stderr, "Usage:   mark_denovo [options] [-O <v|z|b|u>] [-i <in.vcf|in.vcf.gz|in.bcf>] -p in.ped [-o <out.vcf|out.vcf.gz|out.bcf>]\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "           -c INT          Minimum number of reads in all trio members [%d]\n", min_dp);
        fprintf(stderr, "           -a INT          Minimum number of reads supporting the alternate allele [%d]\n", min_alt);
        fprintf(stderr, "           -s INT          Minimum phred-scaled confidence for the parental genotype [%d]\n", par_pl_pen);
        fprintf(stderr, "           -s INT          Minimum phred-scaled confidence for the child's genotype [%d]\n", prb_pl_pen);
        fprintf(stderr, "           -l INT          zlib compression level for the output VCF/BCF [%d]\n", compression);
        fprintf(stderr, "           -i FILE         The input file [stdin]\n");
        fprintf(stderr, "           -p FILE         The input PED file\n");
        fprintf(stderr, "           -o FILE         The output file [stdout]\n");
        fprintf(stderr, "           -d FILE         Output abbreviated information on the identified de novo variants (recommended for large inputs)\n");
        fprintf(stderr, "           -O <v|z|b|u>    v: VCF, z: bgzip compressed VCF, b: BCF, u: uncompressed BCF [v]\n");
        fprintf(stderr, "           -h              Print this help information\n");
        fprintf(stderr, "\n");
        return -1;
    }
    if (!fnped) {
        fprintf(stderr, "Error: Please specify a .ped file with the -p option\n");
        return -1;
    }
    if (!fnin) fnin = "-";
    if (!fnout) fnout = "-";
    if (!(out_mode == 'v' || out_mode == 'b' || out_mode == 'z' || out_mode == 'u')) {
        fprintf(stderr, "Error: Output file mode must be v, b, z, or u\n");
    }
    if (compression < 0) compression = 0;
    if (compression > 9) compression = 9;
    
//    fprintf(stderr, "-1\n");
    if (!(h = read_ped_file(fnped, &fam_ids))) return -1;

//    fprintf(stderr, "0\n");
    fin = bcf_open(fnin, "r");
    hdr = bcf_hdr_read(fin);
    bcf_hdr_append(hdr, "##FORMAT=<ID=DN,Number=1,Type=Integer,Description=\"Allele is de novo in sample, otherwise -1\">");
//    fprintf(stderr, "1\n");
    n_samples = bcf_hdr_nsamples(hdr);
    if (!(dnv_vals = calloc(n_samples, sizeof(int)))) {
        fprintf(stderr, "Error: Memory allocation failure\n");
        return -1;
    }
//    fprintf(stderr, "2\n");
    line = bcf_init();
//    fprintf(stderr, "Reading VCF/BCF\n");
    while ((bcf_read(fin, hdr, line) != -1)) {
        fprintf(stderr, "Reading line chr%" PRId32 ":%" PRId32 "\n", line->rid, line->pos);
        bcf_unpack(line, BCF_UN_ALL);
        fprintf(stderr, "Alleles found: ");
        for (i = 0; i < line->n_allele; ++i) {
            fprintf(stderr, "%s,", line->d.allele[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "finding denovo\n");
        found_dnv = find_denovo(hdr, line, dnv_vals, h, &fam_ids, min_dp, min_alt, par_pl_pen, prb_pl_pen);
        if (found_dnv) {
            printf("Found %d de novo variants at chromsome %" PRId32 ":%" PRId32 " in samples: ", found_dnv, line->rid, line->pos);
            for (i = 0; i < n_samples; ++i) {
                if (dnv_vals[i]) printf("%d,", i);
            }
            printf("\n");
        }
    }

    /*
    // Test the ped file //
    {
        const char *test[5];
        int is_missing;
        test[0] = "11000.p1";
        test[1] = "11007.p1";
        test[2] = "11007.s1";
        test[3] = "11008.fa";
        test[4] = "14698.s1";
        for (i = 0; i < 5; ++i) {
            k = kh_get(ped, h, test[i]);
            fprintf(stderr, "Key %s\n", kh_exist(h, k)? "exists" : "does not exist");
            if ((is_missing = (k == kh_end(h)))) {
                fprintf(stderr, "Key %s not found in hash table with key %u\n", test[i], kh_str_hash_func(test[i]));
            } else {
                res = kh_value(h, k);
                fprintf(stderr, "Key %s found with value %d corresponding to %s\n", test[i], res, fam_ids.s + res);
            }
        }
    }
    */

    /*
    // Print parameters //
    printf("min_dp - %d\n", min_dp);
    printf("min_alt - %d\n", min_alt);
    printf("denovo_penalty - %d\n", denovo_penalty);
    printf("compression - %d\n", compression);
    if (fnin) { printf("fnin - %s\n", fnin); }
    if (fnout) { printf("fnout - %s\n", fnout); }
    if (fndenovo) { printf("fndenovo - %s\n", fndenovo); }
    if (out_mode) { printf("out_mode - %c\n", out_mode); }
    if (fnped) { printf("fnped - %s\n", fnped); }
    */
//    for (k = kh_begin(h); k != kh_end(h); ++k)
//        if (kh_exist(h, k)) free(&kh_key(h, k));
    free(fam_ids.s);
    kh_destroy(ped, h);
    return 0;
}

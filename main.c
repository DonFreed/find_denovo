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
    int ret, res, ln = 0, n = 0, *offsets = 0, max = 0, idx;
    khiter_t k;
    FILE *fp;
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

inline int fmt_to_int(uint8_t **_fmt, int type)
{
    uint8_t *fmt = *_fmt;
    int i, ret=0;
    switch (type) {
        case BCF_BT_INT8: ret = (int)*(fmt++); break;
        case BCF_BT_INT16: for (i = 0; i < 2; ++i) ret |= ((int)*(fmt++) << (8 * i)); break;
        case BCF_BT_INT32: for (i = 0; i < 4; ++i) ret |= ((int)*(fmt++) << (8 * i)); break;
        default: fprintf(stderr, "Error: unknown fmt type %d\n", type); abort(); break;
    }
    *_fmt = fmt;
    return ret;
}

inline int find_ploidy(int size, int type)
{
    switch (type) {
        case BCF_BT_INT8: return size / 1;
        case BCF_BT_INT16: return size / 2;
        case BCF_BT_INT32: return size / 4;
        default: fprintf(stderr, "Error: unknown fmt type %d\n", type); abort(); break;
    }
}

inline int n_pl(n_allele, ploidy)
{
    float res = 1.0;
    int i;
    for (i = 0; i < ploidy; ++i) {
        res *= ((float)(n_allele + i)) / (i + 1);
    }
    return (int)(res + 0.9999);
}

inline void unwind_stk(int *stk, int *idx)
{
    for (; *idx > 0; --(*idx)) {
        if (stk[*idx] < stk[*idx - 1]) {
            stk[*idx]++;
            return;
        } else {
            stk[*idx] = 0;
        }
    }
    *idx = 0;
    stk[0]++;
    return;
}

inline void find_allele_pl2(int allele, int n_allele, int ploidy, int **_allele_pl, int *_max_pl)
{
    int *stk, pl_pos = 0, idx = 0, i, out_idx = 0, *allele_pl = *_allele_pl, max_pl = *_max_pl;
    stk = calloc(ploidy, sizeof(int));
    while (stk[0] < n_allele) {
        if (stk[idx] == allele) {
            for (i = 0; i < n_pl(stk[idx] + 1, ploidy - idx - 1); ++i) {
                if (out_idx >= max_pl) {
                    int *tmp;
                    max_pl = max_pl? max_pl << 1 : 32;
                    if ((tmp = (int*)realloc(allele_pl, sizeof(int) * max_pl))) {
                        allele_pl = tmp;
                    } else {
                        free(allele_pl);
                        allele_pl = 0;
                        return;
                    }
                }
                allele_pl[out_idx++] = pl_pos++;
            }
            unwind_stk(stk, &idx);
        } else if (idx + 1 < ploidy) {
            stk[++idx] = 0;
        } else {
            pl_pos++;
            unwind_stk(stk, &idx);
        }
    }
    *_allele_pl = allele_pl;
    *_max_pl = max_pl;
    return;
}

/*
inline void find_allele_pl(int denovo, int n_allele, int ploidy, int **_allele_pl, int *max_pl)
{
    int *gt, idx = ploidy - 1, pos = 0, sum = 0, i = 0, j = 0, allele_found, *allele_pl = *_allele_pl, n_out = 0;
    gt = calloc(ploidy, sizeof(int));
    n_allele += 1; // Ensure array encompasses end + 1
    while (sum <= ploidy * (n_allele - 1)) {
        allele_found = 0;
        for (i = 0; i < ploidy; ++i) {
            if (gt[0] == denovo) allele_found = 1;
        }
        if (allele_found) {
            if (n_out >= *max_pl) {
                int *tmp;
                *max_pl = *max_pl? *max_pl<<1 : 2;
                if ((tmp = (int*)realloc(allele_pl, sizeof(int) * *max_pl))) {
                    allele_pl = tmp;
                } else {
                    free(allele_pl);
                    allele_pl = 0;
                    return;
                }
            }
            allele_pl[n_out++] = pos;
        }

        if (idx == 0) {
            while (1) {
                if (idx + 1 >= ploidy) {
                    int tmp = gt[idx] + 1;
                    memset(gt, 0, ploidy * sizeof(int));
                    gt[idx] = tmp;
                    idx -= 1;
                    pos += 1;
                    break;
                }
                if (gt[idx + 1] > gt[idx]) {
                    gt[idx] += 1;
                    for (j = 0; j < idx; ++j) gt[j] = 0;
                    idx = idx? idx - 1: 0;
                    pos += 1;
                    break;
                }
                idx += 1;
            }
        } else {
            gt[idx] += 1;
            idx -= 1;
            pos += 1;
        }
        sum = 0;
        for (i = 0; i < ploidy; ++i) sum += gt[i];
    }
    *_allele_pl = allele_pl;
    return;
}
*/

/*
            allele[0] = fmt_to_int(&gt[0], gt_ptr->type);
            //if (bcf_gt_is_missing(allele[0])) break; // missing = 0
            allele[0] = bcf_gt_allele(allele[0]);
*/

inline void count_allele_indivs(bcf_hdr_t *hdr, bcf1_t *line, int **_alleles, int *_max_alleles)
{
    int *alleles = *_alleles, max_alleles = *_max_alleles, tmp_allele, i, last_allele;
    bcf_fmt_t *gt_ptr = bcf_get_fmt(hdr, line, "GT");
    uint8_t *gt = gt_ptr->p, *gt_end;
    for (i = 0; i < line->n_sample; i++) {
        gt_end = gt + gt_ptr->size;
        last_allele = -1;
        while(gt < gt_end) {
        /* Write summary information */

            tmp_allele = fmt_to_int(&gt, gt_ptr->type);
            tmp_allele = bcf_gt_allele(tmp_allele);
            if (tmp_allele != last_allele) {
                if (tmp_allele >= max_alleles) {
                    int *tmp, old = max_alleles;
                    max_alleles = tmp_allele;
                    kroundup32(max_alleles);
                    max_alleles = max_alleles < 32? 32 : max_alleles << 1;
                    if ((tmp = (int*)realloc(alleles, sizeof(int) * max_alleles))) {
                        alleles = tmp;
                    } else {
                        free(alleles);
                        alleles = 0;
                        return;
                    }
                    memset(alleles + old, 0, sizeof(int) * (max_alleles - old));
                }
                alleles[tmp_allele]++;
                last_allele = tmp_allele;
            }
        }
    }

    *_alleles = alleles;
    *_max_alleles = max_alleles;
    return;
}

int find_denovo(bcf_hdr_t *hdr, bcf1_t *line, int *dnv_vals, khash_t(ped) *h, kstring_t *fam_ids, int min_dp, int min_alt, int par_pl_pen, int prb_pl_pen)
{
    int i, j, res, n, max = 0, *offsets = 0, int_id, n_denovo = 0, ploidy, *allele_pl = 0, max_pl = 0;
    int pl_idx, dnv_idx, min_pl;
    uint8_t *gt[3], *gt_end[3], *ad[3], *ad_end[3], *pl[3], *pl_end[3];
    khiter_t k;
    kstring_t s = {0, 0, 0};
    bcf_fmt_t *gt_ptr = bcf_get_fmt(hdr, line, "GT"), *ad_ptr = bcf_get_fmt(hdr, line, "AD"), *pl_ptr = bcf_get_fmt(hdr, line, "PL");
    if ((!(gt_ptr)) || (!(ad_ptr)) || (!(pl_ptr))) return 0;
    ploidy = find_ploidy(gt_ptr->size, gt_ptr->type);
    for (i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        int last_allele = -1, allele[3] = {-1, -1, -1}, denovo = -1, depth[3] = {0, 0, 0}, pass_filt = 1;
        k = kh_get(ped, h, bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
//        fprintf(stderr, "individual %d or %s\n", i, bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
        if ((k == kh_end(h))) {
            continue;
        }

        /* Get child data */
        gt[0] = gt_ptr->p + (i * gt_ptr->size);
        gt_end[0] = gt[0] + gt_ptr->size;
        ad[0] = ad_ptr->p + (i * ad_ptr->size);
        ad_end[0] = ad[0] + ad_ptr->size;
        pl[0] = pl_ptr->p + (i * pl_ptr->size);
        pl_end[0] = pl[0] + pl_ptr->size;

        /* Get parental data */
        s.l = 0;
        res = kh_value(h, k);
        kputs(fam_ids->s + res, &s);
        n = ksplit_core(s.s, '\t', &max, &offsets);
        assert(n == 2);
        for (j = 1; j < 3; ++j) {
            int_id = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, s.s + offsets[j - 1]);
//            fprintf(stderr, "parent %d column %d\n", j - 1, int_id);
            gt[j] = gt_ptr->p + (int_id * gt_ptr->size);
            gt_end[j] = gt[j] + gt_ptr->size;
            ad[j] = ad_ptr->p + (int_id * ad_ptr->size);
            ad_end[j] = ad[j] + ad_ptr->size;
            pl[j] = pl_ptr->p + (int_id * pl_ptr->size);
            pl_end[j] = pl[j] + pl_ptr->size;
        }

        /* Check genotypes */
//        fprintf(stderr, "Child's alleles: ");
        while (gt[0] < gt_end[0]) {
            allele[0] = fmt_to_int(&gt[0], gt_ptr->type);
            //if (bcf_gt_is_missing(allele[0])) break; // missing = 0
            allele[0] = bcf_gt_allele(allele[0]);
//            fprintf(stderr, "%d,", allele[0]);
            if (allele[0] == last_allele) continue;
            for (j = 1; j < 3; ++j) {
                while (allele[0] > allele[j] && gt[j] < gt_end[j]) {
                    allele[j] = bcf_gt_allele(fmt_to_int(&gt[j], gt_ptr->type));
                }
            }
            if (allele[0] != allele[1] && allele[0] != allele[2]) {
                denovo = allele[0];
            }
            last_allele = allele[0];
        }
//        fprintf(stderr, "\n");
        if (denovo < 0) continue;

        fprintf(stderr, "individual %s, column %i has denovo allele\n", bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i), i);

        /* Check depths */
//        fprintf(stderr, "Child depths: ");
        j = 0;
        while (ad[0] < ad_end[0]) {
            int tmp_ad = fmt_to_int(&ad[0], ad_ptr->type);
//            fprintf(stderr, "%d,", tmp_ad);
            /* Check alt depth */
            if (j == denovo) {
                if (tmp_ad < min_alt) { pass_filt = 0; }
            }
            depth[0] += tmp_ad;
            ++j;
        }
//        fprintf(stderr, "\n");
        if (!pass_filt) { fprintf(stderr, "Removed due to child alt depth\n"); continue; }
        /* Check total depth */
        for (j = 1; j < 3; ++j) {
            while (ad[j] < ad_end[j]) {
                depth[j] += fmt_to_int(&ad[j], ad_ptr->type);
            }
        }
        for (j = 0; j < 3; ++j) {
            if (depth[j] < min_dp) pass_filt = 0;
        }
        if (!pass_filt) { fprintf(stderr, "Continue due to family depth %d, %d, %d\n", depth[0], depth[1], depth[2]); continue; }

        /* Check proband PL tags */
        fprintf(stderr, "Calling find allele with %d, %d, %d\n", denovo, line->n_allele, ploidy);
        find_allele_pl2(denovo, line->n_allele, ploidy, &allele_pl, &max_pl);
        fprintf(stderr, "Child PL tags: ");
        pl_idx = 0; dnv_idx = 0; min_pl = 100000;
        while (pl[0] < pl_end[0]) {
            int tmp_pl = fmt_to_int(&pl[0], pl_ptr->type);
            fprintf(stderr, "%d,", tmp_pl);
            if (pl_idx == allele_pl[dnv_idx]) {
                pl_idx++;
                dnv_idx++;
                continue;
            }
            min_pl = tmp_pl < min_pl? tmp_pl : min_pl;
            pl_idx++;
        }
        fprintf(stderr, "\n");
        if (min_pl < prb_pl_pen) { 
            fprintf(stderr, "Continue due to child min pl %d\n", min_pl); 
            for (j = 0; j < max_pl - 1; ++j) { fprintf(stderr, "%d,", allele_pl[j]); }
            fprintf(stderr, "\n");
            continue; 
        }

        /* Check parental PL tags */
        for (j = 1; j < 3; ++j) {
            pl_idx = 0; dnv_idx = 0; min_pl = 100000;
            fprintf(stderr, "Parent %d PL tags: ", j);
            while (pl[j] < pl_end[j]) {
                int tmp_pl = fmt_to_int(&pl[j], pl_ptr->type);
                fprintf(stderr, "%d,", tmp_pl);
                if (pl_idx == allele_pl[dnv_idx]) {
                    min_pl = tmp_pl < min_pl? tmp_pl : min_pl;
                    dnv_idx++;
                }
                pl_idx++;
            }
            fprintf(stderr, "\n");
            if (min_pl < par_pl_pen) { fprintf(stderr, "min_pl %d < par_pl pen %d\n", min_pl, par_pl_pen); pass_filt = 0; }
        }
        if (!pass_filt) { fprintf(stderr, "Continue due to parental min pl\n"); continue; }

        /* Found a de novo variant */
        if (denovo >= 0) {
//            fprintf(stderr, "Found de novo allele %d\n", denovo);
            ++n_denovo;
            dnv_vals[i] = denovo + 1; // 0 is not de novo
        }
    }
    free(allele_pl);
    free(s.s);
    free(offsets);
    return n_denovo;
}

int main(int argc, char *argv[])
{
    int c, min_dp = 20, min_alt = 3, par_pl_pen = 20, prb_pl_pen = 20, help = 0, compression = 7, i, j, res, *dnv_vals = 0, found_dnv, n_samples;
    int *alleles = 0, max_alleles = 0;
    uint32_t max_indiv = 2;
    char *fnin = 0, *fnout = 0, *fndenovo = 0, *fnped = 0;
    char _out_mode = 'v', out_mode[8] = "w";
    kstring_t fam_ids = {0, 0, 0}, alt_out = {0, 0, 0};
    khash_t(ped) *h;
    khiter_t k;
    htsFile *fin, *fout;
    FILE *fdenovo = 0;
    bcf_hdr_t *hdr;
    bcf1_t *line;
    while ((c = getopt(argc, argv, "c:a:s:t:l:i:p:o:d:O:n:h")) >= 0) {
        if (c == 'c') min_dp = atoi(optarg);
        else if (c == 'a') min_alt = atoi(optarg);
        else if (c == 's') prb_pl_pen = atoi(optarg);
        else if (c == 't') par_pl_pen = atoi(optarg);
        else if (c == 'l') compression = atoi(optarg);
        else if (c == 'i') fnin = strdup(optarg);
        else if (c == 'p') fnped = strdup(optarg);
        else if (c == 'o') fnout = strdup(optarg);
        else if (c == 'd') fndenovo = strdup(optarg);
        else if (c == 'O') _out_mode = optarg[0];
        else if (c == 'n') max_indiv = atoi(optarg);
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
        fprintf(stderr, "           -t INT          Minimum phred-scaled confidence for the parental genotype [%d]\n", par_pl_pen);
        fprintf(stderr, "           -s INT          Minimum phred-scaled confidence for the child's genotype [%d]\n", prb_pl_pen);
        fprintf(stderr, "           -l INT          zlib compression level for the output VCF/BCF [%d]\n", compression);
        fprintf(stderr, "           -i FILE         The input file [stdin]\n");
        fprintf(stderr, "           -p FILE         The input PED file\n");
        fprintf(stderr, "           -o FILE         The output file [stdout]\n");
        fprintf(stderr, "           -d FILE         Output abbreviated information on the identified de novo variants (recommended for large inputs)\n");
        fprintf(stderr, "           -O <v|z|b|u>    v: VCF, z: bgzip compressed VCF, b: BCF, u: uncompressed BCF [v]\n");
        fprintf(stderr, "           -n INT          The maximum number of non-sibling individuals with the allele for a de novo call [%d]\n", max_indiv);
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
    if (!(_out_mode == 'v' || _out_mode == 'b' || _out_mode == 'z' || _out_mode == 'u')) {
        fprintf(stderr, "Error: Output file mode must be v, b, z, or u\n");
    }
    if (fndenovo) fdenovo = fopen(fndenovo, "w");
    if (compression < 0) compression = 0;
    if (compression > 9) compression = 9;
    
    if (!(h = read_ped_file(fnped, &fam_ids))) return -1;

    fin = bcf_open(fnin, "r");
    out_mode[1] = _out_mode;
    sprintf(out_mode + 2, "%d", compression);
    fout = bcf_open(fnout, out_mode);
    if (!(hdr = bcf_hdr_read(fin))) {
        fprintf(stderr, "Error: Could not read header from %s\n", fnin);
        return -1;
    }
    bcf_hdr_append(hdr, "##FORMAT=<ID=DN,Number=1,Type=Integer,Description=\"Allele + 1 is de novo in sample, otherwise 0\">");
    if ((bcf_hdr_write(fout, hdr))) {
        fprintf(stderr, "Error: Could not write to file %s\n", fnout);
        return -1;
    }
    n_samples = bcf_hdr_nsamples(hdr);
    if (!(dnv_vals = calloc(n_samples, sizeof(int)))) {
        fprintf(stderr, "Error: Memory allocation failure\n");
        return -1;
    }
    line = bcf_init();
    while ((bcf_read(fin, hdr, line) != -1)) {
        fprintf(stderr, "Reading line chr%" PRId32 ":%" PRId32 "\n", line->rid, line->pos);

        /* Check filter */
        bcf_unpack(line, BCF_UN_FLT);
        if (!(bcf_has_filter(hdr, line, "PASS"))) {
            /* Write to output */
            bcf_write(fout, hdr, line);
            continue;
        }

        bcf_unpack(line, BCF_UN_ALL);

        /* Find putative de novo variants */
        memset(dnv_vals, 0, sizeof(int) * n_samples);
        found_dnv = find_denovo(hdr, line, dnv_vals, h, &fam_ids, min_dp, min_alt, par_pl_pen, prb_pl_pen);
        if (!found_dnv) {
            continue;
        }
        fprintf(stderr, "Found %d de novo variants at chromsome %" PRId32 ":%" PRId32 " in samples: ", found_dnv, line->rid, line->pos);
        for (i = 0; i < n_samples; ++i) {
            if (dnv_vals[i]) fprintf(stderr, "%d,", i);
        }
        fprintf(stderr, "\n");

        /* Remove multi-sample variants */
        if (max_alleles) memset(alleles, 0, sizeof(int) * (max_alleles - 1));
        if (found_dnv && max_indiv < line->n_sample) {
            count_allele_indivs(hdr, line, &alleles, &max_alleles);
            fprintf(stderr, "Found individuals with alleles: ");
            for (i = 0; i < max_alleles; ++i) {
                fprintf(stderr, "%d,", alleles[i]);
            }
            fprintf(stderr, "\n");

            for (i = 0; i < n_samples; ++i) {
                if (dnv_vals[i]) {
                    if (alleles[dnv_vals[i] - 1] > max_indiv) { 
                        fprintf(stderr, "Removing de novo in column %d as %d alleles were found\n", i, alleles[dnv_vals[i] - 1]);
                        dnv_vals[i] = 0;
                        found_dnv--;
                    }
                }
            }
        }

        /* Add de novo information to the format field */
        if (found_dnv) {
            int first = 1, n_vals = 0, res;
            void *vals = 0;
            bcf_update_format_int32(hdr, line, "DN", dnv_vals, n_samples);

            /* Write summary information */
            // Chrom, pos, ref, alt, info, fmt, sampleid, -, sample field
            if (fdenovo) {
                for (i = 1; i < line->n_allele; ++i) {
                    if (alt_out.l) kputsn(",", 1, &alt_out);
                    kputs(line->d.allele[i], &alt_out);
                }
                fprintf(fdenovo, "%s\t%" PRId32"\t%s\t%s\t", bcf_hdr_id2name(hdr, line->rid), line->pos + 1, line->d.allele[0], alt_out.s);

                /* Print info field */
                for (i = 0; i < line->n_info; ++i) {
                    const char *key = hdr->id[BCF_DT_ID][line->d.info[i].key].key;
                    void *p;
                    if (!first) fprintf(fdenovo, ";");
                    fprintf(fdenovo, "%s", key);
                    switch (line->d.info[i].type) {
                        case (BCF_BT_INT8 || BCF_BT_INT16 || BCF_BT_INT32):
                            res = bcf_get_info_values(hdr, line, key, &vals, &n_vals, BCF_HT_INT);
                            p = vals;
                            for (j = 0; j < res; ++j) {
                                uint32_t tmp = *((uint32_t *)p);
                                if (!j) fprintf(fdenovo, ",");
                                fprintf(fdenovo, "%" PRId32, tmp);
                                p += 4;
                            }
                            break;
                        case BCF_BT_FLOAT:
                            res = bcf_get_info_values(hdr, line, key, &vals, &n_vals, BCF_HT_REAL);
                            p = vals;
                            for (j = 0; j < res; ++j) {
                                float tmp = *((float *)p);
                                if (!j) fprintf(fdenovo, ",");
                                fprintf(fdenovo, "%f", tmp);
                                p += sizeof(float);
                            }
                            break;
                        case BCF_BT_CHAR:
                            res = bcf_get_info_values(hdr, line, key, &vals, &n_vals, BCF_HT_STR);
                            fprintf(fdenovo, "%s", ((char *)vals));
                            break;
                        default:
                            break;
                    }
                    first = 0;
                }

                /* Print format field */
                fprintf(fdenovo, "\tID:GT:AD:PL:DN");
                for (i = 0; i < n_samples; ++i) {
                    if (dnv_vals[i]) {
                        bcf_fmt_t *fmt_ptr = bcf_get_fmt(hdr, line, "GT");
                        uint8_t *fmt, *fmt_end;
                        int val;
                        const char *fmt_str[3];
                        fmt_str[0] = "AD"; fmt_str[1] = "PL"; fmt_str[2] = "DN";
                        fprintf(fdenovo, "\t%s:", bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
                        fmt = fmt_ptr->p + (i * fmt_ptr->size);
                        fmt_end = fmt + fmt_ptr->size;
                        first = 1;
                        while (fmt < fmt_end) {
                            val = bcf_gt_allele(fmt_to_int(&fmt, fmt_ptr->type));
                            if (!first) fprintf(fdenovo, "/");
                            fprintf(fdenovo, "%d", val);
                            first = 0;
                        }

                        for (j = 0; j < 3; ++j) {
                            fmt_ptr = bcf_get_fmt(hdr, line, fmt_str[j]);
                            fmt = fmt_ptr->p + (i * fmt_ptr->size);
                            fmt_end = fmt + fmt_ptr->size;
                            first = 1;
                            fprintf(fdenovo, ":");
                            while (fmt < fmt_end) {
                                val = fmt_to_int(&fmt, fmt_ptr->type);
                                if (!first) fprintf(fdenovo, ",");
                                fprintf(fdenovo, "%d", val);
                                first = 0;
                            }
                        }
                    }
                }
            }
            free(vals);
        }

        /* Write to the output VCF */
        bcf_write(fout, hdr, line);
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

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

typedef struct {
    int fa_id;
    int mo_id;
} parent_t;

KHASH_MAP_INIT_INT(ped2, parent_t)

typedef struct {
    int *stk;
    int *allele_pl;
    int max_allele;
} scratch_t;

inline void destroy_scratch(scratch_t *data)
{
    free(data->stk);
    free(data->allele_pl);
    free(data);
}

inline scratch_t *init_scratch() {
    scratch_t *tmp = calloc(1, sizeof(scratch_t));
    return tmp;
}

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
            k = kh_put(ped, h, strdup(line.s + offsets[1]), &ret);
            kh_value(h, k) = idx;
        }
        line.l = 0;
    }

    free(line.s);
    free(offsets);
    return h;
}

khash_t(ped2) *read_ped_file2(const char *fnped, bcf_hdr_t *hdr)
{
    int res, ln = 0, *offsets = 0, max = 0, n, ret = 0;
    FILE *fp;
    kstring_t line = {0, 0, 0};
    khash_t(ped2) *h = kh_init(ped2);
    khiter_t k;

    fp = fopen(fnped, "r");
    while ((res = kgetline(&line, (kgets_func *)fgets, fp)) != EOF) {
        ++ln;
        n = ksplit_core(line.s, '\t', &max, &offsets);
        if (n != 6) {
            fprintf(stderr, "Error: ped file is poorly formated at line %d. Incorrect number of columns\n", ln);
            return 0;
        }
        if (line.s[offsets[2]] != '0') {
            parent_t tmp;
            int child_col = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, line.s + offsets[1]);
            tmp.fa_id = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, line.s + offsets[2]);
            tmp.mo_id = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, line.s + offsets[3]);
            if (child_col >= 0 && tmp.fa_id >= 0 && tmp.mo_id >= 0) {
                k = kh_put(ped2, h, (khint32_t)child_col, &ret);
                kh_value(h, k) = tmp;
            } else {
                fprintf(stderr, "Warning: trio %s, %s, %s found in .ped file but not found in vcf header\n", line.s + offsets[1], line.s + offsets[2], line.s + offsets[3]);
            }
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

inline void find_allele_pl2(int allele, int n_allele, int ploidy, int **_allele_pl, int *_max_pl, int *stk)
{
    int pl_pos = 0, idx = 0, i, out_idx = 0, *allele_pl = *_allele_pl, max_pl = *_max_pl;
    memset(stk, 0, sizeof(int) * ploidy);
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
    allele_pl[out_idx] = -1;
    *_allele_pl = allele_pl;
    *_max_pl = max_pl;
    return;
}

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

int find_denovo(bcf_hdr_t *hdr, bcf1_t *line, int *dnv_vals, khash_t(ped2) *h, int min_dp, int min_alt, int par_pl_pen, int prb_pl_pen, scratch_t *scratch)
{
    int i, j, n_denovo = 0, ploidy;
    int pl_idx, dnv_idx, min_pl;
    uint8_t *gt[3], *gt_end[3], *ad[3], *ad_end[3], *pl[3], *pl_end[3];
    khiter_t k;
    parent_t parents;
    bcf_fmt_t *gt_ptr = bcf_get_fmt(hdr, line, "GT"), *ad_ptr = bcf_get_fmt(hdr, line, "AD"), *pl_ptr = bcf_get_fmt(hdr, line, "PL");
    if ((!(gt_ptr)) || (!(ad_ptr)) || (!(pl_ptr))) return 0;
    ploidy = find_ploidy(gt_ptr->size, gt_ptr->type);
    if (!(scratch->stk)) {
        int *tmp;
        if ((tmp = (int*)malloc(ploidy * sizeof(int)))) {
            scratch->stk = tmp;
        } else {
            return 0;
        }
    }
    for (i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        int last_allele = -1, allele[3] = {-1, -1, -1}, denovo = -1, depth[3] = {0, 0, 0}, pass_filt = 1;
        k = kh_get(ped2, h, i);
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

        /* Get father data */
        parents = kh_value(h, k);
        gt[1] = gt_ptr->p + (parents.fa_id * gt_ptr->size);
        gt_end[1] = gt[1] + gt_ptr->size;
        ad[1] = ad_ptr->p + (parents.fa_id * ad_ptr->size);
        ad_end[1] = ad[1] + ad_ptr->size;
        pl[1] = pl_ptr->p + (parents.fa_id * pl_ptr->size);
        pl_end[1] = pl[1] + pl_ptr->size;
        
        /* Get mother data */
        gt[2] = gt_ptr->p + (parents.mo_id * gt_ptr->size);
        gt_end[2] = gt[2] + gt_ptr->size;
        ad[2] = ad_ptr->p + (parents.mo_id * ad_ptr->size);
        ad_end[2] = ad[2] + ad_ptr->size;
        pl[2] = pl_ptr->p + (parents.mo_id * pl_ptr->size);
        pl_end[2] = pl[2] + pl_ptr->size;


        /* Check genotypes */
        while (gt[0] < gt_end[0]) {
            allele[0] = fmt_to_int(&gt[0], gt_ptr->type);
            allele[0] = bcf_gt_allele(allele[0]);
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
        if (denovo < 0) continue;

        /* Check depths */
        j = 0;
        while (ad[0] < ad_end[0]) {
            int tmp_ad = fmt_to_int(&ad[0], ad_ptr->type);
            /* Check alt depth */
            if (j == denovo) {
                if (tmp_ad < min_alt) { pass_filt = 0; }
            }
            depth[0] += tmp_ad;
            ++j;
        }
        if (!pass_filt) { continue; }

        /* Check total depth */
        for (j = 1; j < 3; ++j) {
            while (ad[j] < ad_end[j]) {
                depth[j] += fmt_to_int(&ad[j], ad_ptr->type);
            }
        }
        for (j = 0; j < 3; ++j) {
            if (depth[j] < min_dp) pass_filt = 0;
        }
        if (!pass_filt) { continue; }

        /* Check proband PL tags */
        find_allele_pl2(denovo, line->n_allele, ploidy, &(scratch->allele_pl), &(scratch->max_allele), scratch->stk);
        pl_idx = 0; dnv_idx = 0; min_pl = 100000;
        while (pl[0] < pl_end[0]) {
            int tmp_pl = fmt_to_int(&pl[0], pl_ptr->type);
            if (pl_idx == scratch->allele_pl[dnv_idx]) {
                pl_idx++;
                dnv_idx++;
                continue;
            }
            min_pl = tmp_pl < min_pl? tmp_pl : min_pl;
            pl_idx++;
        }
        if (min_pl < prb_pl_pen) { 
            continue; 
        }

        /* Check parental PL tags */
        for (j = 1; j < 3; ++j) {
            pl_idx = 0; dnv_idx = 0; min_pl = 100000;
            while (pl[j] < pl_end[j]) {
                int tmp_pl = fmt_to_int(&pl[j], pl_ptr->type);
                if (pl_idx == scratch->allele_pl[dnv_idx]) {
                    min_pl = tmp_pl < min_pl? tmp_pl : min_pl;
                    dnv_idx++;
                }
                pl_idx++;
            }
            if (min_pl < par_pl_pen) { pass_filt = 0; }
        }
        if (!pass_filt) { continue; }

        /* Found a de novo variant */
        if (denovo >= 0) {
            ++n_denovo;
            dnv_vals[i] = denovo + 1; // 0 is not de novo
        }
    }
    return n_denovo;
}

int main(int argc, char *argv[])
{
    int c, min_dp = 20, min_alt = 3, par_pl_pen = 20, prb_pl_pen = 20, help = 0, compression = 7, i, j, *dnv_vals = 0, found_dnv, n_samples;
    int *alleles = 0, max_alleles = 0, n_vals = 0, n_threads = 0;
    void *vals = 0;
    uint32_t max_indiv = 2;
    char *fnin = 0, *fnout = 0, *fndenovo = 0, *fnped = 0;
    char _out_mode = 'v', out_mode[8] = "w";
    kstring_t alt_out = {0, 0, 0};
    khash_t(ped2) *h;
    khiter_t k;
    htsFile *fin, *fout;
    FILE *fdenovo = 0;
    bcf_hdr_t *hdr;
    bcf1_t *line;
    scratch_t *scratch;
    while ((c = getopt(argc, argv, "c:a:s:t:l:i:p:o:d:O:n:x:h")) >= 0) {
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
        else if (c == 'x') n_threads = atoi(optarg);
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
        fprintf(stderr, "           -x INT          number of extra compression/decompression threads [%d]\n", n_threads);
        fprintf(stderr, "           -h              Print this help information\n");
        fprintf(stderr, "\n");
        return -1;
    }
    if (!fnped) {
        fprintf(stderr, "Error: Please specify a .ped file with the -p option\n");
        return -1;
    }
    if (!fnout) fnout = "-";
    if (!(_out_mode == 'v' || _out_mode == 'b' || _out_mode == 'z' || _out_mode == 'u')) {
        fprintf(stderr, "Error: Output file mode must be v, b, z, or u\n");
    }
    if (fndenovo) fdenovo = fopen(fndenovo, "w");
    if (compression < 0) compression = 0;
    if (compression > 9) compression = 9;
    
    if (fnin) {
        fin = bcf_open(fnin, "r");
    } else { 
        fin = bcf_open("-", "r");
    }
    if (n_threads) { hts_set_threads(fin, n_threads); }
    out_mode[1] = _out_mode;
    sprintf(out_mode + 2, "%d", compression);
    if (fnout) {
        fout = bcf_open(fnout, out_mode);
    } else {
        fout = bcf_open("-", out_mode);
    }
    if (n_threads) { hts_set_threads(fout, n_threads); }
    if (!(hdr = bcf_hdr_read(fin))) {
        fprintf(stderr, "Error: Could not read header from %s\n", fnin);
        return -1;
    }

    if (!(h = read_ped_file2(fnped, hdr))) return -1;

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
    scratch = init_scratch();
    line = bcf_init();
    while ((bcf_read(fin, hdr, line) != -1)) {

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
        found_dnv = find_denovo(hdr, line, dnv_vals, h, min_dp, min_alt, par_pl_pen, prb_pl_pen, scratch);
        if (!found_dnv) {
            bcf_write(fout, hdr, line);
            continue;
        }

        /* Remove multi-sample variants */
        if (max_alleles) memset(alleles, 0, sizeof(int) * (max_alleles - 1));
        if (found_dnv && max_indiv < line->n_sample) {
            count_allele_indivs(hdr, line, &alleles, &max_alleles);

            for (i = 0; i < n_samples; ++i) {
                if (dnv_vals[i]) {
                    if (alleles[dnv_vals[i] - 1] > max_indiv) { 
                        dnv_vals[i] = 0;
                        found_dnv--;
                    }
                }
            }
        }

        /* Add de novo information to the format field */
        if (found_dnv) {
            int first = 1, res;
            bcf_update_format_int32(hdr, line, "DN", dnv_vals, n_samples);

            /* Write summary information */
            // Chrom, pos, ref, alt, info, fmt, sampleid, -, sample field
            if (fdenovo) {
                alt_out.l = 0;
                for (i = 1; i < line->n_allele; ++i) {
                    if (alt_out.l) kputsn(",", 1, &alt_out);
                    kputs(line->d.allele[i], &alt_out);
                }
                fprintf(fdenovo, "%s\t%" PRId32"\t%s\t%s\t", bcf_hdr_id2name(hdr, line->rid), line->pos + 1, line->d.allele[0], alt_out.s);

                /* Print info field */
                for (i = 0; i < line->n_info; ++i) {
                    const char *key = hdr->id[BCF_DT_ID][line->d.info[i].key].key;
                    void *p;
                    if (i) fprintf(fdenovo, ";");
                    fprintf(fdenovo, "%s", key);
                    switch (line->d.info[i].type) {
                        case BCF_BT_INT8:
                        case BCF_BT_INT16:
                        case BCF_BT_INT32:
                            res = bcf_get_info_values(hdr, line, key, &vals, &n_vals, BCF_HT_INT);
                            p = vals;
                            for (j = 0; j < res; ++j) {
                                uint32_t tmp = *((uint32_t *)p);
                                if (j) {
                                    fprintf(fdenovo, ",");
                                } else {
                                    fprintf(fdenovo, "=");
                                }
                                fprintf(fdenovo, "%" PRId32, tmp);
                                p += 4;
                            }
                            break;
                        case BCF_BT_FLOAT:
                            res = bcf_get_info_values(hdr, line, key, &vals, &n_vals, BCF_HT_REAL);
                            p = vals;
                            for (j = 0; j < res; ++j) {
                                float tmp = *((float *)p);
                                if (j) {
                                    fprintf(fdenovo, ",");
                                } else {
                                    fprintf(fdenovo, "=");
                                }
                                fprintf(fdenovo, "%.3f", tmp);
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
                fprintf(fdenovo, "\n");
            }
        }

        /* Write to the output VCF */
        bcf_write(fout, hdr, line);
    }

    bcf_close(fin);
    bcf_close(fout);
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    destroy_scratch(scratch);
    free(vals);
    free(fnin);
    free(fnout);
    free(fndenovo);
    free(fnped);
    free(alt_out.s);
    free(alleles);
    free(dnv_vals);
    kh_destroy(ped2, h);
    return 0;
}

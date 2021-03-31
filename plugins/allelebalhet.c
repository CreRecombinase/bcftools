/*  plugins/allelebalhet.c -- allele balance on across sample hets

    Copyright (C) 2021 Genentech Inc.

    Authors: Zia Khan <khanz12@gene.com> and Nicholas Knoblauch <knoblaun@gene.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <htslib/kbitset.h>
#include <htslib/vcfutils.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include "bcftools.h"
#include <htslib/vcf.h>

typedef struct _args_t
{
  bcf_hdr_t *hdr;        /*! VCF file header */
  float het_ab_thresh_snp;   /*! threshold for heterozygous allele balance filter for SNPs */
  float het_ab_thresh_indel;   /*! threshold for heterozygous allele balance filter for INDELs*/
  uint32_t min_het_covg_snp; /*! minimum coverage of heterozygous sites across samples before applying filter */
  uint32_t min_het_covg_indel; /*! minimum coverage of heterozygous sites across samples before applying filter */
  uint32_t min_n_sample;       /*! Number of (het) samples that have to pass the allele balance test*/
  kbitset_t *rm_als;     /*! bitset for keeping track of which alleles to keep and drop */
} args_t;

static args_t args;

const char *about(void)
{
    return
        "Using the AD tag, computes the per-sample coverage of each allele at het sites.\n"
        "Filters variants unless a minumum number of heterozygous samples have both: an allele\n"
        "balance within a specified threshold and a total allele depth above a given threshold. \n"
        "Read depth and allele balance threshold can be specified separately for SNPs and INDELs. \n";
}

const char *usage(void)
{
    return
        "\n"
        "About:   Using the AD tag, this plugin computes the per-sample read depth at het sites, \n"
                  "as well as the proportion of reads that come from each non-ref allele.\n"
        "         Filters variants unless a minumum number of heterozygous samples have both: an allele balance within a specified threshold and a total allele depth above a given threshold. \n"
        "         Depth and allele balance threshold can be specified separately for SNPs and INDELs. \n"      
        "Usage:   bcftools +allelebalhet <multisample.bcf/.vcf.gz> [General Options] -- [Plugin Options] \n"
        "\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "  -s,--mindepth_snp <integer>  min depth of total coverage for het sites across samples for SNPs\n"
        "  -i,--mindepth_indel <integer>  min depth of total coverage for het sites across samples for INDELs\n"
        "  -S,--ab_snp <float>   allele balance threshold for SNPs.  A SNP site must have at least one sample with an allele balance between ab_snp and (1-ab_snp) otherwise the site is excluded\n"
        "  -I,--ab_indel <float>   allele balance threshold for SNPs.  A SNP site must have at least one sample with an allele balance between ab_snp and (1-ab_snp) otherwise the site is excluded\n"
        "  -n,--n_passing <integer>  minimum number of samples that have to pass both the allele balance and allele depth test for the allele to be included (default 1)\n"
        "\n"
        "Example:\n"
        "   bcftools plugin +allelebalhet in.vcf -o out.vcf -- -s 7 -i 10 -S 0.15 -I 0.2 -n 2\n"
        "\n";
}


/*
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  args.het_ab_thresh_snp = 0.15; 
  args.het_ab_thresh_indel = 0.2;
  args.min_het_covg_snp = 7;
  args.min_het_covg_indel = 10;
  args.min_n_sample = 1;

  static struct option loptions[] =  {
      {"mindepth_snp",1,0,'s'},
      {"mindepth_indel",1,0,'i'},
      {"ab_snp",1,0,'S'},
      {"ab_indel",1,0,'I'},      
      {"n_passing",1,0,'n'},      
      {0,0,0,0}
  };
  
  int c;
  while ((c = getopt_long(argc, argv, "s:S:i:I:n:?h",loptions,NULL)) >= 0)    {
    switch (c) {
    case 's': args.min_het_covg_snp = atoi(optarg); break;
    case 'S': args.het_ab_thresh_snp = atof(optarg); break;
    case 'i': args.min_het_covg_indel = atoi(optarg); break;
    case 'n': args.min_n_sample = atoi(optarg); break;
    case 'I': args.het_ab_thresh_indel = atof(optarg); break;      
    case 'h': case '?':
    default: fprintf(stderr,"%s", usage()); exit(1); break;
    }
  }
  //TODO ensure that min_n_sample is <= the number of total samples
  args.het_ab_thresh_snp = fmax(args.het_ab_thresh_snp,1-args.het_ab_thresh_snp);
  args.het_ab_thresh_indel = fmax(args.het_ab_thresh_snp,1-args.het_ab_thresh_indel);
  args.hdr = bcf_hdr_dup(in);  

    
  return 0; // Print header always. 
}


// This is a temporary fix until `bcf_remove_allele_set` is fixed (https://github.com/samtools/htslib/issues/1259)
int bcf_remove_allele_set_2(const bcf_hdr_t *header, bcf1_t *line, const struct kbitset_t *rm_set)
{
    int *map = (int*) calloc(line->n_allele, sizeof(int));
    uint8_t *dat = NULL;

    // create map of indexes from old to new ALT numbering and modify ALT
    kstring_t str = {0,0,0};
    kputs(line->d.allele[0], &str);

    int nrm = 0, i,j;  // i: ori alleles, j: new alleles
    for (i=1, j=1; i<line->n_allele; i++)
    {
        if ( kbs_exists(rm_set, i) )
        {
            // remove this allele
            line->d.allele[i] = NULL;
            nrm++;
            map[i]=-1;  // removed alleles should be set to missing
            continue;
        }
        kputc(',', &str);
        kputs(line->d.allele[i], &str);
        map[i] = j;
        j++;
    }
    if ( !nrm ) goto clean;

    int nR_ori = line->n_allele;
    int nR_new = line->n_allele-nrm;
    if ( nR_new<=0 ) // should not be able to remove reference allele
    {
        hts_log_error("Cannot remove reference allele at %s:%"PRIhts_pos" [%d]",
            bcf_seqname_safe(header,line), line->pos+1, nR_new);
        goto err;
    }
    int nA_ori = nR_ori-1;
    int nA_new = nR_new-1;

    int nG_ori = nR_ori*(nR_ori + 1)/2;
    int nG_new = nR_new*(nR_new + 1)/2;

    bcf_update_alleles_str(header, line, str.s);

    // remove from Number=G, Number=R and Number=A INFO fields.
    int mdat = 0, ndat = 0, mdat_bytes = 0, nret;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *info = &line->d.info[i];
        int vlen = bcf_hdr_id2length(header,BCF_HL_INFO,info->key);

        if ( vlen!=BCF_VL_A && vlen!=BCF_VL_G && vlen!=BCF_VL_R ) continue; // no need to change

        int type = bcf_hdr_id2type(header,BCF_HL_INFO,info->key);
        if ( type==BCF_HT_FLAG ) continue;
        int size = 1;
        if ( type==BCF_HT_REAL || type==BCF_HT_INT ) size = 4;

        mdat = mdat_bytes / size;
        nret = bcf_get_info_values(header, line, bcf_hdr_int2id(header,BCF_DT_ID,info->key), (void**)&dat, &mdat, type);
        mdat_bytes = mdat * size;
        if ( nret<0 )
        {
            hts_log_error("Could not access INFO/%s at %s:%"PRIhts_pos" [%d]",
                bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname_safe(header,line), line->pos+1, nret);
            goto err;
        }
        if ( nret==0 ) continue; // no data for this tag

        if ( type==BCF_HT_STR )
        {
            str.l = 0;
            char *ss = (char*) dat, *se = (char*) dat, s = ss[0];
            if ( vlen==BCF_VL_A || vlen==BCF_VL_R )
            {
                int nexp, inc = 0;
                if ( vlen==BCF_VL_A )
                {
                    nexp = nA_ori;
                    inc  = 1;
                }
                else
                    nexp = nR_ori;
                for (j=0; j<nexp; j++)
                {
                    if ( !*se ) break;
                    while ( *se && *se!=',' ) se++;
                    if ( kbs_exists(rm_set, j+inc) )
                    {
                        if ( *se ) se++;
                        ss = se;
                        continue;
                    }
                    if ( str.l ) kputc(',',&str);
                    kputsn(ss,se-ss,&str);
                    if ( *se ) se++;
                    ss = se;
                }
                if ( j==1 && s == '.' ) continue; // missing
                if ( j!=nexp )
                {
                    hts_log_error("Unexpected number of values in INFO/%s at %s:%"PRIhts_pos"; expected Number=%c=%d, but found %d",
                        bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname_safe(header,line), line->pos+1, vlen==BCF_VL_A ? 'A' : 'R', nexp, j);
                    goto err;
                }
            }
            else    // Number=G, assuming diploid genotype
            {
                int k = 0, n = 0;
                for (j=0; j<nR_ori; j++)
                {
                    for (k=0; k<=j; k++)
                    {
                        if ( !*se ) break;
                        while ( *se && *se!=',' ) se++;
                        n++;
                        if ( kbs_exists(rm_set, j) || kbs_exists(rm_set, k) )
                        {
                            if ( *se ) se++;
                            ss = se;
                            continue;
                        }
                        if ( str.l ) kputc(',',&str);
                        kputsn(ss,se-ss,&str);
                        if ( *se ) se++;
                        ss = se;
                    }
                    if ( !*se ) break;
                }
                if ( n==1 && s == '.' ) continue; // missing
                if ( n!=nG_ori )
                {
                    hts_log_error("Unexpected number of values in INFO/%s at %s:%"PRIhts_pos"; expected Number=G=%d, but found %d",
                        bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname_safe(header,line), line->pos+1, nG_ori, n);
                    goto err;
                }
            }

            nret = bcf_update_info(header, line, bcf_hdr_int2id(header,BCF_DT_ID,info->key), (void*)str.s, str.l, type);
            if ( nret<0 )
            {
                hts_log_error("Could not update INFO/%s at %s:%"PRIhts_pos" [%d]",
                    bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname_safe(header,line), line->pos+1, nret);
                goto err;
            }
            continue;
        }

        if (nret==1) // could be missing - check
        {
            int missing = 0;
            #define BRANCH(type_t, convert, is_missing) { \
                type_t val = convert(info->vptr); \
                if ( is_missing ) missing = 1; \
            }
            switch (info->type) {
                case BCF_BT_INT8:  BRANCH(int8_t, le_to_i8,  val==bcf_int8_missing); break;
                case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, val==bcf_int16_missing); break;
                case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, val==bcf_int32_missing); break;
                case BCF_BT_FLOAT: BRANCH(float,   le_to_float, bcf_float_is_missing(val)); break;
                default: hts_log_error("Unexpected type %d", info->type); goto err;
            }
            #undef BRANCH
            if (missing) continue; // could remove this INFO tag?
        }

        if ( vlen==BCF_VL_A || vlen==BCF_VL_R )
        {
            int inc = 0, ntop;
            if ( vlen==BCF_VL_A )
            {
                if ( nret!=nA_ori )
                {
                    hts_log_error("Unexpected number of values in INFO/%s at %s:%"PRIhts_pos"; expected Number=A=%d, but found %d",
                        bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname_safe(header,line), line->pos+1, nA_ori, nret);
                    goto err;
                }
                ntop = nA_ori;
                ndat = nA_new;
                inc  = 1;
            }
            else
            {
                if ( nret!=nR_ori )
                {
                    hts_log_error("Unexpected number of values in INFO/%s at %s:%"PRIhts_pos"; expected Number=R=%d, but found %d",
                        bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname_safe(header,line), line->pos+1, nR_ori, nret);
                    goto err;
                }
                ntop = nR_ori;
                ndat = nR_new;
            }
            int k = 0;

            #define BRANCH(type_t,is_vector_end) \
            { \
                type_t *ptr = (type_t*) dat; \
                int size = sizeof(type_t); \
                for (j=0; j<ntop; j++) /* j:ori, k:new */ \
                { \
                    if ( is_vector_end ) { memcpy(dat+k*size, dat+j*size, size); break; } \
                    if ( kbs_exists(rm_set, j+inc) ) continue; \
                    if ( j!=k ) memcpy(dat+k*size, dat+j*size, size); \
                    k++; \
                } \
            }
            switch (type)
            {
                case BCF_HT_INT:  BRANCH(int32_t,ptr[j]==bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float,bcf_float_is_vector_end(ptr[j])); break;
            }
            #undef BRANCH
        }
        else    // Number=G
        {
            if ( nret!=nG_ori )
            {
                hts_log_error("Unexpected number of values in INFO/%s at %s:%"PRIhts_pos"; expected Number=G=%d, but found %d",
                    bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname_safe(header,line), line->pos+1, nG_ori, nret);
                goto err;
            }
            int k, l_ori = -1, l_new = 0;
            ndat = nG_new;

            #define BRANCH(type_t,is_vector_end) \
            { \
                type_t *ptr = (type_t*) dat; \
                int size = sizeof(type_t); \
                for (j=0; j<nR_ori; j++) \
                { \
                    for (k=0; k<=j; k++) \
                    { \
                        l_ori++; \
                        if ( is_vector_end ) { memcpy(dat+l_new*size, dat+l_ori*size, size); break; } \
                        if ( kbs_exists(rm_set, j) || kbs_exists(rm_set, k) ) continue; \
                        if ( l_ori!=l_new ) memcpy(dat+l_new*size, dat+l_ori*size, size); \
                        l_new++; \
                    } \
                } \
            }
            switch (type)
            {
                case BCF_HT_INT:  BRANCH(int32_t,ptr[l_ori]==bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float,bcf_float_is_vector_end(ptr[l_ori])); break;
            }
            #undef BRANCH
        }

        nret = bcf_update_info(header, line, bcf_hdr_int2id(header,BCF_DT_ID,info->key), (void*)dat, ndat, type);
        if ( nret<0 )
        {
            hts_log_error("Could not update INFO/%s at %s:%"PRIhts_pos" [%d]",
                bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname_safe(header,line), line->pos+1, nret);
            goto err;
        }
    }

    // Update GT fields, the allele indexes might have changed
    for (i=1; i<nR_ori; i++) 
    if(map[i]!=i)
    {
        mdat = mdat_bytes / 4;  // sizeof(int32_t)
        nret = bcf_get_genotypes(header,line,(void**)&dat,&mdat);
        mdat_bytes = mdat * 4;
        if ( nret>0 )
        {
            nret /= line->n_sample;
            int32_t *ptr = (int32_t*) dat;
            for (i=0; i<line->n_sample; i++)
            {
                for (j=0; j<nret; j++)
                {
                    if ( bcf_gt_is_missing(ptr[j]) ) continue;
                    if ( ptr[j]==bcf_int32_vector_end ) break;
                    int al = bcf_gt_allele(ptr[j]);
                    ptr[j] = (map[al]+1)<<1 | (ptr[j]&1);
                }
                ptr += nret;
            }
            nret = bcf_update_genotypes(header, line, (void*)dat, nret*line->n_sample);
            if ( nret<0 )
            {
                hts_log_error("Could not update FORMAT/GT at %s:%"PRIhts_pos" [%d]",
                    bcf_seqname_safe(header,line), line->pos+1, nret);
                goto err;
            }
        }
    }

    // Remove from Number=G, Number=R and Number=A FORMAT fields.
    // Assuming haploid or diploid GTs
    for (i=0; i<line->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &line->d.fmt[i];
        int vlen = bcf_hdr_id2length(header,BCF_HL_FMT,fmt->id);

        if ( vlen!=BCF_VL_A && vlen!=BCF_VL_G && vlen!=BCF_VL_R ) continue; // no need to change

        int type = bcf_hdr_id2type(header,BCF_HL_FMT,fmt->id);
        if ( type==BCF_HT_FLAG ) continue;

        int size = 1;
        if ( type==BCF_HT_REAL || type==BCF_HT_INT ) size = 4;

        mdat = mdat_bytes / size;
        nret = bcf_get_format_values(header, line, bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), (void**)&dat, &mdat, type);
        mdat_bytes = mdat * size;
        if ( nret<0 )
        {
            hts_log_error("Could not access FORMAT/%s at %s:%"PRIhts_pos" [%d]",
                bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, nret);
            goto err;
        }
        if ( nret == 0 ) continue; // no data for this tag

        if ( type==BCF_HT_STR )
        {
            int size = nret/line->n_sample;     // number of bytes per sample
            str.l = 0;
            if ( vlen==BCF_VL_A || vlen==BCF_VL_R )
            {
                int nexp, inc = 0;
                if ( vlen==BCF_VL_A )
                {
                    nexp = nA_ori;
                    inc  = 1;
                }
                else
                    nexp = nR_ori;
                for (j=0; j<line->n_sample; j++)
                {
                    char *ss = ((char*)dat) + j*size, *se = ss + size, *ptr = ss, s = ss[0];
                    int k_src = 0, k_dst = 0, l = str.l;
                    for (k_src=0; k_src<nexp; k_src++)
                    {
                        if ( ptr>=se || !*ptr) break;
                        while ( ptr<se && *ptr && *ptr!=',' ) ptr++;
                        if ( kbs_exists(rm_set, k_src+inc) )
                        {
                            ss = ++ptr;
                            continue;
                        }
                        if ( k_dst ) kputc(',',&str);
                        kputsn(ss,ptr-ss,&str);
                        ss = ++ptr;
                        k_dst++;
                    }
                    if ( k_src==1 && s == '.' ) continue; // missing
                    if ( k_src!=nexp )
                    {
                        hts_log_error("Unexpected number of values in FORMAT/%s at %s:%"PRIhts_pos"; expected Number=%c=%d, but found %d",
                            bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, vlen==BCF_VL_A ? 'A' : 'R', nexp, k_src);
                        goto err;
                    }
                    l = str.l - l;
                    for (; l<size; l++) kputc(0, &str);
                }
            }
            else    // Number=G, diploid or haploid
            {
                for (j=0; j<line->n_sample; j++)
                {
                    char *ss = ((char*)dat) + j*size, *se = ss + size, *ptr = ss, s = ss[0];
                    int k_src = 0, k_dst = 0, l = str.l;
                    int nexp = 0; // diploid or haploid?
                    while ( ptr<se )
                    {
                        if ( !*ptr ) break;
                        if ( *ptr==',' ) nexp++;
                        ptr++;
                    }
                    if ( ptr!=ss ) nexp++;
                    if ( nexp==1 && s == '.' ) continue; // missing
                    if ( nexp!=nG_ori && nexp!=nR_ori )
                    {
                        hts_log_error("Unexpected number of values in FORMAT/%s at %s:%"PRIhts_pos"; expected Number=G=%d(diploid) or %d(haploid), but found %d",
                            bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, nG_ori, nR_ori, nexp);
                        goto err;
                    }
                    ptr = ss;
                    if ( nexp==nG_ori ) // diploid
                    {
                        int ia, ib;
                        for (ia=0; ia<nR_ori; ia++)
                        {
                            for (ib=0; ib<=ia; ib++)
                            {
                                if ( ptr>=se || !*ptr ) break;
                                while ( ptr<se && *ptr && *ptr!=',' ) ptr++;
                                if ( kbs_exists(rm_set, ia) || kbs_exists(rm_set, ib) )
                                {
                                    ss = ++ptr;
                                    continue;
                                }
                                if ( k_dst ) kputc(',',&str);
                                kputsn(ss,ptr-ss,&str);
                                ss = ++ptr;
                                k_dst++;
                            }
                            if ( ptr>=se || !*ptr ) break;
                        }
                    }
                    else    // haploid
                    {
                        for (k_src=0; k_src<nR_ori; k_src++)
                        {
                            if ( ptr>=se || !*ptr ) break;
                            while ( ptr<se && *ptr && *ptr!=',' ) ptr++;
                            if ( kbs_exists(rm_set, k_src) )
                            {
                                ss = ++ptr;
                                continue;
                            }
                            if ( k_dst ) kputc(',',&str);
                            kputsn(ss,ptr-ss,&str);
                            ss = ++ptr;
                            k_dst++;
                        }
                        if ( k_src!=nR_ori )
                        {
                            hts_log_error("Unexpected number of values in FORMAT/%s at %s:%"PRIhts_pos"; expected Number=G=%d(haploid), but found %d",
                                bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, nR_ori, k_src);
                            goto err;
                        }
                        l = str.l - l;
                        for (; l<size; l++) kputc(0, &str);
                    }
                }
            }
            nret = bcf_update_format(header, line, bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), (void*)str.s, str.l, type);
            if ( nret<0 )
            {
                hts_log_error("Could not update FORMAT/%s at %s:%"PRIhts_pos" [%d]",
                    bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, nret);
                goto err;
            }
            continue;
        }

        int nori = nret / line->n_sample;
        if ( nori==1 && !(vlen==BCF_VL_A && nori==nA_ori) ) // all values may be missing - check
        {
            int all_missing = 1;
            #define BRANCH(type_t, convert, is_missing) { \
                for (j=0; j<line->n_sample; j++) \
                { \
                    type_t val = convert(fmt->p + j*fmt->size); \
                    if ( !(is_missing)) { all_missing = 0; break; } \
                } \
            }
            switch (fmt->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8, val==bcf_int8_missing); break;
                case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, val==bcf_int16_missing); break;
                case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, val==bcf_int32_missing); break;
                case BCF_BT_FLOAT: BRANCH(float,   le_to_float, bcf_float_is_missing(val)); break;
                default: hts_log_error("Unexpected type %d", fmt->type); goto err;
            }
            #undef BRANCH
            if (all_missing) continue; // could remove this FORMAT tag?
        }

        if ( vlen==BCF_VL_A || vlen==BCF_VL_R || (vlen==BCF_VL_G && nori==nR_ori) ) // Number=A, R or haploid Number=G
        {
            int inc = 0, nnew;
            if ( vlen==BCF_VL_A )
            {
                if ( nori!=nA_ori )
                {
                    hts_log_error("Unexpected number of values in FORMAT/%s at %s:%"PRIhts_pos"; expected Number=A=%d, but found %d",
                        bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, nA_ori, nori);
                    goto err;
                }
                ndat = nA_new*line->n_sample;
                nnew = nA_new;
                inc  = 1;
            }
            else
            {
                if ( nori!=nR_ori )
                {
                    hts_log_error("Unexpected number of values in FORMAT/%s at %s:%"PRIhts_pos"; expected Number=R=%d, but found %d",
                        bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, nR_ori, nori);
                    goto err;
                }
                ndat = nR_new*line->n_sample;
                nnew = nR_new;
            }

            #define BRANCH(type_t,is_vector_end) \
            { \
                for (j=0; j<line->n_sample; j++) \
                { \
                    type_t *ptr_src = ((type_t*)dat) + j*nori; \
                    type_t *ptr_dst = ((type_t*)dat) + j*nnew; \
                    int size = sizeof(type_t); \
                    int k_src, k_dst = 0; \
                    for (k_src=0; k_src<nori; k_src++) \
                    { \
                        if ( is_vector_end ) { memcpy(ptr_dst+k_dst, ptr_src+k_src, size); break; } \
                        if ( kbs_exists(rm_set, k_src+inc) ) continue; \
                        memcpy(ptr_dst+k_dst, ptr_src+k_src, size); \
                        k_dst++; \
                    } \
                } \
            }
            switch (type)
            {
                case BCF_HT_INT:  BRANCH(int32_t,ptr_src[k_src]==bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float,bcf_float_is_vector_end(ptr_src[k_src])); break;
            }
            #undef BRANCH
        }
        else    // Number=G, diploid or mixture of haploid+diploid
        {
            if ( nori!=nG_ori )
            {
                hts_log_error("Unexpected number of values in FORMAT/%s at %s:%"PRIhts_pos"; expected Number=G=%d, but found %d",
                    bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, nG_ori, nori);
                goto err;
            }
            ndat = nG_new*line->n_sample;

            #define BRANCH(type_t,is_vector_end) \
            { \
                for (j=0; j<line->n_sample; j++) \
                { \
                    type_t *ptr_src = ((type_t*)dat) + j*nori; \
                    type_t *ptr_dst = ((type_t*)dat) + j*nG_new; \
                    int size = sizeof(type_t); \
                    int ia, ib, k_dst = 0, k_src; \
                    int nset = 0;   /* haploid or diploid? */ \
                    for (k_src=0; k_src<nG_ori; k_src++) { if ( is_vector_end ) break; nset++; } \
                    if ( nset==nR_ori ) /* haploid */ \
                    { \
                        for (k_src=0; k_src<nR_ori; k_src++) \
                        { \
                            if ( kbs_exists(rm_set, k_src) ) continue; \
                            memcpy(ptr_dst+k_dst, ptr_src+k_src, size); \
                            k_dst++; \
                        } \
                        memcpy(ptr_dst+k_dst, ptr_src+k_src, size); \
                    } \
                    else /* diploid */ \
                    { \
                        k_src = -1; \
                        for (ia=0; ia<nR_ori; ia++) \
                        { \
                            for (ib=0; ib<=ia; ib++) \
                            { \
                                k_src++; \
                                if ( is_vector_end ) { memcpy(ptr_dst+k_dst, ptr_src+k_src, size); ia = nR_ori; break; }  \
                                if ( kbs_exists(rm_set, ia) || kbs_exists(rm_set, ib) ) continue; \
                                memcpy(ptr_dst+k_dst, ptr_src+k_src, size); \
                                k_dst++; \
                            } \
                        } \
                    } \
                } \
            }
            switch (type)
            {
                case BCF_HT_INT:  BRANCH(int32_t,ptr_src[k_src]==bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float,bcf_float_is_vector_end(ptr_src[k_src])); break;
            }
            #undef BRANCH
        }
        nret = bcf_update_format(header, line, bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), (void*)dat, ndat, type);
        if ( nret<0 )
        {
            hts_log_error("Could not update FORMAT/%s at %s:%"PRIhts_pos" [%d]",
                bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname_safe(header,line), line->pos+1, nret);
            goto err;
        }
    }

clean:
    free(str.s);
    free(map);
    free(dat);
    return 0;

err:
    free(str.s);
    free(map);
    free(dat);
    return -1;
}
/* 
kbitset_t** realloc_allele_idx(const bcf1_t *rec,kbitset_t **allele_idx,uint32_t *n_allele_idx){
  //Do we need to add alleles?
  if((*n_allele_idx) < rec->n_allele){
  kbitset_t ** old_idx = allele_idx;
  allele_idx =   (kbitset_t **) realloc(old_idx,sizeof(kbitset_t*)*(rec->n_allele));
  //If we moved the pointer we have to delete the old stuff
  if(old_idx!=allele_idx){
    for(int i=0; i<(*n_allele_idx); i++){
      kbs_destroy(old_idx[i]);
    }
    for(int i=0; i<rec->n_allele; i++){
      allele_idx[i]=kbs_init(rec->n_sample);
    }
  }else{ //We didn't move the pointer, so we only need to 1) clear the old indices and 2) initialize the new allele index    
    for(int i=0; i<*n_allele_idx; i++){
      kbs_clear(allele_idx[i]);
    }
    for(int i=(*n_allele_idx); i<rec->n_allele; i++){
      allele_idx[i]=kbs_init(rec->n_sample);
    }
  }
  *n_allele_idx = rec->n_allele;
  }else{ 
     for(int i=0; i<*n_allele_idx; i++){
      kbs_clear(allele_idx[i]);
    }
  }
return allele_idx;
} */


 bool allele_passes_balance_threshold(const bcf1_t *line,const  bcf_fmt_t *fmt,const int ismpl,const int ial,const int jal,const float balance_threshold, const int64_t min_depth){
// float depth_ratio;
  #define BRANCH_INT(type_t,missing,vector_end) {  \
     const type_t *p = (type_t *) (fmt->p + fmt->size*ismpl); \
     if ( p[ial]==vector_end || p[jal]==vector_end ) return false; \
     if ( p[ial]==missing || p[jal]==missing ) return false; \
     if ( !p[ial] && !p[jal] ) return false; \
     const type_t total_depth = p[ial]+p[jal]; \
     if(total_depth< min_depth) return false; \
     float depth_ratio = (float)p[ial]/(total_depth); \
     depth_ratio = fmax(depth_ratio,1-depth_ratio); \
     if ( (depth_ratio > balance_threshold) ) return false; \
  } /*This is the move we have to do because C doesn't have overloads...*/
  switch (fmt->type) {
      case BCF_BT_INT8: BRANCH_INT(int8_t, bcf_int8_missing, bcf_int8_vector_end); break;
      case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
      case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
      default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, fmt->type); exit(1); break;
  }
  //fprintf(stderr, "est=%f\t%d\t%d\n", depth_ratio, total_depth,ial);
  return true;
  #undef BRANCH_INT
}
/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t *process(bcf1_t *rec)
{
  if ( !rec->n_sample ) return rec; // no samples, skip record.

/* Initialize a full bitset with an element for each allele.
As hets pass the ratio and depth test, remove them from the bitset.
Anything that's left after looping over all the samples gets removed
*/
  if ( !args.rm_als )
    args.rm_als = kbs_init(rec->n_allele);                                         
  else if ( args.rm_als->n_max < rec->n_allele )
    kbs_resize(&args.rm_als, rec->n_allele);
  //kbs_insert((args.rm_als),0);
  kbs_insert_all(args.rm_als);


///args.allele_idx = realloc_allele_idx(rec,args.allele_idx,&args.n_allele_idx);
// Start a counter at 0 for each allele;
uint32_t allele_counter[rec->n_allele];
memset(allele_counter,0,sizeof(uint32_t)*rec->n_allele);


 bcf_fmt_t *gt_fmt_ptr;
 if ( !(gt_fmt_ptr = bcf_get_fmt(args.hdr,rec,"GT")) )
   return rec;
 bcf_fmt_t *ad_fmt_ptr;
 if ( !(ad_fmt_ptr = bcf_get_fmt(args.hdr,rec,"AD")))
   return rec;



  // Scan through GT and AD values for each sample.
  bool any_passed = false;
  for (uint32_t i=0; i< rec->n_sample; i++) {
    int ial, jal=0;
    int gt = bcf_gt_type(gt_fmt_ptr, i, &ial, &jal);
    //kbs_insert(args.allele_idx[ial],i);
    //kbs_insert(args.allele_idx[jal],i);
    if(!(gt== GT_HET_AA || gt == GT_HET_RA)) continue; // If we aren't het with two different alts (HET_AA) or het with ref and alt (HET_RA) we can move on
    const bool is_indel = (bcf_get_variant_type(rec,ial) == VCF_INDEL) ||  (bcf_get_variant_type(rec,jal) == VCF_INDEL);
    const float ab_thresh = is_indel ? args.het_ab_thresh_indel : args.het_ab_thresh_snp;
    const int32_t min_depth = is_indel ? args.min_het_covg_indel : args.min_het_covg_snp;

    if(allele_passes_balance_threshold(rec,ad_fmt_ptr,i,ial,jal,ab_thresh,min_depth)){
      
      allele_counter[ial]++;
      allele_counter[jal]++;
      //TODO figure out what to do if allele is balanced between two non-ref alleles
      if(allele_counter[ial]>=args.min_n_sample){
        if(ial>0){
          any_passed=true;
          kbs_delete((args.rm_als),ial);
        }
      }
      if(allele_counter[jal]>=args.min_n_sample){
        if(jal>0) {
          any_passed=true;
          kbs_delete((args.rm_als),jal);
      }
      }
      
      bool all_passed = true;
      for(int k=0; k<rec->n_allele;k++){
        all_passed = all_passed && allele_counter[k]>=args.min_n_sample;
      }
      if(all_passed){
        return rec;
      }
    }
  }
  if(!any_passed) return NULL;
  bcf_unpack(rec,BCF_UN_ALL);
  if ( bcf_remove_allele_set_2(args.hdr, rec, args.rm_als)!=0 ) // Drop the problematic allele(s) from the row and return it
      error("Failed to subset alleles\n");
  // TODO: This is only a good estimate if there are a minimum number of total reads.
  return rec;
}


/*
    Clean up.
*/
void destroy(void)
{
    bcf_hdr_destroy(args.hdr);
    kbs_destroy(args.rm_als);
}



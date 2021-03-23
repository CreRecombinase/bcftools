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
        "Filters variants unless a minumum number of heterozygous samples have both: an allele balance within a specified threshold and a total allele depth above a given threshold. \n"
        "Depth and allele balance threshold can spe specified separately for SNPs and INDELs. \n";
}

const char *usage(void)
{
    return
        "\n"
        "About:   Using the AD tag, computes the per-sample coverage of each allele at het sites.\n"
        "         Filters variants unless a minumum number of heterozygous samples have both: an allele balance within a specified threshold and a total allele depth above a given threshold. \n"
        "         Depth and allele balance threshold can spe specified separately for SNPs and INDELs. \n"      
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
  while ((c = getopt_long(argc, argv, "s:S:i:n:?h",loptions,NULL)) >= 0)    {
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
  kbs_insert_all(args.rm_als);

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
    if(!(gt== GT_HET_AA || gt == GT_HET_RA)) continue; // If we aren't het with two different alts (HET_AA) or het with ref and alt (HET_RA) we can move on
    const bool is_indel = (bcf_get_variant_type(rec,ial) == VCF_INDEL) ||  (bcf_get_variant_type(rec,jal) == VCF_INDEL);
    const float ab_thresh = is_indel ? args.het_ab_thresh_indel : args.het_ab_thresh_snp;
    const int32_t min_depth = is_indel ? args.min_het_covg_indel : args.min_het_covg_snp;

    if(allele_passes_balance_threshold(rec,ad_fmt_ptr,i,ial,jal,ab_thresh,min_depth)){
      
      allele_counter[ial]++;
      allele_counter[jal]++;
      //TODO figure out what to do if allele is balanced between two non-ref alleles
      if(allele_counter[ial]>=args.min_n_sample){
        if(ial>0)
          any_passed=true;
        kbs_delete((args.rm_als),ial);
      }
      if(allele_counter[jal]>=args.min_n_sample){
        if(jal>0) 
          any_passed=true;
        kbs_delete((args.rm_als),jal);
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
  if ( bcf_remove_allele_set(args.hdr, rec, args.rm_als)!=0 ) // Drop the problematic allele(s) from the row and return it
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



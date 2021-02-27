/*  plugins/allelebalhet.c -- allele balannce on hets

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

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <htslib/vcf.h>

typedef struct _args_t
{
  bcf_hdr_t *hdr;        /*! VCF file header */
  uint32_t *gt_arr;     /*! temporary array, to store GTs of current line/record */
  int ngt_arr;          /*! hold the number of current GT array entries */
  uint32_t *ad_arr;     /*! temporary array, to store ADs of current line/record*/
  int nad_arr;          /*! hold the number of current AD array entries */
  float snp_ab_thresh;  
} args_t;

static args_t args;

const char *about(void)
{
    return
        "Filter sites with heterozygous genotype calls that.\n";
}

const char *usage(void)
{
    return
        "\n"
        "About:   \n"
        "         \n"
        "Usage:   bcftools +allelebalhet <multisample.bcf/.vcf.gz> [General Options] -- [Plugin Options] \n"
        "\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "  -ad,---a   \n"
        "             \n"
        "\n"
        "Example:\n"
        "   bcftools +allelebalhet in.vcf -- -adcut 0.1 \n"
        "\n";
}


/*
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  args.hdr = bcf_hdr_dup(in);

  args.gt_arr = NULL; args.ngt_arr = 0;
  args.ad_arr = NULL; args.nad_arr = 0;
    
  return 0; // Print header always. 
}


/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t *process(bcf1_t *rec)
{
  if ( !rec->n_sample ) return rec;
  
  int ngts = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.ngt_arr);
  if ( ngts<0 ) return rec; // no genotypes, ouput record

  if(ngts / rec->n_sample != 2) return rec; // strange ploidy, just leave untouched for now

  int nad = bcf_get_format_int32(args.hdr, rec, "AD", &args.ad_arr, &args.nad_arr);
  if ( nad<0 ) return rec; // no AD field, output record

  if(ngts != nad) return rec; // AD values don't match GTs, output record

  // Scan through GT and AD values for each sample.
  uint32_t i = 0, ad0_het = 0, ad1_het = 0;
  for (i=0; i< rec->n_sample; i++) {
    uint32_t *gt = args.gt_arr + i*2;
    uint32_t *ad = args.ad_arr + i*2;

    // Are we done? 
    if ( gt[0] == bcf_int32_vector_end ) break;
    if ( gt[1] == bcf_int32_vector_end ) break;

    // Not a het if any one is missing.
    if(bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1])) continue; 

    // TODO: How to handle multi-allelics? 
    int allele0 = bcf_gt_allele(gt[0]);
    int allele1 = bcf_gt_allele(gt[1]);

    if(allele0 != allele1) {
      // Found a het site, count reads from the two alleles.
      ad0_het += ad[0];
      ad1_het += ad[1];
    }
  }

  // Total number of reads from het sites must be > 0
  uint32_t total_het = ad0_het + ad1_het;
  if(total_het > 0) {
    float est = (float)ad0_het / total_het;
    int failed_flag = 0;
    fprintf(stderr, "est=%f %d\n", est, total_het);

    // TOOD: This is only a good estimate if there are a minimum number of total reads.
    if(est > 0.7 || est < 0.3) failed_flag = 1;
  
    if(failed_flag) return NULL; else return rec;
  } else {
    return rec;
  }
}


/*
    Clean up.
*/
void destroy(void)
{
    bcf_hdr_destroy(args.hdr);
    free(args.gt_arr);
    free(args.ad_arr);
}



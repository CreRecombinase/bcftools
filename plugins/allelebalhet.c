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
  bcf_hdr_t *hdr;     /*! VCF file header */
  int nsmpl;          /*! number of samples, can be determined from header but is needed in multiple contexts */
  uint32_t *gt_arr;        /*! temporary array, to store GTs of current line/record */
  int ngt_arr;        /*! hold the number of current GT array entries */
  uint32_t *ad_arr;        /*! temporary array, to store ADs of current line/record*/
  int nad_arr;        /*! hold the number of current AD array entries */
} args_t;

static args_t args;

const char *about(void)
{
    return
        "Filter sites with bad allele balance.\n";
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
  args.nsmpl = bcf_hdr_nsamples(in);
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
  fprintf(stderr, "ngts=%d\n", ngts);
  if ( ngts<0 ) return rec; // no genotypes, ouput record

  if(ngts / args.nsmpl != 2) return rec; // strange ploidy, just leave untouched for now

  int nad = bcf_get_format_int32(args.hdr, rec, "AD", &args.ad_arr, &args.nad_arr);
  if ( nad<0 ) return rec; // no AD field, output record
  fprintf(stderr, "nad=%d\n", nad);

  if(ngts != nad) return rec; // AD values don't match GTs
  
  uint32_t i = 0;
  for (i=0; i<args.nsmpl; i++) {
    uint32_t *gt = args.gt_arr + i*2;
    uint32_t *ad = args.ad_arr + i*2;
    
    if ( gt[0]==bcf_int32_vector_end ) break;
    if ( gt[1]==bcf_int32_vector_end ) break;

    // Not a het if any one is missing.
    if(bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1])) continue; 
      
    int allele_index0 = bcf_gt_allele(gt[0]);
    int allele_index1 = bcf_gt_allele(gt[1]);

    if(allele_index0 != allele_index1) {
      fprintf(stderr, "%d(%d),%d(%d) ", allele_index0, ad[0], allele_index1, ad[1]);
    }
  }
  fprintf(stderr, "\n");
  
  
  
  return rec;
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



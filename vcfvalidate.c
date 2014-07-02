/* The MIT License

   Copyright (c) 2014 Genome Research Ltd.
   Authors:  see http://github.com/samtools/bcftools/blob/master/AUTHORS

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <unistd.h>
#include <getopt.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/faidx.h>
#include "bcftools.h"

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    char **argv, *targets_list, *regions_list, *ref_fname;
    faidx_t *faidx;
    int argc, check_counts;
    int *info_ac, info_mac, *fmt_ac, fmt_mac;
}
args_t;

static void destroy_data(args_t *args)
{
    free(args->info_ac);
    free(args->fmt_ac);
}

int vcf_validate(args_t *args, bcf1_t *line)
{
    int i;
    // check for parse errors
    if ( line->errcode ) 
    {
        if (line->errcode&BCF_ERR_CTG_UNDEF)
            fprintf(stderr,"Contig not defined in header found at %s:%d.\n", bcf_seqname(args->hdr,line), line->pos+1);
        if (line->errcode&BCF_ERR_TAG_UNDEF)
            fprintf(stderr,"FILTER, INFO or FORMAT tags not defined in header found at %s:%d.\n", bcf_seqname(args->hdr,line), line->pos+1);
        if (line->errcode&BCF_ERR_NCOLS)
            fprintf(stderr,"Number of columns does not match number of samples in the header at %s:%d (%d vs %d).\n", bcf_seqname(args->hdr,line), line->pos+1, line->n_sample, bcf_hdr_nsamples(args->hdr));
        return -1;
    }
    
    // check INFO field lengths
    for (i=0; i<line->n_info; i++)
    {
        bcf_unpack(line, BCF_UN_INFO);
        bcf_info_t *info = &line->d.info[i];
        int vlen = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key);
        
        int type = bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key);
        if ( type==BCF_HT_FLAG ) {
            if (info->len>0) {
                fprintf(stderr, "INFO/%s is Type=Flag in header, but found %d values at: %s:%d\n", 
                    bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key), info->len, bcf_seqname(args->hdr,line), line->pos+1);        
                return -1;
            }
            continue;
        }
        int size = 1;
        if ( type==BCF_HT_REAL || type==BCF_HT_INT ) size = 4;

        uint8_t *dat = NULL;
        int mdat = 0, mdat_bytes = 0, nfound;
        mdat = mdat_bytes / size;
        nfound = bcf_get_info_values(args->hdr, line, bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key), (void**)&dat, &mdat, type);
        mdat_bytes = mdat * size;
        if ( nfound<0 ) 
        { 
            fprintf(stderr,"Could not access INFO/%s at %s:%d [%d]\n", 
                bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key), bcf_seqname(args->hdr,line), line->pos+1, nfound); 
            return -1;
        }
        if ( type==BCF_HT_STR )
        {
            int j, svlen = 0;
            char *str = (char*) dat;
            for (j=0; j<nfound; j++) {
                if (str[j]==',')
                    svlen++;
            }
            nfound = svlen+1;
        }
        
        if ( vlen==BCF_VL_A && nfound!=line->n_allele-1 )
        {
            fprintf(stderr, "INFO/%s is Number=A in header. Unexpected number of fields at: %s:%d .. %d vs %d\n", 
                bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key), bcf_seqname(args->hdr,line), line->pos+1, nfound, line->n_allele-1);        
            return -1;
        }
        else if ( vlen==BCF_VL_R && nfound!=line->n_allele )
        {
            fprintf(stderr, "INFO/%s is Number=R in header. Unexpected number of fields at: %s:%d .. %d vs %d\n", 
                bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key), bcf_seqname(args->hdr,line), line->pos+1, nfound, line->n_allele);        
            return -1;
        }
        else if ( vlen==BCF_VL_FIXED && nfound!=bcf_hdr_id2number(args->hdr,BCF_HL_INFO,info->key) )
        {
            fprintf(stderr, "INFO/%s is Number=%d in header, but %d fields found at: %s:%d\n", 
                bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key), bcf_hdr_id2number(args->hdr,BCF_HL_INFO,info->key), nfound, bcf_seqname(args->hdr,line), line->pos+1);        
            return -1;
        }
        else if ( vlen==BCF_VL_G )
        {
            // TODO: guess ploidy based on length of GT -> fill in expected number for G
            int ng = line->n_allele*(line->n_allele + 1)/2;
            if (nfound!=ng && nfound!=line->n_allele) { // diploid or haploid for the moment
                fprintf(stderr, "INFO/%s is Number=G in header. Unexpected number of fields at: %s:%d .. %d vs %d\n", 
                    bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key), bcf_seqname(args->hdr,line), line->pos+1, nfound, ng);        
                return -1;                
            }
        }
    }

    // check FORMAT field lengths
    for (i=0; i<line->n_fmt; i++)
    {
        bcf_unpack(line, BCF_UN_FMT);
        bcf_fmt_t *fmt = &line->d.fmt[i];
        // int vlen = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id);

        int type = bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id);
        if ( type==BCF_HT_FLAG ) {
            if (fmt->n>0) {
                fprintf(stderr, "FORMAT/%s is Type=Flag in header, but found %d values at: %s:%d\n", 
                    bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id), fmt->n, bcf_seqname(args->hdr,line), line->pos+1);        
                return -1;
            }
            continue;
        }
        // TODO: the rest of FORMAT
    }

    // check REF matches reference fasta
    if (args->faidx)
    {
        char *ref = line->d.allele[0];
        int ref_len = strlen(ref);
        int fai_ref_len;
        char *fai_ref = faidx_fetch_seq(args->faidx, bcf_seqname(args->hdr,line), line->pos, line->pos+1, &fai_ref_len);
        for (i=0; i<fai_ref_len; i++)
            if ( (int)fai_ref[i]>96 ) fai_ref[i] -= 32;

        for (i=0; i<fai_ref_len && i<ref_len; i++)
            if ( ref[i] != fai_ref[i] && ref[i] - 32 != fai_ref[i] )
            {
                fprintf(stderr, "The reference sequence differs at: %s:%d+%d .. %c vs %c\n", bcf_seqname(args->hdr,line), line->pos+1, i, ref[i], fai_ref[i]);        
                return -1;
            }
    }

    // check INFO/AC and INFO/AN agree with genotypes
    if (args->check_counts && bcf_hdr_nsamples(args->hdr)>0)
    {
        hts_expand(int, line->n_allele, args->info_mac, args->info_ac);
        if (bcf_calc_ac(args->hdr, line, args->info_ac, BCF_UN_INFO)==1) { // AC,AN in INFO if exist
            hts_expand(int, line->n_allele, args->fmt_mac, args->fmt_ac);
            bcf_calc_ac(args->hdr, line, args->fmt_ac, BCF_UN_FMT); // AC,AN in genotypes
            for (i=0; i<line->n_allele; i++)
            {
                if (args->info_ac[i]!=args->fmt_ac[i])
                {
                    fprintf(stderr, "Allele counts differ for counts based on INFO (AC,AN) and GT at: %s:%d .. REF:%d vs %d", bcf_seqname(args->hdr,line), line->pos+1, args->info_ac[i], args->fmt_ac[i]);
                    int j;
                    for (j=1; j<line->n_allele; j++)
                        fprintf(stderr, ", ALT%d:%d vs %d", j, args->info_ac[j], args->fmt_ac[j]);
                    fprintf(stderr, "\n");
                    return -1;
                }
            }
        } 
    }
    return 0;
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   VCF validation. This is a first implementation meant to quickly catch\n");
    fprintf(stderr, "         common errors in VCF files that can prevent conversion to BCF such as\n");
    fprintf(stderr, "         trailing space, INFO/FORMAT fields missing from the header wrong\n");
    fprintf(stderr, "         number of values in vector fields. Not a full validation.\n");
    fprintf(stderr, "Usage:   bcftools validate [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c, --check-counts                  check INFO/AC and INFO/AN counts\n");
    fprintf(stderr, "    -f, --fasta-ref <file>              faidx indexed reference sequence file to check REF allele\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -r, --regions <region>              restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>           restrict to regions listed in a file\n");
    fprintf(stderr, "    -t, --targets <region>              similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>           similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfvalidate(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    int targets_is_file = 0, regions_is_file = 0;

    static struct option loptions[] = 
    {
        {"check-counts",0,0,'c'},
        {"fasta-ref",1,0,'f'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "cf:r:R:t:T:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'c': args->check_counts = 1; break;
            case 'f': args->ref_fname = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage(args);
    }
    else fname = argv[optind];

    // read in the regions from the command line
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }

    if (args->ref_fname)
    {
        args->faidx = fai_load(args->ref_fname);
    }

    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open or the file not indexed: %s\n", fname);
    
    int exit_status = 0;
    args->hdr = args->files->readers[0].header;
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = args->files->readers[0].buffer[0];
        if (vcf_validate(args, line)<0)
            exit_status = 1;
    }
    destroy_data(args);
    bcf_sr_destroy(args->files);
    return exit_status;
}

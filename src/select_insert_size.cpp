#include <stdio.h>
#include "sam.h"
 
void usage(const char *msg) {
    fputs(msg, stderr);
    fputs("\nUSAGE: bam_file insert_size_min_inclusive insert_size_max_inclusive filtered_bam_file\nNote that insert size is negative for overlapping reads!", stderr);
}



int main(int argc, char **argv) {
    if (argc != 4) {
        usage("missing argument");
        return 2;
    }
    int insert_min, insert_max;
    if (sscanf(argv[1], "%i", &insert_min)!=1) {printf ("error argv[1]- not an integer");}
    if (sscanf(argv[2], "%i", &insert_max)!=1) {printf ("error argv[1]- not an integer");}
    samfile_t *bam = samopen(argv[0], "rb", NULL);
    if (bam == NULL) {
        fprintf(stderr, "Could not open %s\n", argv[0]);
        return 2;
    }
    samfile_t *filtered_bam = samopen(argv[3], "wb", NULL);
    if (filtered_bam == NULL) {
        fprintf(stderr, "Could not open %s\n", argv[3]);
        return 2;
    }
    bam1_t *aln = bam_init1();
    int input_count, filtered_count = 0;
    while (samread(bam, aln) >= 0) {
        const bam1_core_t *c = &aln->core;
        input_count++;
        if ((c->flag & BAM_FPROPER_PAIR) == BAM_FPROPER_PAIR &&
            (c->flag & BAM_FSECONDARY) != BAM_FSECONDARY &&
            c->isize >= insert_min && c->isize <= insert_max)
        {
            samwrite(filtered_bam, aln);
            filtered_count++;
        }
    }
    printf("Total reads in input bam file: %i\n", input_count);
    printf("Filtered reads with insert size in [%i,%i]: %i\n", insert_min, insert_max, filtered_count);
    samclose(bam);
    samclose(filtered_bam);
    return 0;
}
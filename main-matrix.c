#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <errno.h>

#include "random.h"
#include "lib.h"
#include "netsim.h"

int intcmp(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

int main(int argc, char *argv[])
{
  FILE *fpkdis;
  char fperrors_name[80];
  int i, j, k;
  Genotype indiv;
  float initProteinConc[NGENES], kdis[NUM_K_DISASSEMBLY];

  int c, directory_success;
  int hold_genotype_constant = 0;
  int curr_seed;

  verbose = 0;

  /* change to get a different genotype */
  dummyrun = 4;

  /* parse command-line options */
  while ((c = getopt (argc, argv, "hvgd:r:p:t:c:")) != -1) {
    switch (c)
      {
      case 'd':
        output_directory = optarg;
        break;
      case 'r':
        dummyrun = atoi(optarg);
        break;
      case 'p':
        current_ploidy = atoi(optarg);
        break;
      case 't':
        tdevelopment = atof(optarg);
        break;
      case 'c':
        critical_size = atof(optarg);
        break;
      case 'g':
        hold_genotype_constant = 1;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        fprintf(stderr, "%s [-d DIRECTORY] [-r DUMMYRUN] [-h] [-g] [-p PLOIDY] [-t DEVELOPMENTTIME] [-c CRITICALSIZE]\n", argv[0]);
        exit(0);
        break;
      default:
        abort();
      }
  }

/* create output directory if needed */
#ifdef __unix__
  directory_success = mkdir(output_directory, S_IRUSR|S_IWUSR|S_IXUSR);
#else 
#ifdef __WIN32__
  directory_success = mkdir(output_directory);
#endif
#endif

  if (directory_success==-1) 
    if (errno == EEXIST) {
      fprintf(stderr, "directory '%s' already exists\n", output_directory);
    } else {
      fprintf(stderr, "directory '%s' cannot be created\n", output_directory);
      exit(-1);
    }

  sprintf(fperrors_name, "%s/netsimerrors.txt", output_directory);
  fperrors = fopen(fperrors_name, "w");

  /* get the kdis.txt values */
  fpkdis = fopen("kdis.txt","r");
  for (j = 0; j < NUM_K_DISASSEMBLY; j++) {
    fscanf(fpkdis,"%f", &kdis[j]);
  }
  fclose(fpkdis); 

  /* change random number */
  for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);


  /**********************************************************************
   *  modifiable part starts here
   **********************************************************************/

  /* initialize protein and mRNA concentrations */
  /* Jasmin: you can use these initial concentrations for computing your kons */
  for (i=0; i<NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
  }

  /* create sequences and binding site matrix */
  initialize_genotype(&indiv, kdis);
  
  /* print binding sites */
  /*print_all_binding_sites(indiv.ploidy, indiv.allBindingSites, indiv.bindSiteCount, 
			  indiv.transcriptionFactorSeq, indiv.cisRegSeq); 
  printf("tfsPerGene = %d", indiv.tfsPerGene);

  /* Jasmin:
     pseudo-code for looping (see also print_all_binding_sites in netsim.c):*/
     //system("PAUSE");
    int sitePos[10];
    int transFactor[10];
    system("PAUSE");
    struct AllTFBindingSites *rearranged;
    rearranged = malloc(indiv.tfsPerGene[0]*sizeof(struct AllTFBindingSites));
    rearranged[0] =indiv.allBindingSites[0];
    printf("%d\n", indiv.allBindingSites[0].leftEdgePos);
    printf("%d\n", rearranged[0].leftEdgePos);
    //int *startPos = calloc(indiv.bindSiteCount, sizeof(int));
    system("PAUSE");
    //printf("tfsPerGene = %d\n", indiv.tfsPerGene[0]);
    int *storeLEP = calloc(10,sizeof(int));
    storeLEP[0]=indiv.allBindingSites[0].leftEdgePos;
    int count =1;
     int front = 1;
     int back = 10;
     while(front < back){
         if(indiv.allBindingSites[front].leftEdgePos<indiv.allBindingSites[back].leftEdgePos){
            if(storeLEP[front-1]>indiv.allBindingSites[front].leftEdgePos){
              storeLEP[count]=storeLEP[count-1];
              storeLEP[count-1]=indiv.allBindingSites[front].leftEdgePos;
              //storeLEP[count+1]=indiv.allBindingSites[back].leftEdgePos;
              printf("%d\n", storeLEP[count-1]);
              printf("%d\n", storeLEP[count]);
              //printf("%d\n", storeLEP[count++]);
            }
            else{
              storeLEP[count]=indiv.allBindingSites[front].leftEdgePos;
              printf("%d\n", storeLEP[count]);
              //storeLEP[count++]=indiv.allBindingSites[back].leftEdgePos;
                 }
            front++;
            count++;
         }else{
            printf("%d\n", indiv.allBindingSites[back].leftEdgePos);
             printf("%d\n", indiv.allBindingSites[front].leftEdgePos);
             back--;
         }
    
       
     }
     system("PAUSE");
     /*int ch;
     for(ch=0; ch<indiv.tfsPerGene[0];ch++){
               printf("%d\n", rearranged[ch].leftEdgePos);
     }*/
    // system("PAUSE");
     /*for(i=0; i < indiv.bindSiteCount; i++){
       if(indiv.allBindingSites[i].cisregID == 0 && i<10){
         printf("numBS=%d\n", i);
         printf("regID=%d\n", indiv.allBindingSites[i].cisregID);
         printf("  position=%d\n", indiv.allBindingSites[i].leftEdgePos);
         printf(" transcription-factor: %3d\n", indiv.allBindingSites[i].tfID);
         //printf("  p=%d, tf=%d\n", indiv.allBindingSites[i].leftEdgePos, indiv.allBindingSites[i].tfID);
         //transFactor[i] = indiv.allBindingSites[i].tfID;
         sitePos[i]=indiv.allBindingSites[i].leftEdgePos;
         }
       }
     
    /*qsort(sitePos, 10, sizeof(int), intcmp);

     int m;
     for(j=0;j<10;j++){
       printf("%d\n", sitePos[j]);
       for(m=0; m<10;m++){
         if(sitePos[j] == indiv.allBindingSites[m].leftEdgePos){
            transFactor[j]= indiv.allBindingSites[m].tfID;
         }
       }
    // printf("%d\n", initProteinConc[j]);
     }
      printf("\n");
     for(j=0;j<10;j++){
     //printf("%d\n", sitePos[j]);
     printf("%d  %d   %.3f\n", sitePos[j], transFactor[j], initProteinConc[(transFactor[j])]);
     
     }
     printf("\n");*/
      
  system("PAUSE");
  /* free dynamically allocated all binding sites list */
  free(indiv.allBindingSites);
  free(rearranged);
  free(storeLEP);
  
  /* close error file */
  fclose(fperrors);
}

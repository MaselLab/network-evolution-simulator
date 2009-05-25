/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Computes probability cloud of occupancy of transcription factors in cis-regulatory regions
 * Authors: Jasmin Uribe, Alex Lancaster
 * Copyright (c) 2008, 2009 Arizona Board of Regents (University of Arizona)
 */

#ifdef __WIN32__
#include <windows.h>
#endif

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

#include "engine.h"

#include "random.h"
#include "lib.h"
#include "netsim.h"

#define BUFSIZE 250

int array_size = 1000;

FILE *sparseMatrixV1;
FILE *statesV1;
FILE *columnV1;
FILE *left_edge_positions;
FILE *bVector;

typedef int (*compfn)(const void*, const void*);

int intcmp(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

/*Reorders binding sites from smallest left_edge_position to largest left_edge_position */
int compare(struct AllTFBindingSites  *elem1, struct AllTFBindingSites *elem2){
    if(elem1->left_edge_pos < elem2->left_edge_pos)
        return -1;
    else if(elem1->left_edge_pos > elem2->left_edge_pos)
         return 1;
    else
         return 0;
}

/*Holds the row number and value of non-zero entry in matrix*/
struct Rowtype {
    long rownum;
    float *kval;
};

/*Holds the column number and all the Rowtypes in that column*/   
struct Ttype {
  int col;
  struct Rowtype *row;
  int rowCount;
};  

/*Converts the given binary configuration (*bits) into a decimal*/
unsigned long convertToDecimal(int *bits, int TFBS){
  int n = TFBS-1;
  int count =0;
  int record[TFBS];
  int rec=0;
  while( n>=0 ){
    if(bits[n]==0){
      count++;
    }else{
      record[rec]=(TFBS-1)-n;
      rec++;
    }
    n--;
  }
  int i;
  unsigned long decimal=0;
  for(i=0; i<rec;i++){
    decimal+= (long)pow(2,record[i]);
  }
  return decimal;
}              

/*Determines whether the given site (bindSite) is hindered given the binary configuration (*bits) 
  and the start positions of all other sites (*startPos)
  returns 1 if site is hindered
  0 if site not hindered*/
int isHindered(int bindSite, int *bits, int *startPos){
  int s = bindSite-1;
  int count=0;
  int check = startPos[bindSite]-(HIND_LENGTH);
  if(check<0){
    check=0;}
  while(s>=0 && startPos[s] <= startPos[bindSite] && startPos[s] >= check){
    if(bits[s] == 1){
      return 1;
    } else {
      s--;
      count++;
    }
  }
  return 0;
}

/*Generates all possible configurations given a start site (bindSite), the number of binding sites (TFBS), and 
  the start positions of all binding sites (*startPos)
  *statesArray is an array of all the configurations */ 
void configure(int bindSite, int *bits, int *numStates, unsigned long *statesArray, int TFBS, int *startPos){
  
  if(bindSite<TFBS-1){
    bits[bindSite] = 0;
    configure(bindSite+1,bits,numStates,statesArray,TFBS,startPos);
    if(!isHindered(bindSite, bits,startPos)){
      bits[bindSite] = 1;
      configure(bindSite+1,bits,numStates,statesArray,TFBS,startPos);
    }
  } else {
    bits[TFBS-1] = 0;
    statesArray[(*numStates)] = convertToDecimal(bits, TFBS);
    (*numStates)++;
    int i;
    for(i=0; i<TFBS; i++){
      printf("%d", bits[i]);
    }
    printf("\n");
    if(!isHindered(bindSite, bits,startPos)){
      bits[TFBS-1]=1;
      //convertToDecimal(bits, TFBS);
      //printf("%d\n", TFBS);
      statesArray[(*numStates)] = convertToDecimal(bits, TFBS);
      //printf("numstates=%d, %lu convert=%lu\n", (*numStates), statesArray[(*numStates)], convertToDecimal(bits, TFBS));
      (*numStates)++;
      for(i=0; i<TFBS; i++){
	printf("%d", bits[i]);
      }
      printf("\n");
    }
  } 
}

/*Populates the diagonal entries in the transition matrix*/
void diagonal(int col, float *diag, struct Ttype *arrayT, int m, int n){
  
  int x;
  float value=0;
  diag[col]=0;
  for (x=0; x<m; x++) {
    value -= *(arrayT[n].row[x].kval);
  }
  diag[col]=value;
  arrayT[n].row[m].rownum = col;
  arrayT[n].row[m].kval = &(diag[col]);     
} 

/*Populates non-daigonal entries in the transition matrix, then calls diagonal() to populate the diagonal entries.
  All rates are stored in arrayT, which is an array of the Ttype structs or an array of column vectors*/  
void transitions(int startSite,int size, unsigned long *viableStates, int TFBSites,
                 struct Ttype *arrayT, float kon[], float koff[5],int *hammDist, float *diag, int *TFon){
  /*States and column files used for debugging and testing. statesV1.txt has the list of states in decimal, 
      columnV1 has the list of columns.*/
  statesV1 = fopen("statesV1.txt", "w");
  columnV1 = fopen("columnV1.txt","w");
  if ((statesV1 = fopen("statesV1.txt", "w"))) {
    if ((columnV1 = fopen("columnV1.txt", "w"))) {
      int i, p, j, m, a;
      unsigned long row;
      int n=0;
      printf("size=%d TFBS= %d\n", size, TFBSites);
      for( i=0;i<size; i++){
        arrayT[n].col = i;
        arrayT[n].row = malloc((array_size)*sizeof(struct Rowtype));
        printf("viableStates:%lu, col num:%d\n",viableStates[i], i);
	    fprintf(statesV1,"%lu \n", viableStates[i]);
	    fprintf(columnV1,"%d\n",i);
        m=0;
        row = viableStates[i];
        for(p=0;p<TFBSites;p++){
          if(!(row & (1L<<p))){   
	              //R-C SWITCH: row = viableStates[i] ^ (1L<<p);    
             row = viableStates[i] | (1L<<p);
	         if(row!=0 && row!=viableStates[i]){
                 //R_C SWITCH: if(row!=viableStates[i]){
	           for(j=0;j<size; j++){
		         if(row==viableStates[j]){
		           arrayT[n].row[m].rownum = j;
		           a = hammDist[p];
		           arrayT[n].row[m].kval = &(koff[a]);
		           m++;
		         }
	           }
	         }
          }else{
	            //R-C SWITCH: row = viableStates[i]|(1<<p);
	            //R-C SWITCH: if(row!=0 && row!=viableStates[i]){
	        row = viableStates[i] ^ (1L<<p);
	        if(row != viableStates[i]){     
	          for( j=0;j<size; j++){   
		        if(row==viableStates[j]){
		          arrayT[n].row[m].rownum = j;
		          a = hammDist[p];
		          arrayT[n].row[m].kval = &(kon[p]);
		          m++;
		        }
              }
            }
          }
        }
	  diagonal(i,diag, arrayT, m, n);
	  m++;
	  arrayT[n].rowCount = m;
	  // printf("rowCount=%d\n",  arrayT[n].rowCount);
	  n++;
    }
  }}
  fclose(columnV1);
  fclose(statesV1);
}

/*Prints the entries in arrayT. Prints the column, row, configuration in binary of column, 
  configuration in binary of row and rate */   
void print_arrayT(struct Ttype *arrayT, int size, unsigned long *viableStates){
  int p, q;  
  printf("%d\n",size); 
  p=0;  
  while (p < size) {
    q=0;
    while (q < arrayT[p].rowCount) {
      printf( "%d  %lu | %lu  %lu  %.2f\n",p,(arrayT[p].row[q].rownum), viableStates[p], viableStates[arrayT[p].row[q].rownum],*arrayT[p].row[q].kval);   
      q++;
    }
    p++; 
    //printf("%d\n", size);  
  }
  printf("\n");
} 

/*Prints column, row, and rate to a text file to pass to Matlab*/
void print_arrayT_MATLAB(struct Ttype *arrayT, int size, unsigned long *viableStates){
  //First part prints to the console
  int p, q;  
  printf("\n Col   Row      kval\n"); 
  p=0;  
  while (p < size) {
    q=0;
    while (q < arrayT[p].rowCount) {
      /*R-C SWITCH : if(p!=4){
	printf("%d     %lu      %f\n",p,arrayT[p].row[q].rownum,   *arrayT[p].row[q].kval); 
	}else{        
	printf("%d     %lu       %d\n",p, arrayT[p].row[q].rownum,  1); 
	}*/
      if(arrayT[p].row[q].rownum!=4){
	printf("%d     %lu      %f\n",p,arrayT[p].row[q].rownum,   *arrayT[p].row[q].kval);                             
      }else{
	printf("%d     %lu       %d\n",p, arrayT[p].row[q].rownum,  1); 
      }
      q++;
    }
    p++;
  }
  int count=0;
  while(count<size){
    printf(" %d     %d    %d\n", count, 4, 1);
    count++;
  }
  //Second part prints to sparseMatrixV1.txt
  sparseMatrixV1 = fopen("sparseMatrixV1.txt", "w");
  if ((sparseMatrixV1 = fopen("sparseMatrixV1.txt", "w"))) {
    int p=0;  
    while (p < size) {
      int q=0;
      while (q < arrayT[p].rowCount) {
	/*R-C SWITCH: if(p!=4){
	  fprintf(sparseMatrixV1, "%ld,   %d,   %f\n" ,arrayT[p].row[q].rownum +1,p+1,  *arrayT[p].row[q].kval);
	  }*/
	if(arrayT[p].row[q].rownum!=4){
	  fprintf(sparseMatrixV1, "%ld,   %d,   %f\n" ,arrayT[p].row[q].rownum +1,p+1,  *arrayT[p].row[q].kval);
	}else{} 
	q++;
      }
      p++;
      }
    int count=0;
    while(count<size){
      //R-C SWITCH: fprintf(sparseMatrixV1," %d     %d    %d\n",  count+1,5, 1);            
      fprintf(sparseMatrixV1," %d     %d    %d\n", 5, count+1, 1);
      count++;
    }
  }
  fclose(sparseMatrixV1);
}

/*Prints the "0" vector. All 0 vector except for one entry is 1 to make sure entries in x vector sum to 1*/
void print_vector_MATLAB(int size){
  int i;
  bVector = fopen("bVector.txt", "w");
  if ((bVector = fopen("bVector.txt", "w"))) {
    for(i=0; i<size; i++){
      if(i!=4){
	fprintf(bVector, "%d\n", 0);
      }else{
	fprintf(bVector, "%d\n", 1);
      }
    }
    fclose(bVector);
  }
}

/*Determines if system is Little Endian or Big Endian
  returns 1 if little
  0 if big*/
int TestByteOrder()
{
  short int word = 0x0001;
  char *byte = (char *) &word;
  return byte[0];
  //return(byte[0] ? LITTLE_ENDIAN : BIG_ENDIAN);
}

/*Converts a text file into a vector of floats*/
void convertFile( char *fileName, float *vect, int size){
  int n;
  FILE *file;
  file = fopen(fileName, "r");
  if(file==NULL){
    printf("Unable to open file!");
  } else {
    for(n=0; n<size; n++){
      int success = fscanf(file, "  %f", &(vect[n]));
      if (!success) {
	printf("The file ran out of input. Check the size parameter to convertFile.\n");              
      }
    }
    fclose(file);
  }
}

/*Counts how many different sites begin at the same base-pair*/
int countMultiples(unsigned long *statesArray, int *startPos, int startSite){
  int start = startPos[0];
  //int count=0;
  //printf("sP[0]=%d, sP[1]=%d, sp[2]=%d, sp[3]=%d\n", startPos[0], startPos[1], startPos[2], startPos[3]);
  //printf("\nStart=%d startPos+1=%d startPos=%d\n", startSite,startPos[1],startPos[0]);
  int count=1;
  while(start==startPos[count]){
    count++;
    startSite++;
  } 
  return count;
}

/*Partitions configurations according to left-most sites*/
void leftMostSites1(int size, unsigned long *states, int mult, int **array, int *count){
  int mask = 01;
  int i,p,h,j;
  int count1;
  j=0;
  
  for(i=0; i<size; i++){
    if((states[i]) % (int)pow(2, mult)==0){
      array[0][j] = states[i];
      printf("%lu\n",states[i]);
      // printf("   %d\n",array[0][j]);
      j++;
    }
  }
  count[0]=j;
  printf("count[0] =%d\n", count[0]);
  
  for(h=0; h<mult; h++){
    count1=0;  
    for(i=0; i<size; i++){
      p = states[i] >> h; 
      if(mask & p){
        array[h+1][count1]=states[i];
	    printf("%d\n", array[h+1][count1]);
	    count1++;
      }
    }
    count[h+1]=count1;
    printf("count1.1=%d\n", count[h+1]);
  }
}

/*Computes final probabilities according to left-most sites partition*/
void probSlide1(unsigned long *statesArray, float *prob, float *outcome, int size, float *previous, int mult, float **array, int *count){
  int mask = 01;
  int i,p,h,j;
  int count1;
  j=0;
  float sumCheck;
  sumCheck=0;
  printf("\n");
  for(i=0; i<size; i++){
    if((statesArray[i]) % (int)pow(2, mult)==0){
      array[0][j] = (prob[i]);
      sumCheck += prob[i];
      //printf("%f\n",prob[i]);
      printf("%d %d   %f\n",0,j,array[0][j]);
      j++;
    }
  }
  count[0]=j;
  printf("count[0] =%d\n", count[0]);
  printf("sumCheck =%f\n", sumCheck);
  
  for(h=0; h<mult; h++){
    count1=0;  
    for(i=0; i<size; i++){
      p = statesArray[i] >> h; 
      if(mask & p){
	    array[h][count1]=prob[i];
	    printf("%d %d  %f\n",(h), count1 ,array[h][count1]);
	    count1++;
      }
    }
    count[h]=count1;
    printf("count1.1=%d\n", count[h]);
  }
  printf("\n");
  //printf("array[1][0]=%f\n", array[0][1]);
  //printf("count[0]=%d\n", count[0]);
  printf("\n");
  
  system("PAUSE");
  
  //Add up probabilities
  int y,x;
  float sum;
  outcome[0] = sumCheck*(*previous);
  printf("previous[0]=%f\n", (*previous));
  printf("sumCheck=%f sumCheck*prev=%f\n", sumCheck, outcome[0]);
  for(y=0; y<mult; y++){
    sum=0;
    for(x=0; x<(count[y]); x++){
      sum += (array[y][x]);
      //printf("%f\n", array[y][x]);
      //printf("    %f\n", sum);
    }
    //printf("outcome[%d] = %f == %f\n", y, outcome[y], outcome[y]*previous[0]);
    outcome[y+1]= sum*(*previous);//normalize to previous
    printf("sum=%f sum*prev=%f\n", sum, outcome[y+1]);
  } 

    *previous = outcome[0];
       
}  

/*Finds the next start site*/
int nextStartSite(int startSite, int *startPos){
  int next;
  int basePairNum;
  basePairNum = startPos[startSite];
  printf("startSite = %d\n basePairNum = %d\n",startSite, basePairNum);
  next = basePairNum;
  int h = startSite+1;
  while(next==startPos[h]){
    h++;
  }
  return h;
}

void populateFinal1(float *outcome,int mult, float *final, int nextSite){
  int i;
  final[0] = outcome[0];
  for(i=0;i<(mult); i++){
    final[nextSite+i] = outcome[i+1];
    //printf("final[%d]=%f\n", nextSite+i, final[nextSite+mult]);
    //printf("outcome[%d]=%f\n", nextSite+i, outcome[i+1]);
  }
}


int main(int argc, char *argv[])
{
  /*NOTE: To make following code compatible with all platfroms and structures need to implement code
    for 4 types: 
    - 32bit Little-Endian
    - 32bit Big-Endian
    - 64bit Little-Endian
    - 64bit Big-Endian
    Current code written for 32bit Little-Endian
  */
  
  //Check Endianess of system
  int testEndian = TestByteOrder();
  if(testEndian==1){printf("%d = LittleEndian\n\n", testEndian);}
  else{printf("%d = BigEndian--WARNING: This code will not work correctly on BigEndian systems!\n\n", testEndian);}
  
  //Check platform of system
  int testPlatform = sizeof(long);
  if(testPlatform == 4){printf("32-bit : Number of Binding Sites cannot exceed 32sites\n\n");}
  else{printf("64-bit\n\n");}
  
  FILE *fpkdis;
  char fperrors_name[80];
  int i, j;
  Genotype indiv;
  float initProteinConc[NGENES], kdis[NUM_K_DISASSEMBLY];

  int c, directory_success;
  int curr_seed;
  int numBp;
  int TFBS;

  numBp = 40;//Size of first sliding window in basse pairs
  
  //Opens Matlab engine
  Engine *ep;

#ifdef __unix__
  const char *startcmd = "/usr/local/MathSoft/bin/matlab -nojvm -nosplash -nodisplay -nodesktop";
#else 
#ifdef __WIN32__
    const char *startcmd = NULL;
#endif
#endif
 
  if(!(ep= engOpen(startcmd))){
      printf("Problem running Matlab.");
      system("PAUSE");
      exit(1);
  } 
  
  verbose = 0;

  /* change to get a different genotype */
  dummyrun = 4;
  current_ploidy = 1;

  /* parse command-line options */
  while ((c = getopt (argc, argv, "hvd:r:p:w:a:")) != -1) {
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
      case 'w':
        numBp = atoi(optarg);
        break;
      case 'a':
        array_size = atoi(optarg);
        break;        
     case 'v':
        verbose = 1;
        break;
      case 'h':
        fprintf(stderr, "%s [-d DIRECTORY] [-r DUMMYRUN] [-h] [-p PLOIDY] [-w WINDOWSIZE] [-a ARRAYSIZE]\n", argv[0]);
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

  if (directory_success==-1) {
    if (errno == EEXIST) {
      fprintf(stderr, "directory '%s' already exists\n", output_directory);
    } else {
      fprintf(stderr, "directory '%s' cannot be created\n", output_directory);
      exit(-1);
    }
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

  /* Initialize protein and mRNA concentrations */
  /* Jasmin: you can use these initial concentrations for computing your kons */
  
  /*Creates array of inial protein concentrations for computig kons*/
  for (i=0; i<NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    //printf("%f\n", initProteinConc[i]);
  }

    
   /*Populate Koff. There are 5 possible koff rates: koff w/0 mismatches, koff w/1 mismatch, koff w/2 mismatches, 
      Koff w/cooperativity on one side, koff w/cooperativity on 2 sides. Cooperativity is not implemented yet.*/
  float Koff[5];
  float Gibbs;
  float RTlnKr;
  float koffCheck;
  float temperature = 293.0;
  RTlnKr = GASCONSTANT * temperature * log(KR);
  
  int jo, jojo;
  for( jo=0; jo<3; jo++){
    Gibbs = (((float) jo)/3.0 - 1.0) * RTlnKr;
    koffCheck = NUMSITESINGENOME*kon*0.25/exp(-Gibbs/(GASCONSTANT*temperature));
    Koff[jo] = koffCheck;
  }
  
  for(jojo=0; jojo<3; jojo++){
    printf("Koff[%d] = %f\n", jojo, Koff[jojo]);
  }
  //*koff = NUMSITESINGENOME*kon*0.25/exp(-Gibbs/(GASCONSTANT*state->temperature));
  //Gibbs = (((float) all_binding_sites[k].hamming_dist)/3.0 - 1.0) * state->RTlnKr;
  //Koff = [0 mismatch, 1 mismatch, 2 mismatch, coop on 1 side, coop on 2 sides]
    
  /* new API requires a Genotype clone */
  Genotype *UNUSED_clone = NULL;
  /* create sequences and binding site matrix */
  initialize_genotype(&indiv, UNUSED_clone, kdis, 0);
  
  /*Sort binding sites from smallest left_edge_position to largest left_edge_position*/
  qsort((void *) &(indiv.all_binding_sites[0]), indiv.sites_per_gene[0],                                 
            sizeof(struct AllTFBindingSites),(compfn)compare );
         printf("\n");
  
  /*Find number of binding sites completely within sliding window*/       
  TFBS =0;
  while(((indiv.all_binding_sites[TFBS].left_edge_pos)+ (HIND_LENGTH-1)) < numBp){
     TFBS++;
  }
  printf("TFBS=%d\n", TFBS);
  system("PAUSE");
  
  /* print binding sites */
  /* print_all_binding_sites(indiv.copies, indiv.all_binding_sites, indiv.bindSiteCount, 
     indiv.transcriptionFactorSeq, indiv.cisRegSeq, indiv.tfsStart); 
			  printf("sites_per_gene = %d", indiv.sites_per_gene); */
  
  /*Create arrays to store information in indiv.all_binding_sites*/
  int *left_edge_pos, *startPos, *hammDist, *TFon; 
  float *diag, *final, previous;
  float Kon[TFBS];//populate Kon. Kons depend on TF concentrations
  
  int *bits = calloc(TFBS, sizeof(int));//array of bits to create configurations
  unsigned long *viableStates;//array of viable states
  struct Ttype *arrayT;//array of column vectors for transition matrix
  
  
  int finalPos =0;//tracks how far final vector has been populated
  int f,lem, bob, n;//counters
  int startSite, posNum, totalBp, nextS, numStates, mult;
  /*startSite - tracks left most binding site or sliding window
    posNum - next avaliable spot in final vector
    totalBp - tracks how many base-pairs have been considered
    nextS - next site to be considered
    numStates - number of possible configurations
    mult - how many sites begin at a binding site*/
  
  int **leftArray;//partitions configurations by left-most sites
  float **probArray;//partitions probability vector
  int *countArray;//tracks how many configurations in each partition
  float *outcome, *vector;/*outcome - probability outcome of each iteration
			    vector - probability vector returned by matlab*/
  
  //allocate memory for these arrays
  left_edge_pos = malloc(indiv.sites_per_gene[0]*sizeof(int));
  startPos=malloc(TFBS*sizeof(int));
  hammDist = malloc(TFBS *sizeof(int));
  diag = malloc(array_size*sizeof(float));
  TFon = malloc(TFBS*sizeof(int));
  final = malloc(array_size*sizeof(float));
  //previous = malloc(1*sizeof(float));
  
  viableStates = malloc((array_size)*sizeof(unsigned long));
  arrayT = malloc(array_size*sizeof(struct Ttype));
  // arrayT = malloc((pow(2,TFBS))*sizeof(struct Ttype));
     
  //First iteration normalizes to 1
  previous=1;
  /*float prev[10] = {1,0,0,0,0,0,0,0,0,0};
  
  for (pi=0; pi<10; pi++){
    previous[pi]=prev[pi];
  }*/
  
  /*Prints all left_edge_positions to file. Used for debugging and testing purposes*/    
  left_edge_positions = fopen("left_edge_positions.txt", "w");
  if ((left_edge_positions = fopen("left_edge_positions.txt", "w"))) {
    for( f=0; f<TFBS; f++) {
      printf("binding site %3d:  ", f);
      printf("%d\n", indiv.all_binding_sites[f].left_edge_pos);
    }
    printf("sites_per_gene = %d", indiv.sites_per_gene[0]);
    system("PAUSE");
	 
    //Stores left_edge_positions of each site on the first gene in left_edge_pos[] array and prints to file. 
    //Commented code is to print other information
    for (i=0; i <indiv.sites_per_gene[0] ; i++) {
      fprintf(left_edge_positions, "binding site %3d:  ", i);
      left_edge_pos[i] =  indiv.all_binding_sites[i].left_edge_pos;
      /*fprintf( left_edge_positions, "binding site %3d:  ", i);
	fprintf(left_edge_positions, "       cis-reg region: %3d",indiv.all_binding_sites[i].cisregID);
	fprintf(left_edge_positions, "         cis-reg copy: %3d", indiv.all_binding_sites[i].geneCopy);
	fprintf(left_edge_positions, " (sequence %.*s)\n", CISREG_LEN, indiv.cisRegSeq[indiv.all_binding_sites[i].cisregID][indiv.all_binding_sites[i].geneCopy]);
	fprintf(left_edge_positions, " transcription-factor: %3d", indiv.all_binding_sites[i].tf_id);
	fprintf(left_edge_positions, " (sequence: %.*s)\n", TF_ELEMENT_LEN, indiv.transcriptionFactorSeq[indiv.all_binding_sites[i].tf_id][indiv.all_binding_sites[i].geneCopy]); 
	fprintf(left_edge_positions, "  L-edge of %2dbp hind: %3d\n", HIND_LENGTH, indiv.all_binding_sites[i].left_edge_pos);        
        
      */
      fprintf(left_edge_positions, "%d\n", indiv.all_binding_sites[i].left_edge_pos);
      
      /*//fprintf(left_edge_positions,  "%d\n", indiv.all_binding_sites[i].left_edge_pos);
	fprintf(left_edge_positions, "  Hind offset position: %3d\n", indiv.all_binding_sites[i].hindPos); 
	fprintf(left_edge_positions, "               strand: %3d\n", indiv.all_binding_sites[i].strand);
	fprintf(left_edge_positions, "         Hamming dist: %3d\n\n", indiv.all_binding_sites[i].hamming_dist); 
      */}
  }
  fclose(left_edge_positions);
  system("PAUSE");
  printf("\n");
  
  startSite=0;//start at left- most site
  
  //Main while loop for sliding window #1
  while (startSite<indiv.sites_per_gene[0]) {
    numStates =0;//initialize number of configurations
    
    //Find number of sites within sliding window. Changes with startSite
    TFBS = 0;
    while(((indiv.all_binding_sites[TFBS].left_edge_pos)+(HIND_LENGTH-1)) < numBp+startSite){
      TFBS++;
    }
    printf("TFBS inside=%d\n", TFBS);
    system("PAUSE");
    
    //populate startpos, hammDist, TFon, Kon for specified number of binding sites 
    for (lem =0; lem<TFBS; lem++) {
      startPos[lem] = indiv.all_binding_sites[lem+startSite].left_edge_pos;
      hammDist[lem] = indiv.all_binding_sites[lem+startSite].hamming_dist;
      TFon[lem] = indiv.all_binding_sites[lem+startSite].tf_id;
      bob = indiv.all_binding_sites[lem+startSite].tf_id;
      Kon[lem] = initProteinConc[bob]*kon;
      
      printf("%d", indiv.all_binding_sites[lem+startSite].left_edge_pos);
      printf(" Hd = %d   tf = %d", hammDist[lem], TFon[lem]);
      printf(" Kon[lem] = %f\n", Kon[lem]);    
    }
    system("PAUSE");
    printf("\n");
     
    //Generate states for given binding sites and hinderances
    configure(0,bits,&numStates,viableStates,TFBS,startPos);
    
     //printf("HERE array=%d\n", array);
    printf("startSite=%d\n", startSite);
    //count how many different sites begin at same place (usually only 1)
    mult=countMultiples(viableStates, startPos, startSite);
    printf("\n MULT = %d\n", mult);
   
    //allocate memory to leftArray, probArray and countArray
    leftArray = malloc((mult+1)*sizeof(int *));
    probArray = malloc((mult+1)*sizeof(float *));
    for(i=0;i<(mult+1);i++){
       leftArray[i]= malloc(100*sizeof(int));
       probArray[i]= malloc(100*sizeof(float));
    }
    countArray= malloc((mult+1)*sizeof(int));
    system("PAUSE"); 
    printf("%d\n", numStates);

    //Print the "0" vector to a text file
    print_vector_MATLAB(numStates);
    system("PAUSE");
     
    //populate transition matrix
    transitions(startSite,numStates,viableStates,TFBS,arrayT, Kon, Koff, hammDist, diag, TFon);
    system("PAUSE");
     
    //print non-sero entries to text file  
    print_arrayT_MATLAB(arrayT,numStates,viableStates);
    //print_arrayT(arrayT, array, viableStates);
    system("PAUSE");
    
    if(!(ep= engOpen(startcmd))){
              printf("Problem running Matlab.");
              system("PAUSE");
              exit(1);
     }
     
    //matlab call to m-file pVectRevised
    engEvalString(ep,"pVectRevised");
    //engClose(ep); 
    printf("matlab done\n");     
    
    /*if(system("matlab -nodisplay -nojvm -nodesktop -nosplash -r \"pVect; exit;\"")!=0){
      printf("Problem running Matlab.");
      exit(1);   
      }else{
      printf("\n");
      }*/
    system("PAUSE");
    
    //allocate memory for outcome and vector
    outcome = malloc( numStates*sizeof(float));
    vector = malloc(numStates*sizeof(float));
    
    //Get solution vector from matlab and store in vector
    convertFile("b.txt", vector, numStates);
    printf("\nvector\n");
    for(n=0; n<numStates; n++){
            printf( "%f\n", (vector[n]));
    }
    
    //Get probabilities of left most sites being bound or unbound and fix these probabilites
    probSlide1(viableStates, vector, outcome, numStates, &previous, mult, probArray, countArray);
    printf("\n");
    /*for(i=0; i<(mult); i++){
          for(j=0; j<countArray[i]; j++){
                   printf("%f\n", probArray[i][j]);
          }
          printf("\n");
       }*/
    printf("finalPos=%d\n", finalPos);
    
    //find next available slot in the final probability vector
    if (finalPos==0) {
       posNum=1;
    } else {
       posNum = 1+finalPos ;
    }
    
    //Prints the outcome of the collapsed probabilities
    printf("\n OUTCOME \n");
    for (n=0; n<mult+1; n++){
      printf("%f\n", outcome[n]);
    }
    printf("\n");
    printf("posNum=%d\n", posNum);//next available position in final vector
    
    
    //Populate the final vector of probabilities with fixed probabilities
    populateFinal1(outcome, mult, final, posNum);
    
    //Prints the final fixed probabilities
    printf("\nFINAL\n");
    for(i=0; i<(posNum+mult); i++){
      printf("%f\n", final[i]);
    }
    printf("\n");
    system("PAUSE");
    
    
    //Find the next starting posittion
    printf("previous Start Site= %d\n", startSite);
    nextS = nextStartSite(startSite, left_edge_pos);
    
    finalPos+=mult;
    
    //Finds the total number of base pairs 
    totalBp=numBp;
    startSite = nextS;
    printf("startSiteHERE = %d finalPos=%d\n", startSite, finalPos);
    totalBp+=startPos[startSite]-1;
    printf("totalBp=%d\n", totalBp);
    printf("\n\n END OF ITERATION!\n\n");
    system("PAUSE");  
  }//end of main while loop
  
  engClose(ep);
  
  /* free dynamically allocated all binding sites list */
  free(indiv.all_binding_sites);
  int d;
  for (d=0; d<numStates; d++) {
    free(arrayT[d].row);
  } 
  for(d=0;d<mult; d++){
    free(leftArray[d]);
    free(probArray[d]);
  }  
  
  
  free(arrayT);
  free(startPos);
  free(viableStates);
  free(bits);
  free(hammDist);
  free(TFon);
  free(diag);
  free(leftArray);
  free(probArray);
  free(outcome);
  free(vector);
  free(countArray);
//  free(previous);
  /* close error file */
  fclose(fperrors);
}

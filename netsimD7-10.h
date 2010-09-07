/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-  */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN

#include <stdio.h>

//#ifndef POP_SIZE
#define POP_SIZE 1         /* default to a single cell if not otherwise defined */
//#endif

#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define EPSILON 1e-6       /* original code used EPSILON 10^-6 */
#define RT_SAFE_EPSILON 1e-6
#define TIME_INFINITY 9.99e10
#define RATE_OPERATIONS 1e6   /* number of operations on a rate before recomputing that rate */

#ifndef MAX_COPIES
#define MAX_COPIES 4       /* each gene can exists with four copies during replication */
#endif

#ifdef  NO_SEPARATE_GENE    /* set to true if selection gene is also a TF */
  #ifndef TFGENES
  #define TFGENES 10          /* number of genes encoding TFs */
  #endif
  #ifndef NGENES
  #define NGENES TFGENES      /* total number of cisreg genes */
  #endif
  #ifndef NPROTEINS
  #define NPROTEINS TFGENES   /* total number of types of proteins, may be >#GENES if extracellular signals */
  #endif
  #define SELECTION_GENE (NGENES-1)    /* index of selection gene */
#else                       /* otherwise, by default assuming selection gene is not a TF */
  #ifndef TFGENES             /* number of genes encoding TFs */
  #define TFGENES 10
  #endif
  #ifndef NGENES
  #define NGENES (TFGENES+1)  /* total number of genes: add the extra (non-TF) selection gene to the total (default case) */
  #endif
  #ifndef NPROTEINS           /* number of proteins: TODO: must be equal to the number of genes currently */
  #define NPROTEINS (TFGENES+1)
  #endif
  #define SELECTION_GENE (NGENES-1)   /* index of selection gene: always the last gene */
#endif

#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 6      /* length of binding element on TF */
#define NUM_K_DISASSEMBLY 133 /* number of differents for PIC disassembly from data file  */

#ifndef HIND_LENGTH
#define HIND_LENGTH 15         /* default length of hindrance (original was 6) */
#endif

/* 
 * define macros for logging output and warning/errors 
 */
#ifndef LOGGING_OFF
  #define LOG(...) { fprintf(fperrors, "[%s: cell %03d] ", __func__, genotype->cell_id); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); }
  #define LOG_NOCELLID(...) { fprintf(fperrors, "[%s] ", __func__); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); }
  #define LOG_ERROR(...) { fprintf(fperrors, "[%s: cell %03d ERROR] ", __func__, genotype->cell_id); \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_ERROR_NOCELLID(...) { fprintf(fperrors, "[%s ERROR] ", __func__); fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_WARNING(...) { fprintf(fperrors, "[%s: cell %03d WARNING] ", __func__, genotype->cell_id); \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_WARNING_NOCELLID(...) { fprintf(fperrors, "[%s WARNING] ", __func__); fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_NOFUNC(...) { fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_VERBOSE_NOFUNC(...) if (verbose) { fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_VERBOSE(...) if (verbose) { \
    if (state!=NULL)                  \
      fprintf(fperrors, "[%s: cell %03d] ", __func__, genotype->cell_id); \
    else \
      fprintf(fperrors, "[%s] ", __func__);  \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_VERBOSE_NOCELLID(...) if (verbose) { \
    fprintf(fperrors, "[%s] ", __func__);                 \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#else
  #define LOG
  #define LOG_ERROR
  #define LOG_NOFUNC
  #define LOG_VERBOSE
#endif

/*
 * primary data structures for model
 */

/* 
 * enumeration for 'konvalues' indices
 * are rates of binding with:
 * element KON_DIFF          is (proteinConc - salphc)/proteindecay
 * element KON_PROTEIN_DECAY is proteindecay
 * element KON_SALPHC        is salphc
 */
//enum { KON_DIFF_INDEX = 0, KON_PROTEIN_DECAY_INDEX = 1, KON_SALPHC_INDEX = 2 };

/*
 * enum for 'CellState'->active indices
 */
enum { OFF_FULL = 1,           /* repr>activ, still nuclesome, and PIC */
       ON_WITH_NUCLEOSOME = 2, /* activ>repr, still nucleosome */
       OFF_NO_PIC = 3,         /* repr>activ, no nucleosome, but no PIC */
       ON_NO_PIC = 4,          /* activ>repr, no nucleosome, but no PIC  */
       OFF_PIC = 5,            /* repr>activ, no nucleosome, and PIC */
       ON_FULL = 6 };          /* activ>repr, no nucleosome, and a PIC: ready to go!  */

enum {STATE_ZERO = 0, 
      STATE_ONE = 1,
      STATE_TWO = 2
};

/*
 * enum for state_change_ids
 */
/*enum { ACETYLATION_STATE = 0, 
       DEACETYLATION_STATE = 1, 
       PICASSEMBLY_STATE = 2,
       TRANSCRIPTINIT_STATE = 3, 
       PICDISASSEMBLY_STATE = 4,
};
*/

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int cisreg_id;     /* cis-reg region */
  int tf_id;         /* transcription factor */
  int strand;        /* strand 0 (forward) or 1 (backward)*/
  int hamming_dist;  /* hamming distance */
  int gene_copy;     /* which copy of gene, 0 to MAX_COPIES-1 */
  int left_edge_pos; /* start position of TF on DNA, always with reference to forward strand */
                     /* note that recognition sitePos, can be computed from: (left_edge_pos+hind_pos) */
  int hind_pos;      /* position of recognition site within the HIND_LENGTH bp hindrance (offset) */
};

typedef struct Genotype Genotype;
struct Genotype {
  int cell_id;                        /* cell ID */
  int founder_id;                     /* keep track of the founder cell */
  int burn_in;                        /* keep track of whether burn-in has finished */
  int divisions;                      /* total number of divisions cell has undergone as mother, reset as daughter */
  char cisreg_seq[NGENES][MAX_COPIES][CISREG_LEN];
  char tf_seq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN];
  int hindrance_positions[TFGENES];     /* offset positions of each TF's hindrance area relative to recognition site*/
  int binding_sites_num;                    /* total number of binding sites */
  AllTFBindingSites *all_binding_sites;
  float mRNAdecay[NGENES];            /* kinetic rates*/
  float proteindecay[NGENES];
  float translation[NGENES];
  int activating[NGENES][MAX_COPIES]; /* 1 is activating, 0 is repressing */
  float pic_disassembly[NGENES][MAX_COPIES];
  int copies[NGENES];                 /* current per-gene ploidy */
  float replication_time[NGENES];     /* per-gene replication time within S phase*/

  /* cached quantities for efficiency, can be recomputed from the above genotype, not part of model */
  int sites_per_gene[NGENES];               /* cache number of TFBSs per gene */
  int site_id_pos[NGENES][MAX_COPIES][2];  /* cache start and end positions of binding siteIDs for each cis-reg sequence for use during replication and division*/
};

/* 
 * transcription/translation delays are sorted linked lists.  Deleting
 * the head each time, and tack new stuff on the end.  Linked lists
 * are easy to create pre-sorted.
 */
typedef struct FixedEvent FixedEvent;
struct FixedEvent {
  int gene_id;
  int copy;
  float time;
  FixedEvent *next;
};

typedef struct CellState CellState;
struct CellState {
  int in_s_phase;                     /* whether cell has entered S (synthesis) phase */
  float division_time;                /* in global time, TIME_INFINITY until start of S phase*/ 
  float cell_size;                    /* size of cell */
  float growth_rate;                  /* total growth rate in the previous deltat */
  int mRNA_cyto_num[NGENES];          /* mRNAs in cytoplasm */
  int mRNA_nuclear_num[NGENES];       /* mRNAs in nucleus */
  int mRNA_transl_cyto_num[NGENES];   /* mRNAs are in the cytoplasm, but only recently */
  FixedEvent *mRNA_transl_time_end;   /* times when mRNAs become fully loaded with ribosomes and start producing protein */
  FixedEvent *mRNA_transl_time_end_last;
  int mRNA_transcr_num[NGENES][MAX_COPIES];  /* mRNAs which haven't finished transcription yet */
  FixedEvent *mRNA_transcr_time_end;  /* times when transcription is complete and an mRNA is available to move to cytoplasm */
  FixedEvent *mRNA_transcr_time_end_last;
  FixedEvent *replication_time_end;   /* times when genes (start and) finish replicating */
  FixedEvent *replication_time_end_last;    

  float protein_conc[NPROTEINS];
  int active[NGENES][MAX_COPIES];
  /* gives the state of each of the genes, according to figure
   *  1 is fully off, 2 meets TF criteria
   *  3 is off but w/o nucleosome, 4 is on but w/o PIC
   *  5 is on but w/o TF criteria, 6 is fully on
   * see the enum definition above
   */

  /* stores corresponding gene_ids ready for [de]acteylation,
     PIC[dis]assembly, transcriptinit, see enum above
     element 0 says which gene copy it is, element 1 says which geneID*/
  int state_change_ids[3][MAX_COPIES*NGENES][2]; 
  float RTlnKr;
  float temperature;

  //merge GillespieRates into CellState
  float transport;
  int transport_operations;
  float mRNAdecay;
  int mRNAdecay_operations;
  float pic_disassembly;
  int pic_disassembly_operations;

   /* number of genes subject to the following transitions, the rates of the
      relevant process are computed from these numbers */
  int acetylation_num[MAX_COPIES];
  int deacetylation_num[MAX_COPIES];
  int pic_assembly_num[MAX_COPIES];
  int transcript_init_num[MAX_COPIES];
  int pic_disassembly_num[MAX_COPIES];

  int stateID_num[3];
  float total;
};

typedef struct TimeCourse TimeCourse;

struct TimeCourse
{
  float concentration;
  float time;
  TimeCourse *next;
};  


typedef struct Dtype Dtype;
struct Dtype{
       float active;
       float repress;
       float ratio;
       int count;
};

typedef struct Wtype Wtype;
struct Wtype{
       int tfbsNum;
       int startPos;
       int hammDist;
       int tfIDon;
       float conc;
       float weight;
}; 

typedef int (*compfn)(const void*, const void*);
//typedef float (*compfnf)(const void*, const void*);






      

/*
 * global variables
 */

/* see netsim.c for documentation for these global constant variables */
extern const int MAXELEMENTS; 
const int MAXBOUND;
const int NMIN;
const float KRNA;
const float TTRANSLATION;
const float TTRANSCRIPTION;
const float PROB_ACTIVATING;
const float TRANSCRIPTINIT; 
const float DEACETYLATE;
const float ACETYLATE;
const float PICASSEMBLY;
const float STARTNUCLEUS;
const float KR;    
const float GASCONSTANT;
const float COOPERATIVITY;
const float COOPERATIVE_DISTANCE; 
const float NUMSITESINGENOME ;
const float NITER;

const float mN ;     

/* see netsim.c for documentation for these global variables */
float kon; 
float kon_after_burnin; 
int burn_in;
float tdevelopment;  
float timemax;       
int current_ploidy;  
int output;
long seed ;          
int dummyrun;        
int recompute_koff;  
int recompute_kon;
float critical_size;
float growth_rate_scaling; 
float Pp;         
float h;
float gmax;
float time_s_phase;          
float time_g2_phase;         
int random_replication_time; 
float protein_aging;
float *active_gene;

/* file output parameters */
char *output_directory ;
int verbose ;
FILE *fperrors;
FILE *fp_cellsize[POP_SIZE];
FILE *recordFile;
FILE *outpuFile;

#if 0 
//FILE *fp_koff[POP_SIZE];
FILE *fp_growthrate[POP_SIZE];
//FILE *fp_tfsbound[POP_SIZE];
FILE *fp_rounding[POP_SIZE];
#endif

/* function prototypes */

extern void initialize_parameters();

extern void initialize_growth_rate_parameters();

extern void initialize_sequence(char [], int, int, int);

extern void print_genotype(Genotype *, int);

extern void print_all_binding_sites(int [NGENES],
                                    AllTFBindingSites *, 
                                    int ,
                                    char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                    char [NGENES][MAX_COPIES][CISREG_LEN],
                                    int [NGENES][MAX_COPIES][2]);

extern void print_tf_occupancy(CellState *,
                               AllTFBindingSites *,
                               float);
extern void initialize_genotype(Genotype *, 
                                Genotype *,
                                float [],
                                int);


extern void mutate(Genotype *,
                   int,
                   int,
                   float);

extern int calc_all_binding_sites_copy(char [NGENES][MAX_COPIES][CISREG_LEN],
                                       char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                       int ,
                                       AllTFBindingSites **,
                                       int *,
                                       int ,
                                       int ,
                                       int [TFGENES]);

extern void calc_all_binding_sites(int [NGENES],
                                   char[NGENES][MAX_COPIES][CISREG_LEN],
                                   char[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                   int *,
                                   AllTFBindingSites **,
                                   int [TFGENES],
                                   int [NGENES],
                                   int [NGENES][MAX_COPIES][2]);

extern int add_fixed_event(int,
                           int,
                           float,
                           FixedEvent **,
                           FixedEvent **);

extern void add_time_point(float,
                           float,
                           TimeCourse **,
                           TimeCourse **);

extern void add_fixed_event_end(int,
                                int,
                                float,
                                FixedEvent **,
                                FixedEvent **);

extern void delete_fixed_event(int,
                               int,
                               int,
                               FixedEvent **,
                               FixedEvent **);

extern void delete_fixed_event_start(FixedEvent **,
                                     FixedEvent **);

extern void initialize_cell(CellState *,
                            Genotype *,
                            int,
                            int [NGENES],
                            float [NGENES],
                            float [NGENES],
                            float [NPROTEINS],
                            int);

extern void initialize_cell_cache(CellState *,
                                  Genotype,
                                  float **,
                                  int,
                                  int);

extern void change_mRNA_cytoplasm(int,
                                  Genotype *,
                                  CellState *);

extern int ready_to_transcribeJU(int,
                                 float *,
                                 CellState *, int);

extern int ready_to_transcribe(int,
                               int,
                               //int *,
                               //int,
                               AllTFBindingSites *,
                               int [NGENES][MAX_COPIES],
                               int *);

extern int is_one_activator(int,
                            int,
                            //int *,
                            //int,
                            AllTFBindingSites *,
                            int [NGENES][MAX_COPIES]);


extern int does_fixed_event_end(FixedEvent *,
                                FixedEvent *,
                                FixedEvent *,
                                float);

extern void end_transcription(float *,
                              float,
                              CellState *,
                              float [NGENES]);


extern void disassemble_PIC(CellState *,
                            Genotype *,
                            int,
                            int);

extern void revise_activity_state(int,
                                  int,
                                  Genotype *,
                                  CellState *,
                                  float *gene_active, float transport);

extern void add_time_points(float,
                            float [NPROTEINS],
                            TimeCourse **,
                            TimeCourse **);

extern void add_integer_time_points(float,
                                    int [NPROTEINS],
                                    TimeCourse **,
                                    TimeCourse **);

extern void reach_s_phase(CellState *, Genotype *, float);

extern float compute_tprime(float, float, float, float);

extern float compute_integral(float, float, float, float, float, float, float, float, float);

extern float compute_growth_rate_dimer(float *,
                                       float , 
                                       float ,
                                       float [NGENES],
                                       int [NGENES],
                                       float,
                                       float,
                                       float, 
                                       float,
                                       float,
                                       float,
                                       float,
                                       int);

extern void update_protein_conc_cell_size(float[],
                                          CellState *,
                                          Genotype *,
                                          float,
                                          float,
                                          TimeCourse **,
                                          TimeCourse **,
                                          float []);

extern int sum_rate_counts(int[MAX_COPIES]);

extern void get_gene(int [MAX_COPIES], int, int *, int *);

extern void transport_event(
                            CellState *,
                            Genotype *,
                            //KonStates *,
                            float [NGENES],
                            TimeCourse **, 
                            TimeCourse **, 
                            float,
                            float,
                            float);


extern void mRNA_decay_event(CellState *, Genotype *,
                             float *, TimeCourse **, TimeCourse **,
                             float, float, float);

extern void histone_acteylation_event(CellState *, Genotype *,
                                      TimeCourse **, TimeCourse **,
                                      float, float);

extern void histone_deacteylation_event(CellState *, Genotype *,
                                        TimeCourse **, TimeCourse **,
                                        float, float);

extern void assemble_PIC_event(CellState *, Genotype *,
                               TimeCourse **, TimeCourse **,
                               float, float);

extern void disassemble_PIC_event(CellState *, Genotype *,
                                  TimeCourse **, TimeCourse **,
                                  float, float, float);

extern void transcription_init_event(CellState *, Genotype *,
                                     TimeCourse **, TimeCourse **,
                                     float, float, float);

extern void replicate_gene(CellState *,
                           Genotype *,
                           int,
                           float);

extern void recalibrate_cell(CellState *,
                             Genotype *,
                             float [NGENES],
                             float [NGENES],
                             float);// float *

extern void initialize_new_cell_genotype(Genotype *, Genotype *);

extern int do_single_timestep(Genotype *, 
                               CellState *,
                               float *,
                               //float *,
                               float [NGENES],
                               float [NGENES],
                               float *,
                               float *,
                               TimeCourse *[NPROTEINS],
                               TimeCourse *[NPROTEINS],
                               int,
                               int,
                               int, 
                               float *) ;
  
extern void init_run_pop(Genotype [POP_SIZE],
                         CellState [POP_SIZE],
                         TimeCourse *[POP_SIZE][NGENES],
                         TimeCourse *[POP_SIZE][NGENES], 
                         float, /* in Kelvin */
                         float [NUM_K_DISASSEMBLY],
                         int,
                         int,
                         int); 

extern void print_time_course(TimeCourse *,
                              int, int);

extern void print_all_protein_time_courses(TimeCourse *[POP_SIZE][NPROTEINS],
                                          TimeCourse *[POP_SIZE][NPROTEINS]);

extern void log_snapshot(
                         CellState *,
                         Genotype *,
                         float [NGENES],
                         float [NGENES],
                         float ,
                         float );
                         
/*extern void active_vect(Genotype *,
                        float [NGENES],
                        float *);*/

//void active_vect(Genotype indiv, float initProteinConc[NGENES], float *gene_active){


#endif /* !FILE_NETSIM_SEEN */

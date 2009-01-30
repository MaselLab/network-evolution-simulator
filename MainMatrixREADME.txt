Basic Outline of Main-Matrix
Jasmin Uribe
Jan 30 2009

main
 -tests Endianness and platform size
 -currently works with 32-bit Little Endian
 -numBp is size of sliding window (can be changed with -w )
 -gets initial protein concentrations from Alex's code (initProteinConc)
 -calculates koff based on number of mismatches 
          koff[0] - 0 mismathces
          koff[1] - 1 mismatch
          koff[2] - 2 mismatches
          koff[3] - cooperativity on one side (NOT IMPLEMENTED)
          koff[4] - cooperativity on 2 sides  (NOT IMPLEMENTED)
 -initializes the genotype using Alex's code
 -sorts binding sites in order of starting position
 -TFBS is number of binding sites completely within window 
 -allocates memory for all arrays
      leftEdgePosition- hold all starting positions of all binding sites on gene
      startPos- starting positions of binding sites within window
      hammDist- hamming distances of all binding sites 
      TFon- transcription factor number (0-9)
      viableStates- array of all viable configurations based on hinderances between sites
      diag- holds the diagonal elements of the transition matrix
      final- the final probabiltiy vector
      previous- the "0" probabiltiy for the previous iteration, used to normalize the current iteration
      
      
      startSite - tracks left most binding site or sliding window
      posNum - next avaliable spot in final vector
      totalBp - tracks how many base-pairs have been considered
      nextS - next site to be considered
      numStates - number of possible configurations
      mult - how many sites begin at a binding site
      
      leftArray- 2 dimensional array that partitions configurations based on left most sites
          rows- different partitions, depends on the number of binding sites starting at the same base-pair
          columns- configurations that belong to each partition
      probArray- 2-d Array that partitions probability vector like leftArray(above)
          rows- different partitions
          columns- probabilites that belong to each partition
      countArray- how many elements in each partition
      
 -prints starting postions of ALL binding sites on gene in leftEdgePositions.txt
 -stores the starting positions of the binding sites in the window on leftEdgePos
       commented stuff is other infomation that comes from the same struct
 -begins with the 0-th binding site (startSite)
 
 -begins Main while-loop for 1st sliding window
 -recomputes number of binding sites within window
 -populates startPos, hammDist, TFon, and Kon for the binding sites within the window
 -finds all possible configurations based on starting positions of binding sites (configure)
 -stores all possible configurations in viableStates
 -counts how different sites begin at the first base-pair (mult)
 -allocates memory for leftArray, probArray and countArray
 -prints the "0" vector to a file for Matlab to read (vector actually has one 1 to keep probabilities adding to 1)
 -Matlab call to pVectRevised that solves Ax = 0, prints x to a text file (b.txt)
 -converts b.txt to a vector in C
 -probSlide1 partitions the x-vector and adds up the probabilities in each partition
 -finalpos is the last occupied position in the final probability vector
 -posNum is the next available position in the final probability vector
 -outcome holds the probabilty of each partition for that iteration
 -populateFinal1 normalizes outcome to the "0" probability of the previous and populates the final probability vector
 -finds the next start site that is not repeat
 -increases the final position
 -increments the total number of binidng sites considered
 -begins Main while-loop again
 
Some stats
 -150 available bp per gene
 -approximately 90 binding sites per gene
 -on the order of 10^8 different configurations for 90 binding sites
 -sparse matrix is less than 20% populated
 -COLUMNS in array must sum to 0
 -each site is 15bp long
 
Still needs to be done
 -Collapsing window
 -once first sliding window gets past 14bp collapse probabilties of final probability vector based on
  number of activators and repressors
  
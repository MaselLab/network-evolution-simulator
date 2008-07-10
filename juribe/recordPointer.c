#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TFBS 5
#define SIZE 1

//static int hind[SIZE][2]={{0,5}};

static int startPos[TFBS]= {0,0,1,4,4};
static int record[4][2]={{0,0},{0,0},{0,0},{0,0}};


int ifPossible(unsigned int e, int hind[SIZE][2]){
    int possible, i, m , b, d;
    unsigned int s, w;
    possible = 1;
    i = SIZE-1;
   
    if (e>0) {
       while (possible==1 && i>=0) {
          d = hind[i][1];
          s = (e >> TFBS - hind[i][0] - hind[i][1]);
          int p= (int) pow(2,d);
          w = s % p;
          if (w>0) {
            m=d;
            b=(1 << d);
            while (m>=0 && w!=b) {
               m--;
               b = b/2;
            }
            if (m<0) {
               possible = 0;
            }
          }
          i--;
       }
    }
    return (possible);   
}

int numStates(int tfbs, int hind[SIZE][2]){
    int i, count;
    
    count=0;
    for(i=0; i<pow(2, tfbs); i++){
      if(ifPossible(i, hind)==1){
          //printf("%d\n", i);
          count++;
      }    
    }
    return(count);
} 

void vectStates(int vect[], int n, int tfbs, int hind[SIZE][2]){
    int i, count, m;
    
    count=0;
    while(count<n){
    for(i=0; i<pow(2, tfbs); i++){
      if(ifPossible(i, hind)==1){
          vect[count]=i;
          count++;
      }    
    }
  }    
}

void baseCase(int *len, int **rec, int *start){


    int n, b, e;
    n=0;
   // while(n<1){
        b=startPos[n];
        e=b+5;
        printf("b=%d e=%d \n", b,e);
        
        int f, count;
      
        count=0;
        f=n;
        while(f<TFBS && (startPos[f]>=b) && (startPos[f]<=e)){
            //printf("startPos=%d\n", startPos[f]);
            count++;
            f++;
        }    
        printf("b=%d e=%d count=%d\n", b,e,count);
        
        rec=realloc(rec,2*count*(sizeof(int)));
        int g, num;
        g=n;
        num = 0;
        while(startPos[g]==startPos[n]){
            num++;
            g++;
        }
        printf("num= %d, count=%d\n", num, count);  
        
        int vectSize, i,p;
        vectSize=0;
        int hinderances[SIZE][2]={{0,count}};
        vectSize=numStates(count, hinderances);
        printf("vectSize=%d\n", vectSize);
        /*int possibleSites[vectSize];
        
        vectStates(possibleSites,vectSize,count,hinderances);*/
        int possibleSites[vectSize];
        
        
        possibleSites[0]=0;
        for(i=0; i<vectSize; i++){
          possibleSites[i+1]= (int)pow(2,i);
        }
    
       for(i=0;i<vectSize;i++){
          printf("%d ", possibleSites[i]);
       }   
      
      int length;
      length = count-num;
      printf("\nlength= %d\n", length); 
      
      int m;
      m=pow(2,length);
      printf("m= %d\n", m);
      
      printf("count= %d\n", count);
      
      int modSites[vectSize];
      for(i=0;i<vectSize;i++){
          modSites[i]= (possibleSites[i])%m;
          printf("%d ", modSites[i]);
      }  
  
      printf("\n\n");
      int check, q, r, v,u,acount;
      //int shortRec[2];
      q=0;
      acount=0;
      
      while(q<vectSize){
          
          check = modSites[q];
          u=0;
          if(check==-1){
              q++;
          } else{    
 
          printf("%d ", check);
          v=0;
          r=0;
         while(r<vectSize){
                        
             if(check==modSites[r]){
                 modSites[r]=-1;
                   v++;
               }
           r++; 
         }
         printf("%d\n", v);
         //printf("q=%d\n", q);
         record[q][0]=check;
         //printf("recq=%d\n", rec[q][0]);
         record[q][1]=v;
          //printf("recq1=%d\n", rec[q][1]);
          system("PAUSE");
       //acount++; 
       q++; 
        //printf("%d ", record[q][1]); 
      }}
     //system("PAUSE");
     printf("acount=%d\n", acount);
     //find next binding site
     q=1;
     check= startPos[n];
     //printf("%d ", check);
     while(q<11){
         if(check==startPos[n+q]){
             q++;
         }
         else{
             break; 
         }    
     } 
    // printf("\n%d ", q); 
        
        n= n+q;
     *len=length;
     *start=n;
}    
int main(int argc, char *argv[])
{
    int a,b;
    int *length, *star;
    int **R;
    a=0;
    b=0;
    R=malloc(5*2*(sizeof(int)));
    length=&a;
    star=&b;
    baseCase(length, R, star);
    
    printf("len=%d star=%d\n", *length, *star);
    printf("R[1]=%d", *(R[1]));
    free(R);
    system("PAUSE");
    }

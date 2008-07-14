#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TFBS 52
#define SIZE 1

//static int hind[SIZE][2]={{0,5}};

static int startPos[TFBS]= {0,0,1,4,4,6,6,7,8,8,9,10,11,14,14,17,19,21,21,22,27,31,31,
                            32,32,34,36,38,38,39,39,41,43,43,46,46,47,48,50,52,57,61,62,
                            64,64,64,65,65,67,67,67,69};
//static int record[4][2]={{0,0},{0,0},{0,0},{0,0}};


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

int numHind(int start){
    int n, b, e;
    n=0;
    b=startPos[start];
    e=b+5;
        
    int f, count;
      
    count=0;
    f=start;
    while(f<TFBS && (startPos[f]>=b) && (startPos[f]<=e)){
        count++;
        f++;
    }
    if(count<5){
       count=5;
   }       
    return(count);
} 
void printR(int *len, int **rec, int *start, int *reclength, int f){
    int h,i;

    printf("len=%d star=%d rlen=%d \n", *len, *start, *reclength);
    printf("f=%d\n",f);
    //printf("rec[f]=%d\n", rec[f-1][1]);
    for(h=0;h<f;h++){
        for(i=0;i<2;i++){
           printf("%d ", rec[h][i]);
       }
       printf("\n");
   }
}    

void baseCase(int *len, int **rec, int *start, int *reclength){


    int n, b, e;
    n=0;
    b=startPos[n];
    e=b+5;
        
    int f, count;
      
    count=0;
    f=n;
    while(f<TFBS && (startPos[f]>=b) && (startPos[f]<=e)){
        count++;
        f++;
    }    
    printf("b=%d e=%d count=%d\n", b,e,count);
        
        
    int g, num;
    g=n;
    num = 0;
    while(startPos[g]==startPos[n]){
         num++;
         g++;
    }
    printf("num= %d, count= %d\n", num, count);  
    system("PAUSE"); 
       
    int vectSize, i,p;
    vectSize=count+1;
    printf("vectSize=%d ", vectSize);
   

   printf("vectSize=%d\n", vectSize);
     system("PAUSE");   
    int *possibleSites;
    possibleSites=malloc(vectSize*sizeof(int));
        
    possibleSites[0]=0;
    for(i=0; i<count; i++){
         possibleSites[i+1]= (int)pow(2,i);
         printf("possibleSites%d = %d\n",i, possibleSites[i]);
    } 
   
    system("PAUSE");
    
    int length;
    length = count-num;
    //printf("\nlength= %d\n", length); 
      
    int m;
    m=pow(2,length);
    //printf("m= %d\n", m);
      
    //printf("count= %d\n", count);
     
     int *modSites;
    modSites=malloc(vectSize*sizeof(int));
     for(i=0;i<vectSize;i++){
         modSites[i]= (possibleSites[i])%m;
        // printf("%d ", modSites[i]);
     }  
  
    // printf("\n\n");
     int check, q, r, v,u,acount;
      
     q=0;
     acount=0;
      
     while(q<vectSize){
          rec[q]=malloc(2*sizeof(int));
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
             rec[q][0]=check;
             rec[q][1]=v;
         
             acount++; 
             q++; 
           }
      }
      //free(hinderances);
      free(modSites);
      free(possibleSites);

     q=1;
     check=startPos[n];
     while(q<10){
         if(check==startPos[n+q]){
             q++;
         }
         else{
             break; 
         }    
     } 
        
     n= n+q;
     *len=length;
     *start=n;
     *reclength=acount;
} 

void genCase(int *prevlen, int **rec, int *start, int *reclength){
    
    int b,e;
    b=startPos[(*start)];
    e=b+5;
    
    int f, count;
    count=0;
    f=*start;
    printf("f=%d, prevlength=%d \n",f, *prevlen); 
    while(f<TFBS && (startPos[f]>=b) && (startPos[f]<=e)){
        count++;
        f++;
    }
    //printf("b=%d e=%d count=%d\n", b,e,count);
    int i;
    printf("REC: "); 
    for(i=0;i<(*reclength);i++){
      printf(" %d ", rec[i][1]);
  }    
    
    
    int length;
    length = count - (*prevlen);
    printf("length=%d\n", length);
    
     //system("PAUSE");  
    int g, num;
    g=*start;
    num = 0;
    while(startPos[g]==startPos[*start]){
       num++;
       g++;
    }
    //printf("num= %d, count=%d\n", num, count);  
        
    int vectSize, p;

    vectSize=count+1;
    printf("vectSize=%d\n", vectSize);
      
     system("PAUSE"); 
    int *possibleSites;
    possibleSites=malloc(vectSize*sizeof(int));
        
    possibleSites[0]=0;
    for(i=0; i<vectSize; i++){
       possibleSites[i+1]= (int)pow(2,i);
    }
 
    for(i=0;i<vectSize;i++){
       //printf("%d ", possibleSites[i]);
    } 
    //printf("\n");
    
    int m;
    int *divSites;
    divSites=malloc(vectSize*sizeof(int));
    m=pow(2,length);
    printf("divsites: ");
      for(i=0;i<vectSize;i++){
          divSites[i]= (possibleSites[i])/m;
          printf("%d ", divSites[i]);
      }
      printf("\n"); 
     
    
    int **tempRec;
    tempRec = malloc(vectSize*sizeof(int));
    int a,z,blah;
    //printf("REC CHECK=%d", rec[vectSize-1][1]);
    for(a=0; a<vectSize; a++){
        tempRec[a]=malloc(2*sizeof(int));
        blah=possibleSites[a]/m;
        z=0;
        while(blah!=rec[z][0] && z < (*reclength)){
            z++;
        }
        if(z<(*reclength)){
           tempRec[a][0]=blah;
           printf("z=%d",z);
           system("PAUSE");
           tempRec[a][1]=rec[z][1];
           system("PAUSE");
           printf("blah=%d, num=%d ", tempRec[a][0], tempRec[a][1]);
           printf("\n");
           //system("PAUSE");
       }    
    }      
    printf("\n");
      //system("PAUSE");
    int nextlen;
    nextlen=count-num;
    //printf("nextlen=%d\n", nextlen);
    

    int *modSites;
    modSites=malloc(vectSize*sizeof(int));
    m=pow(2,nextlen);
    printf("modsites: ");
      for(i=0;i<vectSize;i++){
          modSites[i]= (possibleSites[i])%m;
          printf("%d ", modSites[i]);
      } 
      printf("\n"); 
      
      //system("PAUSE");
    printf("\n");
    int thing, sum,acount; 
    acount=0;
    sum=0;
    //rec=realloc(rec,count*sizeof(int));
    for(a=0;a<nextlen+1; a++){
        rec[a]=malloc(2*sizeof(int));
        thing = modSites[a];
        sum = tempRec[a][1];
        z=a+1;
        while(z<count+1){
            while(thing!=modSites[z] && z<count+1){
                z++;
            }
            if(thing==modSites[z] && z<count+1){
                sum+=tempRec[z][1];
                z++;
            }
            else{
                break;}
                
        }
        rec[a][0]=possibleSites[a];
        rec[a][1]=sum;
        printf("%d %d\n",rec[a][0], rec[a][1]);
        acount++;
    } 
    printf("acoutn=%d\n", acount);
    //rec=realloc(rec,acount*sizeof(int)); 
    //printf("\n");
    free(tempRec);
    {int i;
     for(i=0;i<vectSize;i++){
         free(tempRec[i]);
     }
    }        
   // free(hinderances);
    free(divSites);
    free(modSites);

    int q, check;
     q=1;
     check= startPos[*start];
     while(q<11){
         if(check==startPos[(*start)+q]){
             q++;
         }
         else{
             break; 
         }    
     } 
     
     //printf("%d\n", length);
    *prevlen= nextlen;
    *start=(*start) + q;
    *reclength=acount;          
}                      

  
int main(int argc, char *argv[])
{
    int a,b,c,e,q,f,h, i,j,next;
    int *length, *star, *rlen;
    int **R;
    a=0;
    b=0;
    e=0;
    
    rlen=&e;
    length=&a;
    star=&b;
    c=numHind(0);
    R=malloc(c*(sizeof(int)));
    
   baseCase(length, R, star, rlen);
    
    f=*rlen;
    R=realloc(R,(f)*(sizeof(int)));
    printR(length,R,star,rlen,f);
    next=*star;
    c=numHind(next);
    while(next<47){
       c=numHind(next);
       printf("c=%d next=%d\n",c,next);
       R=realloc(R,(c)*(sizeof(int)));
       printf("BEFORE: len=%d star=%d rlen=%d \n", *length, *star, *rlen);
       genCase(length, R, star, rlen);
       if(*length==0){
          printf("gencase\n");
          f=*rlen;
          printR(length,R,star,rlen,f); 
         
       } else{   
          printf("gencase\n");
          f=*rlen;
          printf("f=%d\n", f);
         
          //R=realloc(R,(f)*(sizeof(int)));
          printR(length,R,star,rlen,f); 
          next=*star;
      }    
   }    
    
  
   /*genCase(length, R, star, rlen);
   printf("gencase\n");
   f=*rlen;
   R=realloc(R,(f)*(sizeof(int)));
   printR(length,R,star,rlen,f);
 
   //system("PAUSE");
   
   genCase(length, R, star, rlen);
   printf("gencase\n");
   f=*rlen;
   printR(length,R,star,rlen,f);
 
   
   genCase(length, R, star, rlen);
   printf("gencase\n");
   f=*rlen;
   printR(length,R,star,rlen,f);*/

   
   for(q=0;q<c;q++){
       free(R[q]);
   }            
    free(R);
    system("PAUSE");
    }

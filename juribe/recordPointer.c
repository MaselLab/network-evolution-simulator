#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TFBS 98
#define SIZE 1
#define HIND_LENGTH 6

//static int hind[SIZE][2]={{0,5}};

static int startPos[TFBS]= /*{1,2,3,3,10,10,11,11,11,13,16,16,17,17,20,20,21,24,24,25,25,25,
                           26,26,26,27,28,28,29,32,36,37,37,41,44,46,48,51,53,53,57,57,57,
                           59,59,61,62,62,63,66,67,69,74,78,79,79,80,80,81,82,84,84,85,86,86,
                           88,88,89,90,90,91,91,92,92,93,94,100,100,100,103,103,104,104,110,
                           116,118,119,119,123,123,125,126,127};

/*{2,3,5,6,6,12,12,13,14,15,15,17,20,27,27,28,30,31,31,31,31,31,
                           32,35,39,41,42,43,44,45,45,46,55,56,57,58,58,59,60,60,60,61,63,63,
                           64,71,76,76,76,77,77,78,80,81,83,84,84,85,90,90,91,91,92,93,94,94,
                           95,96,97,99,100,100,100,102,103,104,104,104,104,106,107,109,110,114,
                           114,114,114,117,117,118,121,121,122,123,123,127,131,132};

/*{0,0,1,4,7,8,11,11,14,14,16,18,19,22,22,23,24};*/

/*{0,0,1,4,4,6,6,7,8,8,9,10,11,14,14,17,19,21,21,22,27,31,31,
                            32,32,34,36,38,38,39,39,41,43,43,46,46,47,48,50,52,57,61,62,
                            64,64,64,65,65,67,67,67,69,71,72,79,79,80,80,81,84,84,86,89,90,
                            91,94,94,94,95,95,95,96,97,98,98,99,100,102,102,105,106,108,108,
                            108,111,112,114,115,116,116,117,118,118,118,121,123,124,126,127,
                            129,130,131,131,139,140,142,144};*/
                            
{5,6,10,10,11,12,12,14,15,17,18,20,20,20,21,21,22,22,25,27,27,28,28,29,29,29,30,30,32,34,37,39,41,43,46,
                           49,49,52,55,56,57,57,57,58,61,61,64,64,66,67,69,72,74,74,76,78,82,83,83,84,84,85,87,88,
                           89,89,90,90,91,92,92,93,93,94,96,97,98,99,101,104,105,107,108,108,112,112,118,118,121,122,
                           124,128,131,134,134,136,141};
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
    int n, b, e, f, count;
    
    n=0;
    b=startPos[start];
    e=b+(HIND_LENGTH-1);
      
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

void printR(int *len, int (*rec)[2], int *start, int *reclength, int f){
    int h,i;

    printf("len=%d star=%d rlen=%d \n", *len, *start, *reclength);
    //printf("f=%d\n",f);

    for(h=0;h<f;h++){
        for(i=0;i<2;i++){
           printf("%d ", rec[h][i]);
       }
       printf("\n");
   }
}    

void baseCase(int *len, int (*rec)[2], int *start, int *reclength){
    int n, b, e, f, g, i, p, q, r, v, u, m, count, num, vectSize, length, check, acount; 
    int *possibleSites;
    int *modSites;
    
    n=0;
    b=startPos[n];
    e=b+(HIND_LENGTH-1);
    
    count=0;
    f=n;
    while(f<TFBS && (startPos[f]>=b) && (startPos[f]<=e)){
        count++;
        f++;
    }        
    //printf("b=%d e=%d count=%d\n", b,e,count);
        
    g=n;
    num = 0;
    while(startPos[g]==startPos[n]){
         num++;
         g++;
    }
    //printf("num= %d, count= %d\n", num, count);     
    
    vectSize=count+1;
    //printf("vectSize=%d ", vectSize);
  
    possibleSites=malloc(vectSize*sizeof(int));      
    possibleSites[0]=0;
    for(i=0; i<count; i++){
         possibleSites[i+1]= (int)pow(2,i);
         printf("possibleSites%d = %d\n",i, possibleSites[i]);
    } 
    
   
    length = count-num;
    m=pow(2,length);  
     
    modSites=malloc(vectSize*sizeof(int));
     for(i=0;i<vectSize;i++){
         modSites[i]= (possibleSites[i])%m;
     }  
      
     q=0;
     acount=0;
     while(q<vectSize){
          check = modSites[q];
          u=0;
          if(check==-1){
              q++;
          } else{    
              v=0;
              r=0;
              while(r<vectSize){       
                  if(check==modSites[r]){
                     modSites[r]=-1;
                     v++;
                  }
                  r++; 
               }
             rec[q][0]=check;
             rec[q][1]=v;
         
             acount++; 
             q++; 
           }
      }

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

void genCase(int *prevlen, int (*rec)[2], int *start, int *reclength){
    int f,i,g,p,m,a,z,q,count,length,num,vectSize,blah,nextlen,thing,sum,acount,check;
    int *possibleSites;
    int *divSites;
    int (*tempRec)[2];
    int *modSites;
    
    {
        int b,e;
        b=startPos[(*start)];
        e=b+(HIND_LENGTH-1);
 
        count=0;
        f=*start;
 
        while(f<TFBS && (startPos[f]>=b) && (startPos[f]<=e)){
             count++;
             f++;
        } 
    }       
    
    /*printf("REC: \n"); 
    for(i=0;i<(*reclength);i++){
      printf("%d %d \n", rec[i][0],rec[i][1]);
    } */   
    
    length = count - (*prevlen);
    //printf("length=%d\n", length);
    
    g=*start;
    num = 0;
    while(startPos[g]==startPos[*start]){
       num++;
       g++;
    }

    vectSize=count+1;
    //printf("vectSize=%d\n", vectSize);
      
    possibleSites=malloc(vectSize*sizeof(int));
        
    possibleSites[0]=0;
    for(i=0; i<vectSize-1; i++){
       possibleSites[i+1]= (int)pow(2,i);
    }
   
    divSites=malloc(vectSize*sizeof(int));
    m=pow(2,length);
   // printf("divsites: ");
      for(i=0;i<vectSize;i++){
          divSites[i]= (possibleSites[i])/m;
          //printf("%d ", divSites[i]);
      }
    //printf("\n"); 
     
    //printf("vectSize=%d", vectSize);
    tempRec = malloc(2*vectSize*sizeof(int));
    
    for(a=0; a<vectSize; a++){
        blah=possibleSites[a]/m;
        z=0;
        while(blah!=rec[z][0] && z < (*reclength)){
            z++;
        }
        if(z<(*reclength) && blah==rec[z][0]){
           tempRec[a][0]=blah;
           tempRec[a][1]=rec[z][1];
           //printf("blah=%d, num=%d ", tempRec[a][0], tempRec[a][1]);
           //printf("\n");
       }    
    }      
    //printf("\n");
    
    nextlen=count-num;
 
    modSites=malloc(vectSize*sizeof(int));
    m=pow(2,nextlen);
    //printf("modsites: ");
      for(i=0;i<vectSize;i++){
          modSites[i]= (possibleSites[i])%m;
          //printf("%d ", modSites[i]);
      } 
      
    printf("\n");
    
    acount=0;
    sum=0;
    rec=realloc(rec,count*2*sizeof(int));
    for(a=0;a<nextlen+1; a++){
        thing = modSites[a];
        sum = tempRec[a][1];
        z=a+1;

        while((z<vectSize)){
	        while((z<vectSize) && thing!=modSites[z]){
                z++;
            }
	        if((z<vectSize) && thing==modSites[z]){
                sum+=tempRec[z][1];
                z++;
            }
            else{
                break;}
                
        }
        rec[a][0]=possibleSites[a];
        rec[a][1]=sum;
        //printf("%d %d\n",rec[a][0], rec[a][1]);
        acount++;
    } 
   // printf("acoutn=%d\n", acount);

    free(tempRec);
    free(divSites);
    free(modSites);

    
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
    int (*R)[2];
    a=0;
    b=0;
    e=0;
    
    rlen=&e;
    length=&a;
    star=&b;
    c=numHind(0);
    R=malloc(c*2*(sizeof(int)));
    
   baseCase(length, R, star, rlen);
    
    f=*rlen;
    R=realloc(R,(f*2)*(sizeof(int)));
    printR(length,R,star,rlen,f);
    next=*star;
    c=numHind(next);
    while(next<70){
       c=numHind(next);
       //printf("c=%d rlen=%d\n",c,*rlen);
       if(c>*rlen){
          R=realloc(R,(c*2)*(sizeof(int)));
      }
      else{
          R=realloc(R,((*rlen)*2)*(sizeof(int)));   
      }     
       //printf("BEFORE: len=%d star=%d rlen=%d \n", *length, *star, *rlen);
       genCase(length, R, star, rlen);
       if(*length==0){
          printf("gencase\n");
          f=*rlen;
          printR(length,R,star,rlen,f); 
         
       } else{   
          printf("gencase\n");
          f=*rlen;
          //printf("f=%d\n", f);
          printR(length,R,star,rlen,f); 
          next=*star;
      }    
   }              
    free(R);
    system("PAUSE");
}

/*
 *	engwindemo.c
 *
 *	This is a simple program that illustrates how to call the MATLAB
 *	Engine functions from a C program for windows
 *
 * Copyright 1984-2003 The MathWorks, Inc.
 */
/* $Revision: 1.10.4.1 $ */
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine.h"

#define BUFSIZE 5000

FILE *error;

//static double Areal[6] = { 1, 2, 3, 4, 5, 6 };

int PASCAL WinMain (HINSTANCE hInstance,
                    HINSTANCE hPrevInstance,
                    LPSTR     lpszCmdLine,
                    int       nCmdShow)

{
	Engine *ep;
	//mxArray *result = NULL, *b = NULL;
	//mxArray *T = NULL, *a = NULL, *d = NULL;
	char buffer[BUFSIZE+1];
	//char output[5000];
//	mxArray *d = NULL;
//	double *dNum;
	
    //double *Dreal, *Dimag;
	//double time[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

	/*
	 * Start the MATLAB engine 
	 */
	if (!(ep = engOpen(NULL))) {
		MessageBox ((HWND)NULL, (LPSTR)"Can't start MATLAB engine", 
			(LPSTR) "Engwindemo.c", MB_OK);
		exit(-1);
	}
	buffer[BUFSIZE] = '\0';
	engOutputBuffer(ep, buffer, 5000);
 //engSetVisible(ep, 1);
	/*
	 * PART I
	 *
	 * For the first half of this demonstration, we will send data
	 * to MATLAB, analyze the data, and plot the result.
	 */

	/* 
	 * Create a variable from our data
	 */
//	T = mxCreateDoubleMatrix(1, 10, mxREAL);
//	memcpy((char *) mxGetPr(T), (char *) time, 10*sizeof(double));
     
     
     
	/*
	 * Place the variable T into the MATLAB workspace
	 */
//	engPutVariable(ep, "T", T);

	/*
	 * Evaluate a function of time, distance = (1/2)g.*t.^2
	 * (g is the acceleration due to gravity)
	 */
//	 engEvalString(ep, "D = .5.*(-9.8).*T.^2;");

/* *************************************************************************
    My code
    
	engEvalString(ep, "pVectCheck;");
    engEvalString(ep, "load sparseMatrixV1.txt");
	engEvalString(ep, "A = spconvert(sparseMatrixV1);");
	if ((result = engGetVariable(ep,"A")) == NULL)
	      printf("Oops! You didn't create a variable X.\n\n");
    system("PAUSE");
	engEvalString(ep, "[a] = textread('bVector.txt', '', 'delimiter', ',');");
	engEvalString(ep, "b = A\a;");
	if((b=engGetVariable(ep,"b")) == NULL){
           sprintf(buffer, "B is NULL");
           printf("Oops! No b calculated!");
    }
    printf("after b\n");
    system("PAUSE");
	engEvalString(ep, "save b.txt b -ascii;");
   **************************************************************************
*/
	
	/*
	 * Plot the result
	 */
	 error = fopen("error.txt", "w");
    if ((error = fopen("error.txt", "w"))) {
       int err = engEvalString(ep, "load sparseMatrixV1.txt;");
	   fprintf(error, "%d\n",err);
    }
    fclose(error);
    engEvalString(ep, "A=spconvert(sparseMatrixV1);");
    engEvalString(ep, "[a] = textread('bVector.txt', '', 'delimiter', ',');");
    engEvalString(ep, "b = A\\a;");
    engEvalString(ep, "save b.txt b -ascii;");
    engEvalString(ep,"pwd");
	/*engEvalString(ep, "title('Position vs. Time for a falling object');");
	engEvalString(ep, "xlabel('Time (seconds)');");
	engEvalString(ep, "ylabel('Position (meters)');");*/

    
    /*
	 * PART II
	 *
	 * For the second half of this demonstration, we will create another mxArray
	 * put it into MATLAB and calculate its eigen values 
	 * 
	 */
	  
//	 a = mxCreateDoubleMatrix(3, 2, mxREAL);         
//	 memcpy((char *) mxGetPr(a), (char *) Areal, 6*sizeof(double));
//	 engPutVariable(ep, "A", a); 

	 /*
	 * Calculate the eigen value
	 */
//	 engEvalString(ep, "d = eig(A*A')");

	 /*
	 * Use engOutputBuffer to capture MATLAB output. Ensure first that
	 * the buffer is always NULL terminated.
	 */
	 
	 //engOutputBuffer(ep, buffer, BUFSIZE);

	 /*
	 * the evaluate string returns the result into the
	 * output buffer.
	 */
//	 engEvalString(ep, "whos");
//	 MessageBox ((HWND)NULL, (LPSTR)buffer, (LPSTR) "MATLAB - whos", MB_OK);
	
	 /*
	 * Get the eigen value mxArray
	 */
	 //d = engGetVariable(ep, "d");
	 	MessageBox ((HWND)NULL, (LPSTR)buffer, (LPSTR)"Engwindemo.c", MB_OK);
	 engClose(ep);
    //sprintf(buffer,"Hello");

    /* sprintf(buffer,"Hello\n");
     if(b==NULL){
        sprintf(buffer,"b=NULL");
       // system("PAUSE");
     }else{
        Dreal = mxGetPr(b);
        sprintf(buffer,"b[0] : %g", Dreal[1]);
     }
   	 MessageBox ((HWND)NULL, (LPSTR)buffer, (LPSTR)"Engwindemo.c", MB_OK);
    
     mxDestroyArray(result);
     mxDestroyArray(b);	*/ 
	
    /* if (d == NULL) {
			MessageBox ((HWND)NULL, (LPSTR)"Get Array Failed", (LPSTR)"Engwindemo.c", MB_OK);
		}
	else {		
		Dreal = mxGetPr(d);
		Dimag = mxGetPi(d);      		
		if (Dimag)
			sprintf(buffer,"Eigenval 2: %g+%gi",Dreal[1],Dimag[1]);
		else
			sprintf(buffer,"Eigenval 2: %g",Dreal[1]);
		MessageBox ((HWND)NULL, (LPSTR)buffer, (LPSTR)"Engwindemo.c", MB_OK);
	    mxDestroyArray(d);
	} */

	
	 /* We're done! Free memory, close MATLAB engine and exit.
	 */
//	mxDestroyArray(T);
//	mxDestroyArray(a);
	
	return(0);
}

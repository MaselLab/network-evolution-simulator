#include <stdio.h>
#include <math.h>

/* program to test pointers to locations in kon array */

void print_matrix(float *matrix[4][4]) {
  int i, j;
  for (i=0; i < 4; i++) {
    for (j=0; j < 4; j++) {
      if (matrix[i][j] != NULL)
	printf("%.4g ", *matrix[i][j]);
      else
	printf("NULL ");
    }
    printf("\n");
  }
}

int main(int argc, char *argv[])
{
  int i, j;

  /* make each element be a NULL ptr by default */
  float *matrix[4][4] = {{NULL, NULL, NULL, NULL}, 
			 {NULL, NULL, NULL, NULL}, 
			 {NULL, NULL, NULL, NULL}, 
			 {NULL, NULL, NULL, NULL}};
  float kon[4] = {0.1, 0.4, 0.2, 0.3};

  /* generate matrix */
  for (i=0; i < 4; i++) {
    for(j=0; j < 4; j++) {
      /* put kon values along the diagonal (this is only for testing!!
	 in reality kon values won't be along diagonal) */
      if (i == j) {
	matrix[i][j] = &(kon[i]);
      }
    }
  }
  /* print out initial matrix */
  printf("initial matrix:\n");
  print_matrix(matrix);

  /* modify kon */
  kon[0] = 0.0789;
  kon[3] = 0.9999;

  printf("\nmodified matrix:\n");

  /* this shows the modified kon values */
  print_matrix(matrix);
}

/* 
   Local Variables:
   compile-command: "gcc -o element-ptr element-ptr.c -lm"
   End: 
*/



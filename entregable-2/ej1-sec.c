#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

double dwalltime();
void mult_x_col(double *dest,double *param1,double *param2);
void suma(double *dest, double *param1,double *param2,double *param3,double escalar);
void check(double *RES,double promedios);
void init_matrices(double *A,double *B, double *C, double *D);
void init_matrices_tri(double *U,double *L, int N);
void print_matrix(double *RES);
void mult_sup(double *dest,double *param1,double *sup);
void mult_low(double *dest,double *low,double *param1);
double promedio(double *dest,int N);

int N;

int main( int argc, char *argv[]){
	int i;
	double *A,*B,*C,*D,*L,*U,*AUX1,*AUX2,*AUX3,*AUX4;
 	double timetick,n_timetick;
 	double promedio_u, promedio_l;

	//Controla los argumentos al programa
 	if ((argc != 2) || ((N = atoi(argv[1])) <= 0) ){
		printf("\nUsar: %s n \n  n: Dimension de la matriz (nxn X nxn)\n ", argv[0]);
		exit(1);
	}

 	//Aloca memoria para las matrices
	double triangular_matrix_size = N*(N+1)/2;
	double square_matrix_size = N*N;

	//Aloca memoria para las matrices
	A=(double*)malloc(sizeof(double)*square_matrix_size);
	C=(double*)malloc(sizeof(double)*square_matrix_size);
	B=(double*)malloc(sizeof(double)*square_matrix_size);
	D=(double*)malloc(sizeof(double)*square_matrix_size);
	L=(double*)malloc(sizeof(double)*triangular_matrix_size);
	U=(double*)malloc(sizeof(double)*triangular_matrix_size);
	AUX1 = (double*)malloc(sizeof(double)*square_matrix_size);
	AUX2 = (double*)malloc(sizeof(double)*square_matrix_size);
	AUX3 = (double*)malloc(sizeof(double)*square_matrix_size);
	AUX4 = (double*)malloc(sizeof(double)*square_matrix_size);

	init_matrices(A,B,C,D);
	init_matrices_tri(U,L,triangular_matrix_size);
  	timetick = dwalltime();

	mult_x_col(AUX1,A,B);


	promedio_u = promedio(L,triangular_matrix_size);
	promedio_l = promedio(U,triangular_matrix_size);

	mult_low(AUX2,L,C);

	mult_sup(AUX3,D,U);

	double p = promedio_u*promedio_l;
	// printf("escalar = %f\n", p);
	suma(AUX4,AUX1,AUX2,AUX3,p);

	n_timetick = dwalltime() - timetick;
    printf("Tiempo de ejecución: %f\n", n_timetick);
    check(AUX4,p);
	free(AUX4);
	free(AUX3);
	free(AUX2);
	free(AUX1);
	free(A);
	free(B);
	free(C);
	free(D);
	free(L);
	free(U);
	return 0;
}

void init_matrices(double *A,double *B, double *C, double *D){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = 1;
			B[i*N+j] = 1;
			C[i*N+j] = 1;
			D[i*N+j] = 1;
		}
	}
}

void init_matrices_tri(double *U,double *L, int N){
	int i,j;
	for(i=0;i<N;i++){
	    U[i]=1;
		L[i]=1;	
	}
}

void mult_x_col(double *dest,double *param1,double *param2){
  int i,j,k;
  double sum;
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<N;k++){
        sum +=  param1[i*N+k] * param2[j*N+k];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void check(double *RES,double promedio){
	int i,j;
	char col_check = 1;
	char row_check = 1;
	double col_tmp,row_tmp;
	row_tmp = RES[0]/promedio;
	for(i=0;i<N;i++){
		col_tmp = (RES[i*N]/promedio);
		row_check = row_check&&row_tmp==(col_tmp - i);
		for(j=0;j<N;j++){
			col_check=col_check&&(RES[i*N+j]/promedio-j)==col_tmp;
		}
    }
	if(col_check && row_check){
	  printf("Multiplicacion de matrices resultado correcto\n");
	}else{
	  printf("Multiplicacion de matrices resultado erroneo\n");
	}
}

void print_matrix(double *RES){
	int i,j;
	for(i=0;i<N;i++){
     for(j=0;j<N;j++){
	  	printf("%f ",RES[i*N+j]);
     }
	 printf("\n");
    } 
}

void print_matrix_pid(double *RES,int pid){
	int i,j;
	for(i=0;i<N;i++){
     for(j=0;j<N;j++){
	  	printf(" %d: [%d,%d] %f",pid,i,j,RES[i*N+j]);
     }
	 printf("\n");
    } 
}

double dwalltime(){
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

void mult_sup(double *dest,double *param1,double *sup){
  int i,j,k;
  double sum;
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=j;k++){
        sum +=  param1[i*N+k] * sup[(row*(2*N - row + 1))/2 + col - row];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void mult_low(double *dest,double *low,double *param1){
  int i,j,k;
  double sum;
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=i;k++){
        sum +=  low[i+j*(j+1)/2] * param1[j*N+k];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void suma(double *dest, double *param1,double *param2,double *param3,double escalar){
  int i,j,index;
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      index = i*N+j; 
      dest[index] =  (param1[index] + param2[index] + param3[index]) * escalar;
    }
  }
}

double promedio(double *param, int size){
	int i;
	double sum=0;
	for(i=0;i<size;i++){
      sum += param[i];
	}
	return sum / (N*N);
}

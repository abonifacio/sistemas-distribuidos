#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define A_TAG 1
#define B_TAG 1
#define C_TAG 1

// int gettimeofday(struct timeval *tv,struct timezone *tz);
double dwalltime();
void transpuesta(double *dest, double *param);
void mult_x_col(double *dest,double *param1,double *param2);
void check(double *RES);
void init_matrices(double *A, double *B);
 	
int CANT_FILAS,N,NUM_PROCESOS;


int main( int argc, char *argv[]){
	int i,myrank;
	MPI_Status status;
	double *A,*B,*Bt,*C;
 	double timetick,n_timetick;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
	MPI_Comm_size(MPI_COMM_WORLD,&NUM_PROCESOS);

	//Controla los argumentos al programa
 	if ((argc != 2) || ((N = atoi(argv[1])) <= 0) ){
		printf("\nUsar: %s n \n  n: Dimension de la matriz (nxn X nxn)\n ", argv[0]);
		exit(1);
	}

	CANT_FILAS = N / NUM_PROCESOS;
 	//Aloca memoria para las matrices
	Bt=(double*)malloc(sizeof(double)*N*N);

  	if (myrank == 0) {

		A=(double*)malloc(sizeof(double)*N*N);
	  	B=(double*)malloc(sizeof(double)*N*N);
		init_matrices(A,B);	  	
	  
  		timetick = dwalltime();
  		transpuesta(Bt,B);
		free(B);
		// mandar la matriz a los procesos
		C=(double*)malloc(sizeof(double)*N*N);
		for(i=1;i<NUM_PROCESOS;i++){
			MPI_Send(&A[(i+CANT_FILAS-1)*N], CANT_FILAS*N, MPI_DOUBLE, i, A_TAG, MPI_COMM_WORLD);
			MPI_Send(Bt, N*N, MPI_DOUBLE, i, B_TAG, MPI_COMM_WORLD);
		}
	}else{
		// recivo la matriz A y Bt
	  	A=(double*)malloc(sizeof(double)*CANT_FILAS*N);
  		C=(double*)malloc(sizeof(double)*CANT_FILAS*N);
	  	MPI_Recv(A, CANT_FILAS*N, MPI_DOUBLE, 0, A_TAG, MPI_COMM_WORLD, &status);
	  	MPI_Recv(Bt, N*N, MPI_DOUBLE, 0, B_TAG, MPI_COMM_WORLD, &status);
	}

	mult_x_col(C,A,Bt);

	if (myrank == 0) {
		// recibo en RES los resultados Ci
		for(i=1;i<NUM_PROCESOS;i++){
	  		MPI_Recv(&C[(i+CANT_FILAS-1)*N], CANT_FILAS*N, MPI_DOUBLE, i, C_TAG, MPI_COMM_WORLD, &status);
		}

  		n_timetick = dwalltime() - timetick;
	    printf("Tiempo de ejecución: %f\n", n_timetick);
	    check(C);
	}else{
		MPI_Send(C, CANT_FILAS*N, MPI_DOUBLE, 0, C_TAG, MPI_COMM_WORLD);
	}
	free(C);
	free(A);
	free(Bt);
	MPI_Finalize();
	return 0;
}

void init_matrices(double *A, double *B){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = 1;
			B[i*N+j] = 1;
		}
	}
}

void transpuesta(double *dest, double *param){
  int i,j;
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      dest[i*N+j] = param[j*N+i];
    }
  }
}

void mult_x_col(double *dest,double *param1,double *param2){
  int i,j,k;
  double sum;
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<N;k++){
        sum +=  param1[i*N+k] * param2[j*N+k];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void check(double *RES){
	char check = 1;
	int i,j;
	for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      	check=check&&(RES[i*N+j]==N);
     }
    } 
	if(check){
	  printf("Multiplicacion de matrices resultado correcto\n");
	}else{
	  printf("Multiplicacion de matrices resultado erroneo\n");
	}
}

double dwalltime(){
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}
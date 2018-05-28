#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define ROOT_PID 0

// int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)

// int gettimeofday(struct timeval *tv,struct timezone *tz);
double dwalltime();
void transpuesta(double *dest, double *param);
void mult_x_col(double *dest,double *param1,double *param2);
void check(double *RES);
void init_matrices(double *A, double *B);
void print_matrix(double *RES);
void print_matrix_pid(double *RES,int pid);
 	
int CANT_FILAS,N,NUM_PROCESOS;


int main( int argc, char *argv[]){
	int i,myrank;
	MPI_Status status;
	double *A,*A_SRC,*B,*Bt,*C,*RES;
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
  	C=(double*)malloc(sizeof(double)*N*CANT_FILAS);
	A=(double*)malloc(sizeof(double)*N*CANT_FILAS);

  	if (myrank == ROOT_PID) {
	  	A_SRC=(double*)malloc(sizeof(double)*N*N);
  		RES=(double*)malloc(sizeof(double)*N*N);
	  	B=(double*)malloc(sizeof(double)*N*N);
		init_matrices(A_SRC,B);	  
  		timetick = dwalltime();
  		transpuesta(Bt,B);
		free(B);
	}
	MPI_Scatter(A_SRC, (N*N)/NUM_PROCESOS, MPI_DOUBLE, A,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Bcast(Bt, N*N, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	mult_x_col(C,A,Bt);
	MPI_Gather(C, (N*N)/NUM_PROCESOS, MPI_DOUBLE, RES,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);


	if (myrank == ROOT_PID) {
  		n_timetick = dwalltime() - timetick;
	    printf("Tiempo de ejecución: %f\n", n_timetick);
	    check(RES);
	    free(RES);
	    free(A_SRC);
	}
	free(A);
	free(C);
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

void print_matrix(double *RES){
	int i,j;
	for(i=0;i<N;i++){
     for(j=0;j<N;j++){
	  	printf("%f \n",RES[i*N+j]);
     }
	 printf("\n");
    } 
}

void print_matrix_pid(double *RES,int pid){
	int i,j;
	for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
	  	printf("[%d] %f",pid,RES[i*N+j]);
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
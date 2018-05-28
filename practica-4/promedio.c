#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define ROOT_PID 0
#define RESULTADO 4

// int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)

// int gettimeofday(struct timeval *tv,struct timezone *tz);
double dwalltime();
double suma(double *param);
double promedio(double *RES);
void check(double RES);
void init_matrices(double *A);
void print_matrix(double *RES);
void print_matrix_pid(double *RES,int pid);
 	
int CANT_FILAS,N,NUM_PROCESOS;


int main( int argc, char *argv[]){
	int i,myrank;
	MPI_Status status;
	double *A,*A_SRC,*RES,tmp;
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
	A=(double*)malloc(sizeof(double)*N*CANT_FILAS);

  	if (myrank == ROOT_PID) {
	  	A_SRC=(double*)malloc(sizeof(double)*N*N);
	  	RES=(double*)malloc(sizeof(double)*NUM_PROCESOS);
		init_matrices(A_SRC);	  
  		timetick = dwalltime();
	}
	MPI_Scatter(A_SRC, (N*N)/NUM_PROCESOS, MPI_DOUBLE, A,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);	
	tmp = suma(A);
	MPI_Gather(&tmp, 1, MPI_DOUBLE, RES, 1, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);


	if (myrank == ROOT_PID) {
		tmp = promedio(RES);
	    free(RES);
  		n_timetick = dwalltime() - timetick;
	    printf("Tiempo de ejecución: %f\n", n_timetick);
	    check(tmp);
	    free(A_SRC);
	}
	free(A);
	MPI_Finalize();
	return 0;
}

void init_matrices(double *A){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = RESULTADO;
		}
	}
}

double promedio(double *RES){
	int i;
	double sum = 0;
	for(i=0;i<NUM_PROCESOS;i++){
		sum+=RES[i];
	}
	return sum / (N*N);
}

double suma(double *RES){
	int i,j;
	double sum = 0;
	for(i=0;i<CANT_FILAS;i++){
		for(j=0;j<N;j++){
			sum+=RES[i*N+j];
		}
	}
	return sum;
}

void check(double RES){
	if(RES==RESULTADO){
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
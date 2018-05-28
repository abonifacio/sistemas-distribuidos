#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define ROOT_PID 0

// int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)

// int gettimeofday(struct timeval *tv,struct timezone *tz);
double dwalltime();
void mult_x_col(double *dest,double *param1,double *param2);
void suma(double *dest, double *param1,double *param2,double *param3);
void check(double *RES);
void init_matrices(double *A,double *B, double *L, double *C, double *D,double *U);
void print_matrix(double *RES);
void print_matrix_pid(double *RES,int pid);
void mult_sup(double *dest,double *param1,double *sup);
void mult_low(double *dest,double *param1,double *low);

int CANT_FILAS,N,NUM_PROCESOS;


int main( int argc, char *argv[]){
	int i,myrank,matrix_size;
	MPI_Status status;
	double *A,*B,*A_B,*L,*C,*L_C,*D_U,*U,*D,*RES;
 	double timetick,n_timetick;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
	MPI_Comm_size(MPI_COMM_WORLD,&NUM_PROCESOS);

	//Controla los argumentos al programa
 	if ((argc != 2) || ((N = atoi(argv[1])) <= 0) ){
		printf("\nUsar: %s n \n  n: Dimension de la matriz (nxn X nxn)\n ", argv[0]);
		exit(1);
	}

	if (myrank == ROOT_PID) {
  		matrix_size = N*N;
  	}else{
  		matrix_size = N*N;
  	}

	CANT_FILAS = N / NUM_PROCESOS;
 	//Aloca memoria para las matrices
	A = (double*)malloc(sizeof(double)*matrix_size);
	B = (double*)malloc(sizeof(double)*N*N);
	L = (double*)malloc(sizeof(double)*matrix_size);
	C = (double*)malloc(sizeof(double)*N*N);
	D = (double*)malloc(sizeof(double)*matrix_size);
	U = (double*)malloc(sizeof(double)*N*N);

  	if (myrank == ROOT_PID) {
		init_matrices(A,B,L,C,D,U);  
  		timetick = dwalltime();
	}
	MPI_Bcast(B, N*N, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Bcast(C, N*N, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Bcast(U, N*N, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Scatter(A, (N*N)/NUM_PROCESOS, MPI_DOUBLE, A,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Scatter(L, (N*N)/NUM_PROCESOS, MPI_DOUBLE, L,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Scatter(D, (N*N)/NUM_PROCESOS, MPI_DOUBLE, D,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	printf("A(%d)\n", myrank);
	print_matrix_pid(A,myrank);
	printf("B(%d)\n", myrank);
	print_matrix_pid(B,myrank);
	A_B = (double*)malloc(sizeof(double)*matrix_size);
	mult_x_col(A_B,A,B);
	// print_matrix(A_B);
	free(A);
	free(B);
	
	L_C = (double*)malloc(sizeof(double)*matrix_size);
	mult_low(L_C,C,L);
	free(C);
	free(L);
	
	D_U = (double*)malloc(sizeof(double)*matrix_size);
	mult_sup(D_U,D,U);
	free(D);
	free(U);
	
	MPI_Gather(A_B, (N*N)/NUM_PROCESOS, MPI_DOUBLE, A_B,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Gather(L_C, (N*N)/NUM_PROCESOS, MPI_DOUBLE, L_C,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Gather(D_U, (N*N)/NUM_PROCESOS, MPI_DOUBLE, D_U,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	
	RES = (double*)malloc(sizeof(double)*matrix_size);
	suma(RES,A_B,L_C,D_U);
	free(A_B);
	free(L_C);
	free(D_U);

	MPI_Gather(RES, (N*N)/NUM_PROCESOS, MPI_DOUBLE, RES,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);

	if (myrank == ROOT_PID) {
  		n_timetick = dwalltime() - timetick;
	    printf("Tiempo de ejecución: %f\n", n_timetick);
	    check(RES);
	}
	free(RES);
	MPI_Finalize();
	return 0;
}

void init_matrices(double *A,double *B, double *L, double *C, double *D,double *U){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = 1;
			B[i*N+j] = 1;
			C[i*N+j] = 1;
			D[i*N+j] = 1;
			L[i*N+j]=0;
		    U[i*N+j]=0;
		    if(i<=j){
		      U[i*N+j]=1;
		    }
		    if(i>=j){
		      L[i*N+j]=1;
		    }
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
	int i,j;
	char check = 1;
	double tmp;
	for(i=0;i<N;i++){
		tmp = RES[i*N];
		for(j=0;j<N;j++){
	  		printf("%f ",RES[i*N+j]);
			check=check&&(RES[i*N+j]==tmp);
		}
		printf("\n");
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
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=i;k<N;k++){
        sum +=  param1[i*N+k] * sup[k*N+j];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void mult_low(double *dest,double *param1,double *low){
  int i,j,k;
  double sum;
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=i;k++){
        sum +=  param1[i*N+k] * low[k*N+j];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void suma(double *dest, double *param1,double *param2,double *param3){
  int i,j,index;
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      index = i*N+j; 
      dest[index] =  param1[index] + param2[index] + param3[index];
    }
  }
}
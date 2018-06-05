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
void suma2(double *dest, double *param1,double *param2);
void check(double *RES);
void init_matrices(double *A,double *B, double *L, double *C, double *D);
void init_matriz_u(double *U, int N);
void print_matrix(double *RES);
void print_matrix_pid(double *RES,int pid);
int get_index(int row, int col);
void mult_sup(double *dest,double *param1,double *sup);
void print_tri_rank(double *RES);
void mult_low(double *dest,double *low,double *param1);

int CANT_FILAS,N,NUM_PROCESOS,myrank;


int main( int argc, char *argv[]){
	int i,matrix_size,upper_matrix_size;
	MPI_Status status;
	double *A,*B,*A_B,*L,*C,*L_C,*D_U,*U,*D,*M,*RES;
	double *A_SRC,*L_SRC,*D_SRC;
 	double timetick,n_timetick;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
	MPI_Comm_size(MPI_COMM_WORLD,&NUM_PROCESOS);

	//Controla los argumentos al programa
 	if ((argc != 2) || ((N = atoi(argv[1])) <= 0) ){
		printf("\nUsar: %s n \n  n: Dimension de la matriz (nxn X nxn)\n ", argv[0]);
		exit(1);
	}



	// if (myrank == ROOT_PID) {
 //  		matrix_size = N*N;
 //  	}else{
	// 	matrix_size = (N*N)/NUM_PROCESOS;
 //  	}

	CANT_FILAS = N / NUM_PROCESOS;
	matrix_size = CANT_FILAS*N;
	upper_matrix_size = N*(N+1)/2;
 	//Aloca memoria para las matrices
	A = (double*)malloc(sizeof(double)*matrix_size);
	B = (double*)malloc(sizeof(double)*N*N);
	L = (double*)malloc(sizeof(double)*matrix_size);
	C = (double*)malloc(sizeof(double)*N*N);
	D = (double*)malloc(sizeof(double)*matrix_size);
	U = (double*)malloc(sizeof(double)*upper_matrix_size);

  	if (myrank == ROOT_PID) {
		A_SRC = (double*)malloc(sizeof(double)*N*N);
		L_SRC = (double*)malloc(sizeof(double)*N*N);
		D_SRC = (double*)malloc(sizeof(double)*N*N);
		init_matrices(A_SRC,B,L_SRC,C,D_SRC);
		init_matriz_u(U,upper_matrix_size);
  		timetick = dwalltime();
	}
	MPI_Bcast(B, N*N, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Bcast(C, N*N, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Bcast(U, upper_matrix_size, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Scatter(A_SRC, (N*N)/NUM_PROCESOS, MPI_DOUBLE, A,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Scatter(L_SRC, (N*N)/NUM_PROCESOS, MPI_DOUBLE, L,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Scatter(D_SRC, (N*N)/NUM_PROCESOS, MPI_DOUBLE, D,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);

	if (myrank == ROOT_PID) {
		free(A_SRC);
		free(L_SRC);
		free(D_SRC);
	}
	A_B = (double*)malloc(sizeof(double)*matrix_size);
	mult_x_col(A_B,A,B);
	free(A);
	free(B);
	
	L_C = (double*)malloc(sizeof(double)*matrix_size);
	mult_low(L_C,L,C);
	free(C);
	free(L);

	D_U = (double*)malloc(sizeof(double)*matrix_size);
	mult_sup(D_U,D,U);
	free(D);
	free(U);
	
	M = (double*) malloc(sizeof(double)*matrix_size);
	suma(M,A_B,L_C,D_U);
	free(A_B);
	free(L_C);
	free(D_U);

	if(myrank==ROOT_PID){
		RES = (double*)malloc(sizeof(double)*N*N);
	}
	MPI_Gather(M, (N*N)/NUM_PROCESOS, MPI_DOUBLE,RES,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	
	if (myrank == ROOT_PID) {
  		n_timetick = dwalltime() - timetick;
	    printf("Tiempo de ejecución: %f\n", n_timetick);
	    check(RES);
	    // print_matrix(RES);
		free(RES);
	}
	free(M);
	MPI_Finalize();
	return 0;
}

void init_matrices(double *A,double *B, double *L, double *C, double *D){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = 1;
			B[i*N+j] = 1;
			C[i*N+j] = 1;
			D[i*N+j] = 1;
			L[i*N+j]=0;
		    if(i>=j){
		      L[i*N+j]=1;
		    }
		}
	}
}

void init_matriz_u(double *U, int N){
	int i,j;
	for(i=0;i<N;i++){
	    U[i]=1;
	}
}

void print_tri_rank(double *RES){
	int i,j,inicio = myrank*CANT_FILAS, fin = (myrank + 1)*CANT_FILAS;
	for(i=inicio;i<fin;i++){
     for(j=i;j<N;j++){
	  	printf("(%d) %f ",myrank, RES[get_index(i,j)]);
     }
	 printf("\n");
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
	char col_check = 1;
	char row_check = 1;
	double col_tmp,row_tmp;
	row_tmp = RES[0];
	for(i=0;i<N;i++){
		col_tmp = RES[i*N];
		row_check = row_check&&row_tmp==(col_tmp - i);
		for(j=0;j<N;j++){
			col_check=col_check&&(RES[i*N+j]-j)==col_tmp;
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
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=j;k++){
        sum +=  param1[i*N+k] * sup[get_index(k,j)];
      }
      dest[i*N+j] =  sum;
    }
  }
}

int get_index(int row, int col){
	return (row*(2*N - row + 1))/2 + col - row;
}

void mult_low(double *dest,double *low,double *param1){
  int i,j,k,inicio;
  double sum;
  inicio = myrank*CANT_FILAS;
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=(inicio+i);k++){
        sum +=  low[i*N+k] * param1[k*N+j];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void suma2(double *dest, double *param1,double *param2){
  int i,j,index;
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      index = i*N+j; 
      dest[index] =  param1[index] + param2[index];
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
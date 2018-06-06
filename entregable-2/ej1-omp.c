#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <omp.h>

#define ROOT_PID 0

// int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)

// int gettimeofday(struct timeval *tv,struct timezone *tz);
double dwalltime();
void mult_x_col(double *dest,double *param1,double *param2);
void suma(double *dest, double *param1,double *param2,double *param3,double escalar);
void check(double *RES,double promedios);
void init_matrices(double *A,double *B, double *L, double *C, double *D);
void init_matriz_u(double *U, int N);
void print_matrix(double *RES);
void print_matrix_pid(double *RES,int pid);
int get_index(int row, int col);
void mult_sup(double *dest,double *param1,double *sup);
void print_tri_rank(double *RES);
void mult_low(double *dest,double *low,double *param1);
double sumar_low(double *dest);
double sumar_sup(double *dest,int N);

int CANT_FILAS,N,NUM_PROCESOS,myrank,NUM_THREADS;


int main( int argc, char *argv[]){
	int i,matrix_size,upper_matrix_size;
	MPI_Status status;
	double *A,*B,*A_B,*L,*C,*L_C,*D_U,*U,*D,*M,*RES;
	double *A_SRC,*L_SRC,*D_SRC;
 	double timetick,n_timetick;
 	double *promedios,*promedios_res;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
	MPI_Comm_size(MPI_COMM_WORLD,&NUM_PROCESOS);

	//Controla los argumentos al programa
 	 if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((NUM_THREADS = atoi(argv[2])) <= 0) ){
	    printf("\nUsar: %s n t\n  n: Dimension de la matriz (nxn X nxn) t: Cantidad de threads \n", argv[0]);
	    exit(1);
	 }
	 omp_set_num_threads(NUM_THREADS);



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
	promedios = (double*)malloc(sizeof(double)*2);
	promedios_res = (double*)malloc(sizeof(double)*2);

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


	for(i=0;i<NUM_THREADS;i++){
	    printf("Proceso %d thread %d\n", myrank,i);
	}

	
	A_B = (double*)malloc(sizeof(double)*matrix_size);
	mult_x_col(A_B,A,B);
	free(A);
	free(B);

	promedios[0] = sumar_low(L) / (N*N);
	promedios[1] = sumar_sup(U,upper_matrix_size) / (N*N);

	// printf("-------------- Promedio L %f \n", promedios[0]);
	// printf("-------------- Promedio U %f \n", promedios[1]);
	
	L_C = (double*)malloc(sizeof(double)*matrix_size);
	mult_low(L_C,L,C);
	free(C);
	free(L);

	D_U = (double*)malloc(sizeof(double)*matrix_size);
	mult_sup(D_U,D,U);
	free(D);
	free(U);

	MPI_Allreduce(promedios, promedios_res, 2, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
	double p = promedios_res[0]*promedios_res[1];
	// printf("escalar = %f\n", p);
	M = (double*) malloc(sizeof(double)*matrix_size);
	suma(M,A_B,L_C,D_U,p);
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
	    check(RES,p);
	    //print_matrix(RES);
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
  #pragma omp parallel for schedule(static) private(i,j,k,sum)
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

double sumar_low(double *dest){
	int i,j,inicio;
	double sum = 0;
	inicio = myrank*CANT_FILAS;
	#pragma omp parallel for schedule(static) private(i,j,inicio) collapse(2)
		for(i=0;i<CANT_FILAS;i++){
		 	for(j=0;j<=(inicio+i);j++){
				sum += dest[i*N+j];
			}
		}
	return sum;
}

double sumar_sup(double *dest,int N){
	int i,j,desde,hasta,CE;
	CE = N/NUM_PROCESOS;
	desde = myrank * CE;
	hasta = (myrank + 1) * CE;
	if((N-hasta)<CE){
		hasta = N;
	}
	double sum = 0;
	#pragma omp parallel for schedule(static) private(i)
	for(i=desde;i<hasta;i++){
		sum += dest[i];
	}
	return sum;
}

void mult_sup(double *dest,double *param1,double *sup){
  int i,j,k;
  double sum;
  #pragma omp parallel for schedule(static) private(i,j,k,sum) collapse(2)
  	for(i=0;i<CANT_FILAS;i++){
  	   for(j=0;j<N;j++){
  	    sum=0;
  	    for(k=0;k<=j;k++){
  	      sum +=  param1[i*N+k] * sup[get_index(j,k)];
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
  #pragma omp parallel for schedule(static) private(i,j,k,sum) collapse (2)
  	for(i=0;i<CANT_FILAS;i++){
  	   for(j=0;j<N;j++){
  	    sum=0;
  	    for(k=0;k<=(inicio+i);k++){
  	      sum +=  low[i*N+k] * param1[j*N+k];
  	    }
  	    dest[i*N+j] =  sum;
  	  }
  	}
}

void suma(double *dest, double *param1,double *param2,double *param3,double escalar){
  int i,j,index;
  #pragma omp parallel for schedule(static) private(i,j,index) collapse (2)
	  for(i=0;i<CANT_FILAS;i++){
	     for(j=0;j<N;j++){
	      index = i*N+j; 
	      dest[index] =  (param1[index] + param2[index] + param3[index]) * escalar;
	    }
	  }
}

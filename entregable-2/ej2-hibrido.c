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
void mult_sup(double *dest,double *param1,double *sup);
void mult_low(double *dest,double *low,double *param1);
double promedio(double *dest, int N);
double promedio_low(double *dest);

int CANT_FILAS,N,NUM_PROCESOS,myrank,NUM_THREADS;


int main( int argc, char *argv[]){
	int i,matrix_size,upper_matrix_size;
	MPI_Status status;
	double *A,*B,*A_B,*L,*C,*L_C,*D_U,*U,*D,*M,*RES;
	double *A_SRC,*L_SRC,*D_SRC;
 	double timetick_total,timetick_comunicacion,timetick_procesamiento,tiempo_procesamiento = 0,tiempo_comunicacion = 0;
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

	CANT_FILAS = N / NUM_PROCESOS;
	matrix_size = CANT_FILAS*N;
	upper_matrix_size = N*(N+1)/2;
 	//Aloca memoria para las matrices
	A = (double*)malloc(sizeof(double)*matrix_size);
	A_B = (double*)malloc(sizeof(double)*matrix_size);
	L = (double*)malloc(sizeof(double)*matrix_size);
	B = (double*)malloc(sizeof(double)*N*N);
	C = (double*)malloc(sizeof(double)*N*N);
	D = (double*)malloc(sizeof(double)*matrix_size);
	U = (double*)malloc(sizeof(double)*upper_matrix_size);
	L_C = (double*)malloc(sizeof(double)*matrix_size);
	D_U = (double*)malloc(sizeof(double)*matrix_size);
	M = (double*) malloc(sizeof(double)*matrix_size);
	promedios = (double*)malloc(sizeof(double)*2);
	promedios_res = (double*)malloc(sizeof(double)*2);

  	if (myrank == ROOT_PID) {
		A_SRC = (double*)malloc(sizeof(double)*N*N);
		L_SRC = (double*)malloc(sizeof(double)*N*N);
		D_SRC = (double*)malloc(sizeof(double)*N*N);
		RES = (double*)malloc(sizeof(double)*N*N);
		init_matrices(A_SRC,B,L_SRC,C,D_SRC);
		init_matriz_u(U,upper_matrix_size);
  		timetick_comunicacion = dwalltime();
  		timetick_total = dwalltime();
	}
	MPI_Bcast(B, N*N, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Bcast(C, N*N, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Bcast(U, upper_matrix_size, MPI_DOUBLE, ROOT_PID, MPI_COMM_WORLD);
	MPI_Scatter(A_SRC, (N*N)/NUM_PROCESOS, MPI_DOUBLE, A,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Scatter(L_SRC, (N*N)/NUM_PROCESOS, MPI_DOUBLE, L,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	MPI_Scatter(D_SRC, (N*N)/NUM_PROCESOS, MPI_DOUBLE, D,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);

	if (myrank == ROOT_PID) {
  		tiempo_comunicacion += dwalltime() - timetick_comunicacion;
	}
	timetick_procesamiento = dwalltime();
	mult_x_col(A_B,A,B);

	promedios[0] = promedio_low(L);
	promedios[1] = promedio(U,upper_matrix_size);

	mult_low(L_C,L,C);

	mult_sup(D_U,D,U);
	if(myrank == ROOT_PID){
		timetick_comunicacion = dwalltime();
	}
	tiempo_procesamiento += dwalltime() - timetick_procesamiento;
	MPI_Allreduce(promedios, promedios_res, 2, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
	if(myrank == ROOT_PID){
  		tiempo_comunicacion += dwalltime() - timetick_comunicacion;
	}
	timetick_procesamiento = dwalltime();
	double p = promedios_res[0]*promedios_res[1];
	suma(M,A_B,L_C,D_U,p);

	tiempo_procesamiento += dwalltime() - timetick_procesamiento;
	if(myrank==ROOT_PID){
		timetick_comunicacion = dwalltime();
	}
	MPI_Gather(M, (N*N)/NUM_PROCESOS, MPI_DOUBLE,RES,(N*N)/NUM_PROCESOS, MPI_DOUBLE,ROOT_PID,MPI_COMM_WORLD);
	
	printf("Proceso #%d => Tiempo de procesamiento: %f \n", myrank, tiempo_procesamiento);
	if (myrank == ROOT_PID) {
  		tiempo_comunicacion += dwalltime() - timetick_comunicacion;
	    printf("Tiempo de comunicacion: %f \n", tiempo_comunicacion);
		printf("Tiempo total ROOT: %f \n", (dwalltime() - timetick_total));
	    check(RES,p);
	    // print_matrix(RES);
		free(RES);
		free(A_SRC);
		free(L_SRC);
		free(D_SRC);
	}
	free(C);
	free(L);
	free(D);
	free(U);
	free(A);
	free(B);
	free(A_B);
	free(L_C);
	free(D_U);
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

void mult_x_col(double *dest,double *param1,double *param2){
  int i,j,k;
  double sum;
  #pragma omp parallel for private(i,j,k,sum)
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

double dwalltime(){
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

double promedio(double *dest,int SIZE){
	int i,j,desde,hasta,CE;
	CE = SIZE/NUM_PROCESOS;
	desde = myrank * CE;
	hasta = (myrank + 1) * CE;
	if((SIZE-hasta)<CE){
		hasta = SIZE;
	}
	double sum = 0;
	#pragma omp parallel for private(i)
	for(i=desde;i<hasta;i++){
		sum += dest[i];
	}
	return sum / (N*N);
}

double promedio_low(double *dest){
	int i,j,inicio;
	double sum = 0;
	inicio = myrank*CANT_FILAS;
	#pragma omp parallel for private(i,j)
	for(i=0;i<CANT_FILAS;i++){
	 	for(j=0;j<=(inicio+i);j++){
			sum += dest[i*N+j];
		}
	}
	return sum / (N*N);
}

void mult_sup(double *dest,double *param1,double *sup){
  int i,j,k,inicio;
  double sum;
  inicio = myrank*NUM_PROCESOS;
  #pragma omp parallel for private(i,j,k,sum)
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=j;k++){
        sum +=  param1[i*N+k] * sup[(j*(2*N - j + 1))/2 + k - j];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void mult_low(double *dest,double *low,double *param1){
  int i,j,k,inicio;
  double sum;
  inicio = myrank*CANT_FILAS;
  #pragma omp parallel for private(i,j,k,sum)
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
  #pragma omp parallel for private(i,j,index) collapse (2)
  for(i=0;i<CANT_FILAS;i++){
     for(j=0;j<N;j++){
      index = i*N+j; 
      dest[index] =  (param1[index] + param2[index] + param3[index]) * escalar;
    }
  }
}

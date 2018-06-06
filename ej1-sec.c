#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

double dwalltime();
void mult_x_col(double *dest,double *param1,double *param2);
void suma(double *dest, double *param1,double *param2,double *param3,double escalar);
void suma2(double *dest, double *param1,double *param2);
void check(double *RES,double promedios);
void init_matrices(double *A,double *B, double *L, double *C, double *D);
void init_matriz_u(double *U, int N);
void print_matrix(double *RES);
void print_matrix_pid(double *RES,int pid);
int get_index(int row, int col);
void mult_sup(double *dest,double *param1,double *sup);
void mult_low(double *dest,double *low,double *param1);
double sumar_low(double *dest);
double sumar_sup(double *dest,int N);

int N;


int main( int argc, char *argv[]){
	int i,upper_matrix_size;
	double *A,*B,*A_B,*L,*C,*L_C,*D_U,*U,*D,*M;
 	double timetick,n_timetick;
 	double promedio_u, promedio_l;

	//Controla los argumentos al programa
 	if ((argc != 2) || ((N = atoi(argv[1])) <= 0) ){
		printf("\nUsar: %s n \n  n: Dimension de la matriz (nxn X nxn)\n ", argv[0]);
		exit(1);
	}

	upper_matrix_size = N*(N+1)/2;
 	//Aloca memoria para las matrices
	A = (double*)malloc(sizeof(double)*N*N);
	B = (double*)malloc(sizeof(double)*N*N);
	L = (double*)malloc(sizeof(double)*N*N);
	C = (double*)malloc(sizeof(double)*N*N);
	D = (double*)malloc(sizeof(double)*N*N);
	U = (double*)malloc(sizeof(double)*upper_matrix_size);

	init_matrices(A,B,L,C,D);
	init_matriz_u(U,upper_matrix_size);
  	timetick = dwalltime();

	A_B = (double*)malloc(sizeof(double)*N*N);
	mult_x_col(A_B,A,B);
	free(A);
	free(B);

	promedio_u = sumar_low(L) / (N*N);
	promedio_l = sumar_sup(U,upper_matrix_size) / (N*N);

	L_C = (double*)malloc(sizeof(double)*N*N);
	mult_low(L_C,L,C);
	free(C);
	free(L);

	D_U = (double*)malloc(sizeof(double)*N*N);
	mult_sup(D_U,D,U);
	free(D);
	free(U);

	double p = promedio_u*promedio_l;
	// printf("escalar = %f\n", p);
	M = (double*) malloc(sizeof(double)*N*N);
	suma(M,A_B,L_C,D_U,p);
	free(A_B);
	free(L_C);
	free(D_U);

	n_timetick = dwalltime() - timetick;
    printf("Tiempo de ejecución: %f\n", n_timetick);
    check(M,p);
	free(M);
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

double sumar_low(double *dest){
	int i,j;
	double sum = 0;
	for(i=0;i<N;i++){
	 	for(j=0;j<=i;j++){
			sum += dest[i*N+j];
		}
	}
	return sum;
}

double sumar_sup(double *dest,int N){
	int i,j;
	double sum = 0;
	for(i=0;i<N;i++){
		sum += dest[i];
	}
	return sum;
}

void mult_sup(double *dest,double *param1,double *sup){
  int i,j,k;
  double sum;
  for(i=0;i<N;i++){
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
  int i,j,k;
  double sum;
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=i;k++){
        sum +=  low[i*N+k] * param1[j*N+k];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void suma2(double *dest, double *param1,double *param2){
  int i,j,index;
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      index = i*N+j; 
      dest[index] =  param1[index] + param2[index];
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

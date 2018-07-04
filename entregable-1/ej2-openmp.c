#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<sys/time.h>

//Dimension por defecto de las matrices
int N=100;
int NUM_THREADS = 1;

//Para calcular tiempo
double dwalltime(){
  double sec;
  struct timeval tv;

  gettimeofday(&tv,NULL);
  sec = tv.tv_sec + tv.tv_usec/1000000.0;
  return sec;
}


void mult_sup(double *dest,double *sup,double *param1,double escalar){
  int i,j,k;
  double sum;
  #pragma omp parallel for private(i,j,k,sum)
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=i;k++){
        sum += sup[k*N+j-k*(k+1)/2] * param1[j*N+k];
      }
      dest[j*N+i] =  sum * escalar;
    }
  }
}

void mult_low(double *dest,double *low,double *param1,double escalar){
  int i,j,k;
  double sum;
  #pragma omp parallel for private(i,j,k,sum)
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=i;k++){
        sum +=  low[i+j*(j+1)/2] * param1[j*N+k];
      }
      dest[i*N+j] =  sum * escalar;
    }
  }
}


void mult_x_col(double *dest,double *param1,double *param2){
  int i,j,k;
  double sum;
  #pragma omp parallel for private(i,j,k,sum)
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

void transpuesta_escalar(double *dest, double *param,double escalar){
  int i,j;
  #pragma omp parallel for schedule(static) private(i,j)
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      dest[i*N+j] = param[j*N+i] * escalar;
    }
  }
}

void transpuesta(double *dest, double *param){
  int i,j;
  #pragma omp parallel for schedule(static) private(i,j)
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      dest[i*N+j] = param[j*N+i];
    }
  }
}

void suma(double *dest,double *param1,double *param2,double *param3){
	int i,j,index;
  #pragma omp parallel for schedule(static) private(i,j,index)
	for(i=0;i<N;i++){
	   for(j=0;j<N;j++){
	   	index = i*N+j; 
	   	dest[index] =  param1[index] + param2[index] + param3[index];
	  }
	}
}

double promedio(double *param, int size){
  int i;
  double sum=0;
  #pragma omp parallel for private(i) reduction( + : sum )
  for(i=0;i<size;i++){
      sum += param[i];
  }
  return sum / (N*N);
}

void print_matrix(double *M){
  int i,j;
  for (i = 0; i < N; ++i){
    for (j = 0; j < N; ++j){
      printf("M[%d,%d] = %f\n",i,j, M[i*N+j]);
    }
  }
}

int main(int argc,char*argv[]){
 double *A,*At,*B,*C,*D,*E,*F,*L,*U,*AUX1,*AUX2,*AUX3,*AUX4;
 int i,j,k;
 int check = 1;
 double timetick,n_timetick,tmp,sum,u,l,b,pri;

 //Controla los argumentos al programa
 if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((NUM_THREADS = atoi(argv[2])) <= 0) )
  {
    printf("\nUsar: %s n t\n  n: Dimension de la matriz (nxn X nxn) t: Cantidad de threads \n", argv[0]);
    exit(1);
  }
  omp_set_num_threads(NUM_THREADS);

 //Aloca memoria para las matrices
  double triangular_matrix_size = N*(N+1)/2;
  double square_matrix_size = N*N;

 //Aloca memoria para las matrices
  A=(double*)malloc(sizeof(double)*square_matrix_size);
  C=(double*)malloc(sizeof(double)*square_matrix_size);
  B=(double*)malloc(sizeof(double)*square_matrix_size);
  D=(double*)malloc(sizeof(double)*square_matrix_size);
  E=(double*)malloc(sizeof(double)*square_matrix_size);
  F=(double*)malloc(sizeof(double)*square_matrix_size);
  L=(double*)malloc(sizeof(double)*triangular_matrix_size);
  U=(double*)malloc(sizeof(double)*triangular_matrix_size);
  At = (double*)malloc(sizeof(double)*square_matrix_size);
  AUX1 = (double*)malloc(sizeof(double)*square_matrix_size);
  AUX2 = (double*)malloc(sizeof(double)*square_matrix_size);
  AUX3 = (double*)malloc(sizeof(double)*square_matrix_size);
  AUX4 = (double*)malloc(sizeof(double)*square_matrix_size);

  for(i=0;i<N;i++){
   for(j=0;j<N;j++){
    A[i*N+j]=1;
    B[i*N+j]=1;
    C[i*N+j]=1;
    D[i*N+j]=1;
    E[i*N+j]=1;
    F[i*N+j]=1;
   }
  }
  for(i=0;i<triangular_matrix_size;i++){
      U[i]=1;
      L[i]=1;
  }  

  timetick = dwalltime();

  u = promedio(U,triangular_matrix_size);
  l = promedio(L,triangular_matrix_size);
  b = promedio(B,square_matrix_size);

  transpuesta_escalar(At,A,u*l);

  //AxAt
  mult_x_col(AUX2,A,At);

  //AAxC
  mult_x_col(AUX1,AUX2,C);

  //LxB
  mult_low(AUX3,L,B,b);

  //LBxE
  mult_x_col(AUX2,AUX3,E);

  //UxF
  mult_sup(AUX4,U,F,b);

  //DxUF
  mult_x_col(AUX3,D,AUX4);

  // AUX 1 = ulAAC
  // AUX 2 = bLBE
  // AUX 3 = bDUF
  suma(AUX4,AUX1,AUX2,AUX3);

  n_timetick = dwalltime() - timetick;

  printf("Tiempo en segundos %f\n", n_timetick);

  for(i=0;i<N;i++){
    pri = AUX4[i*N];
   for(j=0;j<N;j++){
    check=check&&(AUX4[i*N+j]==pri);
   }
  }

  if(check){
   printf("Multiplicacion de matrices resultado correcto\n");
  }else{
   printf("Multiplicacion de matrices resultado erroneo\n");
  }

 free(AUX4);
 free(AUX3);
 free(AUX2);
 free(AUX1);
 free(A);
 free(B);
 free(C);
 free(D);
 free(E);
 free(F);
 free(L);
 free(U);
 free(At);

 return(0);
}

#include<stdio.h>
#include<stdlib.h>
#include <pthread.h>
#include <sys/time.h>

//Dimension por defecto de las matrices
int N=100;
int TRI_SIZE = 100;
int QUAD_SIZE = 100;
int NUM_THREADS = 1;
int FILAS_POR_THREAD;

double *A,*B,*C,*D,*E,*F,*L,*U,*At,*AUX1,*AUX2,*AUX2,*AUX3,*AUX4;
double u,l,b;
double *u_array,*l_array,*b_array;
pthread_barrier_t barrera;

void print_matrix(double *M);

//Para calcular tiempo
double dwalltime(){
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

void mult_sup(double *dest,double *sup,double *param1,int inicio, int fin, double escalar){
  int i,j,k;
  double sum;
  for(i=inicio;i<fin;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=i;k++){
        sum += sup[k*N+j-k*(k+1)/2] * param1[j*N+k];
      }
      dest[j*N+i] =  sum * escalar;
    }
  }
}

void mult_low(double *dest,double *low,double *param1,int inicio, int fin, double escalar){
  int i,j,k;
  double sum;
  for(i=inicio;i<fin;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<=i;k++){
        sum +=  low[i+j*(j+1)/2] * param1[j*N+k];
      }
      dest[i*N+j] =  sum * escalar;
    }
  }
}


void* mult_x_col(double *dest, double *param1, double *param2, int inicio, int fin){
  int i,j,k;
  double sum;
  for(i=inicio;i<fin;i++){
     for(j=0;j<N;j++){
      sum=0;
      for(k=0;k<N;k++){
        sum +=  param1[i*N+k] * param2[j*N+k];
      }
      dest[i*N+j] =  sum;
    }
  }
}

void transpuesta_escalar(double *dest, double *param,int inicio, int fin,double escalar){
  int i,j;
  for(i=inicio;i<fin;i++){
     for(j=0;j<N;j++){
      dest[i*N+j] = param[j*N+i] * escalar;
    }
  }
}

void* suma(double *dest, double *param1,double *param2,double *param3,int inicio, int fin){
  int i,j,index;
  for(i=inicio;i<fin;i++){
     for(j=0;j<N;j++){
      index = i*N+j; 
      dest[index] =  param1[index] + param2[index] + param3[index];
    }
  }
}

double promedio(double *param,int thread, int size){
  int i;
  int inicio = (thread)*(size/NUM_THREADS); 
  int fin;
  if((thread+1)==NUM_THREADS){
    fin = size;
  }else{
    fin = (thread+1)*(size/NUM_THREADS); 
  }
  double suma=0;
  for(i=inicio;i<fin;i++){
    suma += param[i];
  }
  return suma;
}

double merge_promedio(double *array){
  double dest = 0;
  for ( int i = 0; i < NUM_THREADS; ++i ) {
    dest+=array[i];
  }
  return dest/QUAD_SIZE;
}

void ejecutar_root(void (*func)()){
  int barrier_result = pthread_barrier_wait(&barrera);
  if(barrier_result == PTHREAD_BARRIER_SERIAL_THREAD){
    func();
  } 
  pthread_barrier_wait(&barrera);
}

void post_u(){
  u = merge_promedio(u_array);
}

void post_l(){
  l = merge_promedio(l_array);
}

void post_b(){
  b = merge_promedio(b_array);
}

void* worker(void *arg){
  int thread_id = *(int *)(arg);
  int inicio = thread_id*FILAS_POR_THREAD;
  int fin = (thread_id+1)*FILAS_POR_THREAD;

  u_array[thread_id] = promedio(U,thread_id,TRI_SIZE);
  ejecutar_root(post_u);

  l_array[thread_id] = promedio(L,thread_id,TRI_SIZE);
  ejecutar_root(post_l);
  
  b_array[thread_id] = promedio(B,thread_id,QUAD_SIZE);
  ejecutar_root(post_b);
  
  transpuesta_escalar(At,A,inicio,fin,u*l);
  pthread_barrier_wait(&barrera);
  
  //AxAt
  mult_x_col(AUX2,A,At,inicio,fin);
  pthread_barrier_wait(&barrera);

  //AAxC
  mult_x_col(AUX1,AUX2,C,inicio,fin);

  //LxB
  mult_low(AUX3,L,B,inicio,fin,b);
  pthread_barrier_wait(&barrera);

  //LBxE
  mult_x_col(AUX2,AUX3,E,inicio,fin);

  //UxF
  mult_sup(AUX4,U,F,inicio,fin,b);  
  pthread_barrier_wait(&barrera);

  //DxUF
  mult_x_col(AUX3,D,AUX4,inicio,fin);
  pthread_barrier_wait(&barrera);

  suma(AUX4,AUX1,AUX2,AUX3,inicio,fin);

  pthread_exit(NULL);

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
 int i,j,k;
 int check = 1;
 double timetick,n_timetick,tmp,sum,pri;
 pthread_t *threads;

 //Controla los argumentos al programa
 if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((NUM_THREADS = atoi(argv[2])) <= 0) )
  {
    printf("\nUsar: %s n t\n  n: Dimension de la matriz (nxn X nxn) t: Cantidad de threads \n", argv[0]);
    exit(1);
  }

  threads = (pthread_t *) malloc( NUM_THREADS * sizeof(pthread_t) );
  u_array = (double*)malloc(sizeof(double)*NUM_THREADS);
  l_array = (double*)malloc(sizeof(double)*NUM_THREADS);
  b_array = (double*)malloc(sizeof(double)*NUM_THREADS);
  
  FILAS_POR_THREAD = N / NUM_THREADS;
  TRI_SIZE = N*(N+1)/2;
  QUAD_SIZE = N*N;

 //Aloca memoria para las matrices
  A=(double*)malloc(sizeof(double)*QUAD_SIZE);
  C=(double*)malloc(sizeof(double)*QUAD_SIZE);
  B=(double*)malloc(sizeof(double)*QUAD_SIZE);
  D=(double*)malloc(sizeof(double)*QUAD_SIZE);
  E=(double*)malloc(sizeof(double)*QUAD_SIZE);
  F=(double*)malloc(sizeof(double)*QUAD_SIZE);
  L=(double*)malloc(sizeof(double)*TRI_SIZE);
  U=(double*)malloc(sizeof(double)*TRI_SIZE);
  At = (double*)malloc(sizeof(double)*QUAD_SIZE);
  AUX1 = (double*)malloc(sizeof(double)*QUAD_SIZE);
  AUX2 = (double*)malloc(sizeof(double)*QUAD_SIZE);
  AUX3 = (double*)malloc(sizeof(double)*QUAD_SIZE);
  AUX4 = (double*)malloc(sizeof(double)*QUAD_SIZE);

  for(i=0;i<QUAD_SIZE;i++){
    A[i]=1;
    B[i]=1;
    C[i]=1;
    D[i]=1;
    E[i]=1;
    F[i]=1;
  }

  for(i=0;i<TRI_SIZE;i++){
      U[i]=1;
      L[i]=1;
  }

  pthread_barrier_init(&barrera, NULL, NUM_THREADS);
  
  for (i = 0; i <NUM_THREADS;++i) {
    int *tid;
    tid = (int *) malloc( sizeof(int) );
    *tid = i;
    pthread_create( &threads[i], NULL, worker, (void *)tid );
  }

  timetick = dwalltime();


  for ( i = 0; i < NUM_THREADS; ++i ) {
    pthread_join( threads[i], NULL );
  }

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

#include<stdio.h>
#include<stdlib.h>
#include <pthread.h>

//Dimension por defecto de las matrices
int N=100;
int NUM_THREADS = 1;
double *A,*AT,*C; 
int CANT_FILAS;
pthread_barrier_t barrera;

//Para calcular tiempo
double dwalltime(){
  double sec;
  struct timeval tv;

  gettimeofday(&tv,NULL);
  sec = tv.tv_sec + tv.tv_usec/1000000.0;
  return sec;
}


void * worker( void *arg )
{
  int i, j, k, thread_id, fila_inicio, fila_fin;
  double sum;
  
  thread_id = *(int *)(arg);
  fila_inicio = thread_id * CANT_FILAS;
  fila_fin = (thread_id+1) * CANT_FILAS;

  for(i=fila_inicio;i<fila_fin;i++){
     for(j=0;j<N;j++){
      AT[i*N+j] = A[j*N+i];
    }
  }

  pthread_barrier_wait(&barrera);

  for (i = fila_inicio; i < fila_fin; ++i) {
    for (j = 0; j < N; ++j) { 
      sum = 0; 
      for (k = 0; k < N; ++k) { 
        sum += A[i*N+k] * AT[j*N+k];
      }
      C[i*N+j] = sum;
    }
  }
  pthread_exit( NULL );
}


int main(int argc,char*argv[]){
 int i,j,k;
 int check=1;
 double timetick,n_timetick;
 pthread_t * threads;

 //Controla los argumentos al programa
 if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((NUM_THREADS = atoi(argv[2])) <= 0) )
  {
    printf("\nUsar: %s n t\n  n: Dimension de la matriz (nxn X nxn) t: Cantidad de threads \n", argv[0]);
    exit(1);
  }

  threads = (pthread_t *) malloc( NUM_THREADS * sizeof(pthread_t) );

  CANT_FILAS = N / NUM_THREADS;

  pthread_barrier_init(&barrera, NULL, NUM_THREADS);



 //Aloca memoria para las matrices
  A=(double*)malloc(sizeof(double)*N*N);
  AT=(double*)malloc(sizeof(double)*N*N);
  C=(double*)malloc(sizeof(double)*N*N);

 //Inicializa las matrices A y B en 1, el resultado sera una matriz con todos sus valores en N
  for(i=0;i<N;i++){
   for(j=0;j<N;j++){
    A[i*N+j]=1;
   }
  } 

  for ( i = 0; i < NUM_THREADS; ++i ) {
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

    //Verifica el resultado
       
    for(i=0;i<N;i++){
     for(j=0;j<N;j++){
      check=check&&(C[i*N+j]==N);
     }
    }   

  printf("Tiempo en segundos %f\n", n_timetick);




  if(check){
   printf("Multiplicacion de matrices resultado correcto\n");
  }else{
   printf("Multiplicacion de matrices resultado erroneo\n");
  }



 free(A);
 free(C);
 return(0);
}

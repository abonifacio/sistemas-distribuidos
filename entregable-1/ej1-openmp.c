#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

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

int main(int argc,char*argv[]){
 double *A,*AT,*C;
 int i,j,k;
 int check=1;
 double timetick,n_timetick,temp,sum;

 //Controla los argumentos al programa
 if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((NUM_THREADS = atoi(argv[2])) <= 0) )
  {
    printf("\nUsar: %s n t\n  n: Dimension de la matriz (nxn X nxn) t: Cantidad de threads \n", argv[0]);
    exit(1);
  }

 //Aloca memoria para las matrices
  A=(double*)malloc(sizeof(double)*N*N);
  AT=(double*)malloc(sizeof(double)*N*N);
  C=(double*)malloc(sizeof(double)*N*N);
  omp_set_num_threads(NUM_THREADS);

 //Inicializa las matrices A y B en 1, el resultado sera una matriz con todos sus valores en N
  for(i=0;i<N;i++){
   for(j=0;j<N;j++){
    A[i*N+j]=1;
   }
  }   

  timetick = dwalltime();

  #pragma omp parallel for schedule(static) private(i,j)
  for(i=0;i<N;i++){
   for(j=0;j<N;j++){
    AT[i*N+j]= A[j*N+i];
   }
  }


  #pragma omp parallel for schedule(static) private(i,j,k,sum)
  for(i=0;i<N;i++){
   for(j=0;j<N;j++){
    sum=0;
    for(k=0;k<N;k++){
      sum += A[i*N+k] * AT[j*N+k];
    }
    C[i*N+j]=sum;
   }
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

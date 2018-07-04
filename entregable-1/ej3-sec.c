#include<stdio.h>
#include<stdlib.h>

//Dimension por defecto de las matrices
int N=100;

//Para calcular tiempo
double dwalltime(){
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

int main(int argc,char*argv[]){
 int *A;
 int i;
 double timetick,n_timetick,pares = 0;

 //Controla los argumentos al programa
 if ((argc != 2) || ((N = atoi(argv[1])) <= 0) )
  {
    printf("\nUsar: %s n\n  n: Dimension del vector (nxn X nxn)\n", argv[0]);
    exit(1);
  }

 //Aloca memoria para las matrices
  A=(int*)malloc(sizeof(int)*N);

 //Inicializa las matrices A y B en 1, el resultado sera una matriz con todos sus valores en N
  for(i=0;i<N;i++){
    A[i]=(i+1);
  }   

  timetick = dwalltime();
  for(i=0;i<N;i++){
    if(!(A[i]&1)){
      pares++;
    }
  }
  n_timetick = dwalltime() - timetick;


  printf("Tiempo en segundos %f\n", n_timetick);


  if(pares == (N/2)){
   printf("Conteo de numeros pares resultado correcto (%f)\n",pares);
  }else{
   printf("Conteo de numeros pares resultado erroneo (%f)\n",pares);
  }



 free(A);
 return(0);
}

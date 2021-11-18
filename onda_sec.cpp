#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>

#include <mpi.h>

#define MAXPOINTS 1000
#define MAXSTEPS 1000
#define MINPOINTS 20
#define PI 3.14159265

void init_param(void);
void init_line(void);
void update (void);
void printfinal (void);

int nsteps,                 	/* numero de intervalos de tiempo */
    tpoints; 	     		/* cantidad de puntos en la onda */
double values[MAXPOINTS+2], 	/* valores del campo en el tiempo t */
       oldval[MAXPOINTS+2], 	/* valores en el tiempo (t-dt) */
       newval[MAXPOINTS+2]; 	/* valores en el tiempo (t+dt) */

double t0 = 0.0, tf = 0.0;
int rank, size;
double val_l, val_r;

/***************************************************************************
 *	Input 
 ***************************************************************************/
void init_param(void)
   {
   char tchar[8];

   /* cantidad de puntos e iteraciones */
	tpoints=800;
	nsteps=1000;

   if (rank==0) printf("Usando %d puntos, y  %d pasos\n", tpoints, nsteps);

   }

/***************************************************************************
 *     Inicializacion de puntos en la onda
 **************************************************************************/
void init_line(void)
   {
   int i, j;
   double x, fac, k, tmp;

   /* Valores iniciales con la funcion seno */
   fac = 2.0 * PI;
   k = rank * (tpoints / size) ; 
   tmp = tpoints - 1;
   for (j = rank * (tpoints / size); j < ((tpoints / size) * (rank + 1)); j++) 
   {
      x = k/tmp;
      values[j+1] = sin (fac * x);
      k = k + 1.0;
   } 

   /* Valores iniciales guardados en oldval*/
   for (i = rank * (tpoints / size); i < ((tpoints / size) * (rank + 1)); i++) 
      oldval[i+1] = values[i+1];
   }

/***************************************************************************
 *      Valores nuevos con la ecuacion hiperbolica de onda
 **************************************************************************/
void do_math(int i) //>>REVISAR
   {
   double dtime, c, dx, tau, sqtau;

   dtime = 0.3;
   c = 1.0;
   dx = 1.0;
   tau = (c * dtime / dx);
   sqtau = tau * tau;
   
   newval[i] = (2.0 * values[i]) - oldval[i] 
               + (sqtau * (values[i-1] - (2.0 * values[i]) + values[i+1]));
   }

/***************************************************************************
 *     iteracion en una cantidad determinada de pasos
 **************************************************************************/
void update()
{
int i, j;

/* Actualizacion de valores en cada paso */
for (i = 1; i<= nsteps; i++) {
   /* Actualizacion de valores para cada punto de la onda*/
   for (j = rank * (tpoints / size); j < ((tpoints / size) * (rank + 1)); j++) {
      /* puntos internos y de frontera */
      if ((j == 0) || (j  == tpoints - 1))
         newval[j+1] = 0.0;
      else
         do_math(j+1);
   }

   /* Actualizacion de valores antiguos y nuevos */
   for (j = rank * (tpoints / size); j < ((tpoints / size) * (rank + 1)); j++) {
      oldval[j+1] = values[j+1];
      values[j+1] = newval[j+1];
      }
   }
}

/***************************************************************************
 *     Output
 **************************************************************************/
void printfinal()
{
   int i,x;

   for (i = 1; i <= tpoints; i++) {
      printf("%i: %6.4f \n", i, values[i]);
   //  if (i%10 == 0)
      //   printf("\n");
   }
}
/***************************************************************************
 *	Programa principal
 **************************************************************************/
int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

   t0 = MPI_Wtime();  // tiempo inicial

   if (rank==0) printf("Onda secuencial ...\n");

   init_param();
   
   if (rank==0) printf("Inicializando puntos de la onda...\n");
   
   init_line();

   // SYN extremo izq
   if (rank == 0)
   {
      MPI_Send(&values[(rank+1) * (tpoints / size)], 1, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
   }
   else if (rank == size-1)
   {
      MPI_Recv(&values[rank * (tpoints / size)], 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }
   else
   {
      MPI_Send(&values[(rank+1) * (tpoints / size)], 1, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
      MPI_Recv(&values[rank * (tpoints / size)], 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   // SYN extremo Der
   if (rank == size-1)
   {
      MPI_Send(&values[(rank * (tpoints / size))+1], 1, MPI_DOUBLE, rank-1, rank, MPI_COMM_WORLD);
      
   }
   else if (rank == 0)
   {
      MPI_Recv(&values[((rank+1) * (tpoints / size))+1], 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }
   else
   {
      MPI_Send(&values[(rank * (tpoints / size))+1], 1, MPI_DOUBLE, rank-1, rank, MPI_COMM_WORLD);
      MPI_Recv(&values[((rank+1)* (tpoints / size))+1], 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }
   MPI_Barrier(MPI_COMM_WORLD);


   if (rank==0) printf("Actualizando puntos para la cantidad total de pasos...\n");
   
   update();

   // SYN COMM
   for (int j = rank * (tpoints / size); j < ((tpoints / size) * (rank + 1)); j++)
   {
      if (rank != 0) MPI_Send(&values[j+1], 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
      
      if (rank == 0)
      {
         for (int z = 1; z < size; z++)
         {
            MPI_Recv(&values[((z * (tpoints / size))+1)+j], 1, MPI_DOUBLE, z, z, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         }
      }
   }
   
   MPI_Barrier(MPI_COMM_WORLD);

   if (rank==0) 
   {
      printf("Imprimiendo resultados...\n");
      printfinal();
      printf("\nFin.\n\n");
      tf = MPI_Wtime();  // tiempo final
      printf("Tiempo de ejecucion: %1.6f ms, con %d cores\n", 1000*(tf-t0), size);
   }

   MPI_Finalize();
}

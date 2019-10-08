#include <stdio.h>
#include <math.h>
#include "mpi.h"

/* Author: Shir Chen                                 */
/* This program solves the heat equation for the     */
/* following geometric figure, on 4 processors only. */
/* inner boundaries: 20, outer boundaries: 60        */
/*                                                   */
/*                       200                         */
/*               --------------------                */
/*              |                    |               */
/*           100|                    |100            */
/*              |        100         |               */
/*  ------------       --------       -------------  */
/* |                  |        |                   | */
/* | 100          100 |        | 100           100 | */
/* |                  |        |                   | */
/*  ------------       --------       -------------  */
/*     150      |                    |     150       */
/*           100|                    |100            */
/*              |                    |               */
/*               --------------------                */
/*                       200                         */

#define X	200 // x dimension of local matrix
#define Y	100 // y dimension of local matrix
#define TOTALX	500 // x dimension of total problem
#define TOTALY	300 // y dimension of total problem
#define OUTERBoundary	60 // boundary value
#define INNERBoundary	20 // boundary value
#define STEPS	100 // number of steps

void initialize( double M[X+2][Y+2],double totalM[TOTALX+2][TOTALY+2],int rank  );
void local2matrix( double M[X+2][Y+2],double totalM[TOTALX+2][TOTALY+2],int rank );
void matrix2file( double M[TOTALX+2][TOTALY+2] );

int main(int argc,char *argv[] )
{
    int        rank, value, size, i, j, itcnt;
    MPI_Status status;
    double     diffnorm, gdiffnorm;
    double     localM[X+2][Y+2];
	double     transferM[50];
    double     newM[X+2][Y+2];
	double     totalM[TOTALX+2][TOTALY+2];
	double startwtime = 0.0, endwtime;
	
	// get my rank(id) and how many possesses are running
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
	// this problem is defined for only 4 processes
    if (size != 4) MPI_Abort( MPI_COMM_WORLD, 1 );
	// first, initialize the local matrix and the total problem matrix
	startwtime = MPI_Wtime();
	initialize(localM,totalM,rank);
    itcnt = 0; // for iteration count
	do {
		// for sending, we take the relevant data from the local matrix and transfer it to a temp matrix-transferMatrix
		// for receiving, we need to transfer back the data from transferMatrix to the process local matrix
		switch(rank)
		{
			case 0:
				for(i=0;i<50;i++)
					transferM[i]=localM[i+1][1];
				MPI_Send( transferM, 50, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD );
				for(i=0;i<50;i++)
					transferM[i]=localM[i+150+1][1];			
				MPI_Send( transferM, 50, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD );
				MPI_Recv( transferM, 50, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, &status );
				for(i=0;i<50;i++)
					localM[i+1][0]=transferM[i];		
				MPI_Recv( transferM, 50, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, &status );
				for(i=0;i<50;i++)
					localM[i+150+1][0]=transferM[i];				
				break;
			case 1:
				for(i=0;i<50;i++)
					transferM[i]=localM[i+1][100+1];
				MPI_Send( transferM, 50, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD );
				for(i=0;i<50;i++)
					transferM[i]=localM[i+150+1][100+1];			
				MPI_Send( transferM, 50, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD );
				MPI_Recv( transferM, 50, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, &status );
				for(i=0;i<50;i++)
					localM[i+1][100]=transferM[i];		
				MPI_Recv( transferM, 50, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, &status );
				for(i=0;i<50;i++)
					localM[i+150+1][100]=transferM[i];	
				break;
			case 2:
				for(i=0;i<50;i++)
					transferM[i]=localM[i+150+1][100+1];
				MPI_Send( transferM, 50, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
				for(i=0;i<50;i++)
					transferM[i]=localM[i+150+1][1];			
				MPI_Send( transferM, 50, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD );
				MPI_Recv( transferM, 50, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status );
				for(i=0;i<50;i++)
					localM[i+150+1][100]=transferM[i];		
				MPI_Recv( transferM, 50, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status );
				for(i=0;i<50;i++)
					localM[i+150+1][0]=transferM[i];	
				break;
			case 3:
				for(i=0;i<50;i++)
					transferM[i]=localM[i+1][100+1];
				MPI_Send( transferM, 50, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
				for(i=0;i<50;i++)
					transferM[i]=localM[i+1][1];			
				MPI_Send( transferM, 50, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD );
				MPI_Recv( transferM, 50, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status );
				for(i=0;i<50;i++)
					localM[i+1][100]=transferM[i];		
				MPI_Recv( transferM, 50, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status );
				for(i=0;i<50;i++)
					localM[i+1][0]=transferM[i];	
				break;
		}
	/* Compute new values (but not on boundary) */
	itcnt ++;
	diffnorm = 0.0;
	for (i=1; i<=200; i++) 
	    for (j=1; j<=100; j++) {
		newM[i][j] = (localM[i][j+1] + localM[i][j-1] +
			      localM[i+1][j] + localM[i-1][j]) / 4.0;
		diffnorm += (newM[i][j] - localM[i][j]) * 
		            (newM[i][j] - localM[i][j]);
	    }
	/* Only transfer the interior points */
	for (i=1; i<=200; i++) 
	    for (j=1; j<=100; j++) 
		localM[i][j] = newM[i][j];

	MPI_Allreduce( &diffnorm, &gdiffnorm, 1, MPI_DOUBLE, MPI_SUM,
		       MPI_COMM_WORLD );
	gdiffnorm = sqrt( gdiffnorm );
	if (rank == 0) printf( "At iteration %d, diff is %e\n", itcnt, 
			       gdiffnorm );
    } while (gdiffnorm > 1.0e-2 && itcnt < STEPS);
	// master process is incharge of collecting all data from workers and updating the 
	// totalMatrix with the correct data.
	if(rank==0){
		// update totalM with master data
		local2matrix(localM,totalM,rank);
		// for each worker update totalM with its data after receiving its matrix
		for(i=1;i<4;i++){
			MPI_Recv(localM, (X+2)*(Y+2), MPI_DOUBLE, i,2, MPI_COMM_WORLD, &status);
			local2matrix(localM,totalM,i);
		}
		// write total problem matrix to file
		matrix2file(totalM);
		printf("output updated\n");
		endwtime = MPI_Wtime();
        printf("execution time = %f\n", endwtime-startwtime);
	}else{
		// for workers only send their data to master process
		MPI_Send(localM, (X+2)*(Y+2), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
	}
		
	
    MPI_Finalize( );
    return 0;
}
/* Initialize the matrix and fill the data as specified, */
/* inner and outer boundaries.                           */
/* */
void initialize( double M[X+2][Y+2],double totalM[TOTALX+2][TOTALY+2],int rank ){
	int i,j;
	// total problem matrix
	//firstly, all matrix to outer boundaries
	for (i=0; i<TOTALX+2; i++) 
		for (j=0; j<TOTALY+2; j++)
			totalM[i][j] = OUTERBoundary;
	//matrix to inner boundaries-middle square
	for (i=200+1; i<=300; i++) 
		for (j=100+1; j<=200; j++)
			totalM[i][j] = INNERBoundary;
    // local matrix	
	//firstly, all matrix to outer boundaries
    for (i=0; i<X+2; i++) 
		for (j=0; j<Y+2; j++)
			M[i][j] = OUTERBoundary;				
	//matrix to inner boundaries according to rank
		switch(rank)
	{
		case 0://upper rectangle
			for (i=50+1; i<=150; i++) 
				M[i][0] = INNERBoundary;
			break;
		case 1://lower rectangle
			for (i=50+1; i<=150; i++) 
				M[i][100+1] = INNERBoundary;
			break;
		case 2://left rectangle
			for (j=1; j<=100; j++) 
				M[200+1][j] = INNERBoundary;
			break;
		case 3://right rectangle
			for (j=1; j<=100; j++) 
				M[0][j] = INNERBoundary;
			break;
	}
					

}
/* write M matrix to a text file as output*/
void matrix2file( double M[TOTALX+2][TOTALY+2] ){
	FILE *outfile; 
    int i,j;  
    // open file for writing 
    outfile = fopen ("outputMatrix.txt", "w"); 
    // write to file  
    for (i=0; i<TOTALX+2; i++)
	{
		for (j=0; j<TOTALY+2; j++) {
			fprintf(outfile, "%f,",M[i][j]);
			printf("%f,",M[i][j]);
		}
			
		fprintf(outfile,"\n");
		printf("\n");
	}		
	// close file 
    fclose (outfile);
}
/* according to rank updating the totla problem matrix-totalM with local matrix-M data*/
void local2matrix( double M[][Y+2],double totalM[TOTALX+2][TOTALY+2],int rank){
	int i,j;
		switch(rank)
	{
		case 0://upper rectangle
			for (i=0; i<X+2; i++)
				for (j=0; j<Y+2; j++) 
					totalM[150+i][200+j]=M[i][j];
			break;
		case 1://lower rectangle
			for (i=0; i<X+2; i++)
				for (j=0; j<Y+2; j++) 
					totalM[150+i][j]=M[i][j];
			break;
		case 2://left rectangle
			for (i=0; i<X+1; i++)
				for (j=0; j<Y+1; j++) 
					totalM[i][100+j]=M[i][j];
			break;
		case 3://right rectangle
			for (i=1; i<X+2; i++)
				for (j=1; j<Y+2; j++) 
					totalM[300+i][100+j]=M[i][j];
			break;
			
	}	

}
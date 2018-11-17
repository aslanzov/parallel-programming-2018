#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "mpi.h" 
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) 
{
	int n=3000;
	int ProcNum, ProcRank;
	double paralltime1, paralltime2, linetime1, linetime2;
	int count, ost;
	int *res;
	int *reslin;
	int *vector;
	int *matrix_transport; 
	int **matrix;  
	int *tmp;
	int *tmplong;
	int *SendCountMatrix; //количество столбцов в каждый процесс
	int *SendIndexMatrix; //смещение элементов матрицы при транспонировани
	int *SendCountVector; //количество элементов вектора в каждый процесс
	int *SendIndexVector; //смещение элементов вектора
	int *rbufMatrix;
	int *rbufVector;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	SendCountMatrix = new int[ProcNum];
	SendIndexMatrix = new int[ProcNum];
	SendCountVector = new int[ProcNum];
	SendIndexVector = new int[ProcNum];
	vector = new int[n];
	matrix_transport = new int[n*n];
	res = new int[n];
	reslin = new int[n];
	matrix = new int*[n];
	for (int i = 0; i<n; i++)
		matrix[i] = new int[n];


	if (ProcRank == 0)
	{
		for (int i = 0; i<n; i++)
			for (int j = 0; j < n; j++)
			{
				matrix[i][j] = rand() % 10;
			}

		int k = 0;
		for (int j = 0; j<n; j++)
			for (int i = 0; i<n; i++, k++)
				matrix_transport[k] = matrix[i][j];

		for (int i = 0; i<n; i++)
		{
			vector[i] = rand() % 10;
			res[i] = 0;
		}
	}
	
	paralltime1 = MPI_Wtime();
	
	count = n / ProcNum;
	ost = n%ProcNum;

	for (int i = 0; i<ost; i++)
	{
		SendCountMatrix[i] = n*(count + 1);
		SendIndexMatrix[i] = i*SendCountMatrix[i];
		SendCountVector[i] = count + 1;
		SendIndexVector[i] = i*SendCountVector[i];
	}
	for (int i = ost; i<ProcNum; i++)
	{
		SendCountMatrix[i] = n*count;
		SendIndexMatrix[i] = n*(ost + i*count);
		SendCountVector[i] = count;
		SendIndexVector[i] = ost + i*count;
	}


	rbufMatrix = new int[SendCountMatrix[ProcRank]];
	rbufVector = new int[SendCountVector[ProcRank]];

	
	MPI_Scatterv(matrix_transport, SendCountMatrix, SendIndexMatrix, MPI_INT, rbufMatrix, SendCountMatrix[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(vector, SendCountVector, SendIndexVector, MPI_INT, rbufVector, SendCountVector[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);

	tmplong = new int[SendCountMatrix[ProcRank]];
	tmp = new int[n];
	
	for (int i = 0; i<n; i++)
		tmp[i] = 0;
	for (int i = 0; i<SendCountMatrix[ProcRank]; i++)
		tmplong[i] = 0;

	int k = 0;
	for (int i = 0; i<SendCountVector[ProcRank]; i++)
		for (int j = 0; j<n; j++, k++)
			tmplong[k] = rbufMatrix[j + i*n] * rbufVector[i];
	if (SendCountMatrix[ProcRank] == n)
		tmp = tmplong;
	if (SendCountMatrix[ProcRank]>n)
	{
		for (int i = 0; i<n; i++)
			for (int j = 0; j<SendCountVector[ProcRank]; j++)
			{
				tmp[i] += tmplong[i + j*n];
			}
	}

	MPI_Reduce(tmp, res, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_Comm MPI_COMM_WORLD); 
	paralltime2 = MPI_Wtime();

	if (ProcRank == 0)
	{
		for (int i = 0; i<n; i++)
			reslin[i] = 0;

		linetime1 = MPI_Wtime();

		for (int j = 0; j<n; j++)
			for (int i = 0; i<n; i++)
				reslin[i] += matrix[i][j] * vector[j];

		linetime2 = MPI_Wtime();

		int prov = 0;
		for (int i = 0; i<n; i++)
			if (res[i] != reslin[i])
				prov++;
		if (prov == 0)
			cout << "Correct" << endl;
		else
			cout << "Error" << endl;

		cout << "SingleTime: " << linetime2 - linetime1 << endl;
		cout << "MultiTime: " << paralltime2 - paralltime1 << endl;
	}

	MPI_Finalize();
		
	return 0;
}
#include <iostream>
#include "mpi.h"
#include <time.h>

int* createArray(int &length)
{
	int *arr;
	std::cout << "Input an array size:" << std::endl;
	std::cin >> length;
	arr = new int[length];
	srand(time(0));
	for (int i = 0; i < length; i++)
		arr[i] = rand();
	return arr;
}

int findMin(int *arr, int length)
{
	int min;
	min = arr[0];
	for (int i = 1; i < length; i++)
	{
		if (arr[i] < min)
			min = arr[i];
	}
	return min;
}

int main(int argc, char* argv[])
{
	int n;
	double startTime;
	int errCode;
	MPI_Status Status;

	if ((errCode = MPI_Init(&argc, &argv)) != 0)
	{
		return errCode;
	}

	int ProcNum, ProcRank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
	{
		int *array, length;
		
		array = createArray(length);
		
		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Send(&length, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		n = length / ProcNum;
		startTime = MPI_Wtime();
		for (int i = 1; i < ProcNum - 1; i++)
		{
			MPI_Send(array + n *(i - 1), n, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		MPI_Send(array + n * (ProcNum - 1 - 1), n, MPI_INT, ProcNum - 1, 0, MPI_COMM_WORLD);


		int min, tmp;

		min = array[0];

		for (int i = n * (ProcNum - 1) - 1; i < length; i++)
		{
			if (array[i] < min)
				min = array[i];
		}

		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Recv(&tmp, 1, MPI_INT, MPI_ANY_SOURCE,
				MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			if (tmp < min)
				min = tmp;
		}
	

		std::cout << "MultiProcess result = " << min << std::endl;
		std::cout<< "Time:" << MPI_Wtime()-startTime<<std::endl;
		startTime = MPI_Wtime();
		min = findMin(array, length);

		std::cout << "SingleProcess result = " << min << std::endl;
		std::cout << "Time:" << MPI_Wtime() - startTime << std::endl;
	}
	else
	{
		int* buf, length;

		MPI_Recv(&length, 1, MPI_INT, 0,
			0, MPI_COMM_WORLD, &Status);

		n = length / ProcNum;

		buf = new int[n];

		MPI_Recv(buf, n, MPI_INT, 0,
			0, MPI_COMM_WORLD, &Status);

		int min = buf[0];

		for (int i = 1; i < n; i++)
		{
			if (buf[i] < min)
				min = buf[i];
		}

		MPI_Send(&min, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

	}

	MPI_Finalize();



	return 0;
}

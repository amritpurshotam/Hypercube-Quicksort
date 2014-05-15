//MPI and OpenMP Hybrid

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

int size = 100; //current number of elements in list
int ne = 100; //number of elements required in list
int origsize = 100; //original number of elements in list
int ns = 0; //number of elements sent to partner node
int nr = 0; //number of elements received from partner node
int d = 2; // dimension of hypercube
int *list;

// prints the list
void printlist()
{
	int i = 0;	
	for (i = 0 ; i < size ; i++)
	{
		printf("%d ", list[i]);
	}
	printf("\n");
}

// calculates num^exp
int power(int num, int exp)
{
	int i = 0;
	int ans = num;
	for (i = 1 ; i < exp ; i++)
	{
		ans = ans * num;
	}
	return ans;
}

//rounds off num
int roundoff(double num)
{
	int ans = (int) num;
	num = num - ans;
	if (num >= 0.5) { ans = ans + 1; }
	return ans;
}

// calculates the mean of the list
int calcMean()
{
	int i = 0;
	double mean = 0;
	for (i = 0 ; i < size ; i++)
	{
		mean = mean + list[i];
	}
	mean = roundoff(mean / size);
	int ans = (int) mean;
	return ans;
}

int medianofthree()
{
	int median[3];
	median[0] = list[0];
	median[1] = list[(int) (size/2)];
	median[2] = list[size-1];
	int i = 0;
	for (i = 0 ; i < 2 ; i++)
	{
		int j = i+1;
		for (j = i+1 ; j < 3 ; j++)
		{
			if (median[i] > median[j])
			{
				int temp = median[i];
				median[i] = median[j];
				median[j] = temp; 
			}
		} 
	}
	return median[1];
}

// normal quicksort algorithm
void quicksort(int left, int right)
{	
	if (left < right)
	{
		int i = left;
		int j = right;
		int pivot = list[right];
		int temp = 0;
		
		while (i < j)
		{
			while (list[i] < pivot)
			{
				i++;
			}
			while ((j > i) && (list[j] >= pivot))
			{
				j--;
			}
			if (j > i)
			{
				temp = list[i];
				list[i] = list[j];
				list[j] = temp;
			}
		}
		temp = list[i];
		list[i] = list[right];
		list[right] = temp;
		quicksort(left, i-1);
		quicksort(i+1, right);
	}
}

void main(int argc, char *argv[])
{
	// dynamically allocating memory to list	
	list = (int *) malloc(size * sizeof(int));

	//MPI variable declarations and initialisations
	int numNodes = power(2, d);
	int myid, p;
	int source, dest;
	int tag = 0;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	//if (myid == 0) {printf("List Size: %d\n", origsize*numNodes);}
	struct timeval startTime;
	struct timeval endTime;
	gettimeofday(&startTime, NULL);

	// building the groups and their communicators
	MPI_Group cubeg [d+1][power(2, d-1)];
	MPI_Comm cubec [d+1][power(2, d-1)]; 
	MPI_Comm_group(MPI_COMM_WORLD, &cubeg[d][0]);
	cubec[d][0] = MPI_COMM_WORLD;
	int *rank;

	int r = 0; // row number
	for (r = 1 ; r < d ; r++)
	{
		int len = numNodes / power(2, d-r); // length of rank
		rank = (int *) malloc(len * sizeof(int));
		int n = 0; // node number
		while (n < numNodes)
		{
			int l;
			for (l = 0 ; l < len ; l++)
			{
				rank[l] = n;
				n++;
			}
			MPI_Group_incl(cubeg[d][0], len, rank, &cubeg[r][(n/(power(2,r)))-1]);
			MPI_Comm_create(cubec[d][0], cubeg[r][(n/(power(2,r)))-1], &cubec[r][(n/(power(2,r)))-1]);
		}
	}
	
	
	// creating a list of numbers and distributing it across all the nodes (code finished and working)
	if (myid == 0)
	{
		int i = 0;
		int randnum = 0;
		int globallist[size*numNodes];
		
		#pragma omp parallel for
		for (i = 0; i < size*numNodes ; i++)
		{
			randnum = rand() % (10000000/*size*numNodes*/ + 1);
			globallist[i] = randnum;
			if (i < size) { list[i] = randnum; }
		}
		
		#pragma omp parallel for
		for (i = 1 ; i < numNodes ; i++)
		{
			MPI_Send(&globallist[size*i], size, MPI_INT, i, tag, MPI_COMM_WORLD);
		}
	}
	else { MPI_Recv(&list[0], size, MPI_INT, 0, tag, MPI_COMM_WORLD, &status); }
	
	// starting the hypercube quicksort algorithm
	int bitvalue = power(2, d-1);
	int mask = power(2, d) - 1;
	int pivot = 0;
	int l = 0;
	int partner = 0;
	int i, j;
	int nprocs_cube;
	for (l = d ; l >= 1 ; l--)
	{
		nprocs_cube = power(2, l);
		if ((myid & mask) == 0)
		{
			//pivot = medianofthree();
			pivot = calcMean();
		}
		
		MPI_Bcast(&pivot, 1, MPI_INT, 0, cubec[l][myid/nprocs_cube]);
		
		// partitions list so that list[0:j] <= pivot < list[j+1, size] (working)	
		i = 0;
		j = size-1;
		int temp = 0;
		while (i < j)
		{
			while (list[i] < pivot)
			{
				i++;
			}
			while (list[j] > pivot)
			{
				j--;
			}
			if (j > i)
			{
				if (list[i] == list[j]) 
				{ 
					i++; 
				}
				else
				{	
					temp = list[i];
					list[i] = list[j];
					list[j] = temp;
				}
			}
		}
	
		partner = myid ^ bitvalue; // XOR
		if ((myid & bitvalue) == 0)
		{
			ns = size - (j + 1); // length of right sublist
			MPI_Send(&ns, 1, MPI_INT, partner, tag, MPI_COMM_WORLD); 
			MPI_Recv(&nr, 1, MPI_INT, partner, tag, MPI_COMM_WORLD, &status);
			ne = size - ns + nr;
			if (ne > size)
			{
				list = (int *) realloc (list, ne * sizeof(int));
			}
			size = ne;
			//printf("myid: %d, l: %d, ns: %d, nr: %d, ne: %d, size: %d list before send\n", myid, l, ns, nr, ne, size);
			MPI_Send(&list[j+1], ns, MPI_INT, partner, tag, MPI_COMM_WORLD);
			//printf("myid: %d, l: %d, list sent\n", myid, l);
			MPI_Recv(&list[j+1], nr, MPI_INT, partner, tag, MPI_COMM_WORLD, &status);
		}
		else
		{
			ns = (j + 1); // length of left sublist
			MPI_Send(&ns, 1, MPI_INT, partner, tag, MPI_COMM_WORLD);
			MPI_Recv(&nr, 1, MPI_INT, partner, tag, MPI_COMM_WORLD, &status);
			ne = size - ns + nr;
			if (ne > size)
			{
				list = (int *) realloc (list, ne * sizeof(int));
			}
			//printf("myid: %d, l: %d, ns: %d, nr: %d, ne: %d, size: %d list before send\n", myid, l, ns, nr, ne, size);
			MPI_Send(&list[0], ns, MPI_INT, partner, tag, MPI_COMM_WORLD);
			//printf("myid: %d, l: %d, list sent\n", myid, l);
			//shift right sublist to the left
			for (i = 0 ; i < size - ns ; i++)
			{
				list[i] = list[j+1+i];
			}
			MPI_Recv(&list[size-ns], nr, MPI_INT, partner, tag, MPI_COMM_WORLD, &status);
			size = ne;
		}
		mask = mask ^ bitvalue;
		bitvalue = bitvalue / 2;
	}
	quicksort(0, size-1);
	//printf("myid: %d, l: %d, local list sorted\n", myid, l);

	//assemble back into 1 big list
	if (myid != 0)
	{
		MPI_Send(&size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
		MPI_Send(&list[0], size, MPI_INT, 0, tag, MPI_COMM_WORLD);
	}
	else
	{
		int i; // index of globallist
		int j; // index of numNodes
		int globallist[origsize * numNodes];
		for (i = 0 ; i < size ; i++) { globallist[i] = list[i]; }
		
		for (j = 1 ; j < numNodes ; j++)
		{
			i = size;
			MPI_Recv(&size, 1, MPI_INT, j, tag, MPI_COMM_WORLD, &status);
			int templist[size];
			MPI_Recv(&templist[0], size, MPI_INT, j, tag, MPI_COMM_WORLD, &status);
			size = i + size;
			int k = 0; // index of templist
			while (i < size) 
			{
				globallist[i] = templist[k];
				k++;
				i++;
			}
		}
		
		gettimeofday(&endTime, NULL);
		double tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
		double tE = endTime.tv_sec*1000000 + (endTime.tv_usec);
		printf("Total Time Taken: %f seconds\n", (tE - tS)/1000000);	
		
		//Prints the entire sorted list 
		/*for (i = 0 ; i < origsize*numNodes ; i++)
		{
			printf("%d ", globallist[i]);
		}
		printf("\n");*/
	}
	MPI_Finalize();
}

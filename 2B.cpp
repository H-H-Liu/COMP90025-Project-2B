// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/
// with many fixes and changes for multiple sequence alignment and to include an MPI driver
#include <mpi.h>
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include "sha512.hh"
#include <cmath>
#include <cstdlib>
#include<omp.h>

using namespace std;

std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap, int *penalties, int root);
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);
void do_MPI_task(int rank, int root);

/*
Examples of sha512 which returns a std::string
sw::sha512::calculate("SHA512 of std::string") // hash of a string, or
sw::sha512::file(path) // hash of a file specified by its path, or
sw::sha512::calculate(&data, sizeof(data)) // hash of any block of data
*/

// Return current time, for performance measurement
uint64_t GetTimeStamp()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;

// Driver code
int main(int argc, char **argv)
{
	int rank;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &prov);
	MPI_Comm_rank(comm, &rank);
	if (rank == root)
	{
		int misMatchPenalty;
		int gapPenalty;
		int k;
		std::cin >> misMatchPenalty;
		std::cin >> gapPenalty;
		std::cin >> k;
		std::string genes[k];
		for (int i = 0; i < k; i++)
			std::cin >> genes[i];

		int numPairs = k * (k - 1) / 2;

		int penalties[numPairs];

		uint64_t start = GetTimeStamp();

		// return all the penalties and the hash of all allignments
		std::string alignmentHash = getMinimumPenalties(genes,
														k, misMatchPenalty, gapPenalty,
														penalties, root);

		// print the time taken to do the computation
		printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - start));

		// print the alginment hash
		std::cout << alignmentHash << std::endl;

		for (int i = 0; i < numPairs; i++)
		{
			std::cout << penalties[i] << " ";
		}
		std::cout << std::endl;
	}
	else
	{
		// do stuff for MPI tasks that are not rank==root
		do_MPI_task(rank, root);
	}
	MPI_Finalize();
	return 0;
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/

struct data
{
	int ID;
	int cuml_size;
	std::string S1;
	std::string S2;
	int result;
	std::string alignmentHash;
};

int compare(const void *node1, const void *node2)
{
	string name1 = ((const struct data *)node1)->cuml_size;
	string name2 = ((const struct data *)node2)->cuml_size;

	if (name1 <= name2)
	{
		return -1;
	}
	else if (name1 > name2)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int min3(int a, int b, int c)
{
	if (a <= b && a <= c)
	{
		return a;
	}
	else if (b <= a && b <= c)
	{
		return b;
	}
	else
	{
		return c;
	}
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
int **new2d(int width, int height)
{
	int **dp = new int *[width];
	size_t size = width;
	size *= height;
	int *dp0 = new int[size];
	if (!dp || !dp0)
	{
		std::cerr << "getMinimumPenalty: new failed" << std::endl;
		exit(1);
	}
	dp[0] = dp0;
	for (int i = 1; i < width; i++)
		dp[i] = dp[i - 1] + height;

	return dp;
}

// called by the root MPI task only
// this procedure should distribute work to other MPI tasks
// and put together results, etc.
std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap,
								int *penalties, int root)
{
	double tmp = (k - 1) * k;
	int size = (int)ceil(tmp / 2);
	//	int kek = floor(50);
	struct data *database = (struct data *)malloc(sizeof(struct data) * (size + 2));

	int probNum = 0;
	for (int i = 1; i < k; i++)
	{
		for (int j = 0; j < i; j++)
		{
			std::string gene1 = genes[i];
			std::string gene2 = genes[j];
			int m = gene1.length(); // length of gene1
			int n = gene2.length(); // length of gene2
			int l = m + n;
			struct data sd;
			sd.ID = probNum;
			sd.S1 = gene1;
			sd.S2 = gene2;
			sd.result = -999;
			sd.cuml_size = l;
			database[probNum] = sd;
			probNum++;
		}
	}
	int proc = MPI_Comm_size(MPI_COMM_WORLD, &proc);
	// int highest_rank = size - 1;

	std::qsort(database, size, sizeof(data), compare);

	for (int i = 0; i < probNum; i++)
	{
		std::cout << database[i] << " " << std::endl;
	}

	int pivot = 0;
	int index = 0;

	int *items = new int[size];
	items = {0};
	std::string alignmentHash = "";

	#pragma omp parallel
	{
		#pragma omp master
		{
			while (true)
			{
				if (pivot >= size)
				{
					pivot = 0;
				}
				int red;
				#pragma omp atomic read
				red = items[pivot];

				if (red == 0)
				{
					#pragma omp atomic update
					items[pivot]++;
					#pragma omp single firstprivate(index, pivot)
					{
						MPI_Send(&database[index].S1.length(), 1, MPI_INT, pivot, 7, MPI_COMM_WORLD);
						MPI_Send(&database[index].S2.length(), 1, MPI_INT, pivot, 8, MPI_COMM_WORLD);
						MPI_Send(database[index].S1.c_str(), database[index].S1.length(), MPI_CHAR, pivot, 1, MPI_COMM_WORLD);
						MPI_Send(database[index].S2.c_str(), database[index].S2.length(), MPI_CHAR, pivot, 2, MPI_COMM_WORLD);
						// MPI_Send(database[index].,database[index].S1.length(),MPI_CHAR,pivot,3);
						//	MPI_Send(database[index].S1.length(),1,MPI_INT,pivot,4,MPI_COMM_WORLD);
						//	MPI_Send(database[index].S2.length(),1,MPI_INT,pivot,5,MPI_COMM_WORLD);

						MPI_Send(&pxy, 1, MPI_INT, pivot, 4, MPI_COMM_WORLD);
						MPI_Send(&pgap, 1, MPI_INT, pivot, 5, MPI_COMM_WORLD);

						char *buf1;
						char *buf2;
						int pen, len;
						int tag1;
						int tag2;
						int tag3;
						MPI_Status sta1;
						MPI_Status sta2;

						MPI_Recv(&pen, 1, MPI_INT, pivot, tag1, MPI_COMM_WORLD, sta1);
						MPI_Recv(&len, 1, MPI_INT, pivot, tag2, MPI_COMM_WORLD, sta1);
						MPI_Recv(&buf1, , MPI_CHAR, pivot, tag3, MPI_COMM_WORLD, sta1);
						//	MPI_Recv(&buf2, database[index].S2.length(), MPI_CHAR, pivot, tag2, MPI_COMM_WORLD, sta2);

						std::string a(buf1);

						#pragma omp atomic write
						alignmentHash = sw::sha512::calculate(alignmentHash.append(a));

						#pragma omp atomic update
						items[pivot]--;
					}
					index++;
				}
				pivot++;
			}
		}
	}
}

// called for all tasks with rank!=root
// do stuff for each MPI task based on rank
void do_MPI_task(int rank, int root)
{
	char *buf1, *buf2;
	int l1, l2, pxy, pgap;
	int t1, t2, t3, t4, t5, t6;
	MPI_Status sta1, sta2, sta3, sta4, sta5, sta6;
	MPI_Recv(l1, 1, MPI_CHAR, root, t1, MPI_COMM_WORLD, sta1);
	MPI_Recv(l2, 1, MPI_CHAR, root, t2, MPI_COMM_WORLD, sta2);
	MPI_Recv(buf1, l1, MPI_CHAR, root, t3, MPI_COMM_WORLD, sta3);
	MPI_Recv(buf2, l2, MPI_CHAR, root, t4, MPI_COMM_WORLD, sta4);
	MPI_Recv(pxy, 1, MPI_INT, root, t5, MPI_COMM_WORLD, sta5);
	MPI_Recv(pgap, 1, MPI_INT, root, t6, MPI_COMM_WORLD, sta6);

	std::string a(buf1);
	std::string b(buf2);

	int *xans, *yans;
	int penalty = getMinimumPenalty(a, b, pxy, pgap, xans, yans);

	int l = a.length() + b.length();

	int id = 1;
	int a;
	for (a = l; a >= 1; a--)
	{
		if ((char)yans[a] == '_' && (char)xans[a] == '_')
		{
			id = a + 1;
			break;
		}
	}
	std::string align1 = "";
	std::string align2 = "";
	for (a = id; a <= l; a++)
	{
		align1.append(1, (char)xans[a]);
	}
	for (a = id; a <= l; a++)
	{
		align2.append(1, (char)yans[a]);
	}
	std::string align1hash = sw::sha512::calculate(align1);
	std::string align2hash = sw::sha512::calculate(align2);

	std::string problemhash = sw::sha512::calculate(a.append(b));

	MPI_Send(penalty, 1, MPI_INT, root, 1, MPI_COMM_WORLD);
	MPI_Send(problemhash, problemhash.length(), MPI_CHAR, root, 2, MPI_COMM_WORLD);
	// MPI_Send(align2has, align2has.length(), MPI_CHAR, root, 3, MPI_COMM_WORLD);
}

// function to find out the minimum penalty
// return the minimum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans)
{

	int i, j; // intialising variables

	int m = x.length(); // length of gene1
	int n = y.length(); // length of gene2

	// table for storing optimal substructure answers
	int **dp = new2d(m + 1, n + 1);
	size_t size = m + 1;
	size *= n + 1;
	memset(dp[0], 0, size);

	// intialising the table
	for (i = 0; i <= m; i++)
	{
		dp[i][0] = i * pgap;
	}
	for (i = 0; i <= n; i++)
	{
		dp[0][i] = i * pgap;
	}

	// calcuting the minimum penalty
	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			if (x[i - 1] == y[j - 1])
			{
				dp[i][j] = dp[i - 1][j - 1];
			}
			else
			{
				dp[i][j] = min3(dp[i - 1][j - 1] + pxy,
								dp[i - 1][j] + pgap,
								dp[i][j - 1] + pgap);
			}
		}
	}

	// Reconstructing the solution
	int l = n + m; // maximum possible length

	i = m;
	j = n;

	int xpos = l;
	int ypos = l;

	while (!(i == 0 || j == 0))
	{
		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--;
			j--;
		}
		else if (dp[i - 1][j - 1] + pxy == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--;
			j--;
		}
		else if (dp[i - 1][j] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[i][j - 1] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)'_';
			yans[ypos--] = (int)y[j - 1];
			j--;
		}
	}
	while (xpos > 0)
	{
		if (i > 0)
			xans[xpos--] = (int)x[--i];
		else
			xans[xpos--] = (int)'_';
	}
	while (ypos > 0)
	{
		if (j > 0)
			yans[ypos--] = (int)y[--j];
		else
			yans[ypos--] = (int)'_';
	}

	int ret = dp[m][n];

	delete[] dp[0];
	delete[] dp;
}
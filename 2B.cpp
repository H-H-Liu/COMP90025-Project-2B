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
#include <omp.h>
#include <algorithm>

using namespace std;

std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap, int *penalties);
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);
void do_MPI_task(int rank);

/*
Examples of sha512 which returns a std::string
sw::sha512::calculate("SHA512 of std::string") // hash of a string, or
sw::sha512::file(path) // hash of a file specified by its path, or
sw::sha512::calculate(&data, sizeof(data)) // hash of any block of data
*/

// Return current time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;

// Driver code
int main(int argc, char **argv){
	int rank;
	int prov;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov); // Please note, this line and the one above must be adde for threadding to work as intended
	MPI_Comm_rank(comm, &rank);
	if(rank==root){
		int misMatchPenalty;
		int gapPenalty;
		int k;
		std::cin >> misMatchPenalty;
		std::cin >> gapPenalty;
		std::cin >> k;	
		std::string genes[k];
		for(int i=0;i<k;i++) std::cin >> genes[i];

		int numPairs= k*(k-1)/2;

		int penalties[numPairs];
		
		uint64_t start = GetTimeStamp ();

		// return all the penalties and the hash of all allignments
		std::string alignmentHash = getMinimumPenalties(genes,
			k,misMatchPenalty, gapPenalty,
			penalties);
		
		// print the time taken to do the computation
		printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
		
		// print the alginment hash
		std::cout<<alignmentHash<<std::endl;

		for(int i=0;i<numPairs;i++){
			std::cout<<penalties[i] << " ";
		}
		std::cout << std::endl;
	} else {
		// do stuff for MPI tasks that are not rank==root
		do_MPI_task(rank);
	}
	MPI_Finalize();
	return 0;
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/
class Array2D
{
private:
	size_t m_index;
	size_t xaxis;
	size_t yaxis;
	int *data;

public:
	Array2D(size_t x, size_t y) : m_index(x * y),
								  xaxis(x),
								  yaxis(y)
	{
		this->data = new int[this->m_index];
		//	std::fill(this->data, this->data+this->m_index, 0);
		// memset(this->data, 0, sizeof(int)*this->m_index);
		//	memset(this->data, 0, this->m_index);
		if (!this->data)
		{
			std::cerr << "Array2D: allocation failed, size: " << this->m_index << std::endl;
			exit(1);
		}
	}
	inline void re_alloc(size_t x, size_t y)
	{
		delete[] this->data;
		this->m_index = x * y;
		this->xaxis = x;
		this->yaxis = y;
		this->data = new int[this->m_index];
		// printf("index = %d\n",this->m_index);
		// memset(this->data, 0, this->m_index);
		// std::fill(this->data, this->data+this->m_index, 0);
		// memset(this->data, 0, sizeof(int)*this->m_index);
		// memset(this->data, 0, sizeof(int)*this->m_index);
		if (!this->data)
		{
			std::cerr << "Array2D: allocation failed, size: " << this->m_index << std::endl;
			exit(1);
		}
	}

	inline int &operator()(int x, int y)
	{

		return this->data[x * this->yaxis + y];
	}
	inline void set0(size_t x, size_t y)
	{
		this->m_index = x * y;
		this->xaxis = x;
		this->yaxis = y;
		// printf("index = %d, x=\n",this->m_index);
		// memset(this->data, 0, this->m_index);
		// std::fill(this->data, this->data+this->m_index, 0);
	}

	inline void free()
	{
		delete[] this->data;
	}
};

int getMinimumPenaltyOMP(char *x, char *y, int m, int n, int pxy, int pgap, char xans[], char yans[], Array2D array);

struct data
{
	int ID;
	int cuml_size;
	int string_id1;
	int string_id2;
	int l1;
	int l2;
	int order;
};

bool compare1(data const &lhs, data const &rhs)
{
	return lhs.cuml_size > rhs.cuml_size;
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
								int *penalties)

{
	int root = 0;
	// printf("HELLO WORLD\n");
	double tmp = (k - 1) * k;
	int size = tmp / 2 + 1;
	//	int kek = floor(50);
	// printf("size = %d\n",size);
	data *database = new data[size];

	int probNum = 0;
	for (int i = 1; i < k; i++)
	{
		for (int j = 0; j < i; j++)
		{
			// printf("YES\n");
			std::string gene1 = genes[i];
			std::string gene2 = genes[j];

			int m = gene1.length(); // length of gene1
			int n = gene2.length(); // length of gene2
			int l = m + n;
			//		std::cout << m << " " << n << std::endl;
			database[probNum].ID = probNum;
			//	std::string s1;
			//	std::string s2;

			database[probNum].string_id1 = i;
			database[probNum].string_id2 = j;
			database[probNum].cuml_size = l;
			database[probNum].l1 = m;
			database[probNum].l2 = n;
			database[probNum].order = probNum;
			probNum++;
		}
	}

	std::sort(database, database + probNum, &compare1);

	std::string gene;
	gene = genes[0];

	for (int i = 1; i < k; i++)
	{
		gene.append(" ").append(genes[i]);
	}
	// std::cout << "Gene = " << gene << std :: endl;

	int proc;
	MPI_Comm_size(MPI_COMM_WORLD, &proc);
	// int highest_rank = size - 1;

	int pivot = 1;
	int index = 0;

	int items[proc + 5] = {0};
	//  memset(items,0,proc+5);
	// items = {0};
	std::string alignmentHashs[probNum];

	// std::cout << proc << std::endl;
	for (int p = 1; p < proc; p++)
	{
		int val;
		MPI_Recv(&val, 1, MPI_INT, p, 420, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	int b1[2];
	b1[0] = k;
	b1[1] = gene.length();
	MPI_Bcast(b1, 2, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast(const_cast<char *>(gene.data()), b1[1], MPI_CHAR, 0, MPI_COMM_WORLD);

	//	printf("ALL CLEAR\n");
#pragma omp parallel
	{
		omp_set_num_threads(proc);

		// printf("Threads = %d\n",omp_get_num_threads());
#pragma omp master
		{
			// printf("start parallel\n\n\n");
			while (true)
			{
				if (index >= probNum)
				{
					break;
				}

				if (pivot >= proc)
				{
					pivot = 1;
				}
				int red;
#pragma omp atomic read
				red = items[pivot];

				//	printf("OH NO\n");
				if (red == 0)
				{
#pragma omp atomic update
					items[pivot]++;
#pragma omp task firstprivate(index, pivot)
					{

						int sta = 0;
						MPI_Send(&sta, 1, MPI_INT, pivot, 69, MPI_COMM_WORLD);

						int len[6];
						len[0] = database[index].l1;
						len[1] = database[index].l2;
						len[2] = pxy;
						len[3] = pgap;
						len[4] = database[index].string_id1;
						len[5] = database[index].string_id2;

						MPI_Send(len, 6, MPI_INT, pivot, 0, MPI_COMM_WORLD);

						MPI_Status sta1;
						MPI_Status sta3;
						int result[2];
						MPI_Recv(result, 2, MPI_INT, pivot, 3, MPI_COMM_WORLD, &sta1);
						char buf1[result[1]];
						MPI_Recv(buf1, result[1], MPI_CHAR, pivot, 4, MPI_COMM_WORLD, &sta3);

						penalties[database[index].order] = result[0];

						std::string a(buf1, 128);

						alignmentHashs[database[index].order] = a;

#pragma omp atomic update
						items[pivot]--;
					}
					index++;
				}
				pivot++;
			}
		}
	}
	//   delete[] items;
	//	printf("COMPLETE\n");

#pragma omp parallel for
	for (int i = 1; i < proc; i++)
	{
		int sta = 1;
		MPI_Send(&sta, 1, MPI_INT, i, 69, MPI_COMM_WORLD);
	}
	// printf("Stopping all processes\n");

	//	printf("END\n");
	// Todo, see if you can parallise this

	std::string alignmentHash = "";

	for (int i = 0; i < probNum; i++)
	{
		//	std::cout << "Problem: " << i << " " << alignmentHashs[i] << " " << alignmentHashs[i].length() << std::endl;
		alignmentHash = sw::sha512::calculate(alignmentHash.append(alignmentHashs[i]));
	}

	return alignmentHash;
}

// called for all tasks with rank!=root
// do stuff for each MPI task based on rank
void do_MPI_task(int rank)
{
	int root = 0;

	int conf;
	conf = 1;
	MPI_Send(&conf, 1, MPI_INT, root, 420, MPI_COMM_WORLD);
	//printf("START\n");

	int len = 89 * 2 + 2;
	Array2D array(89 + 1, 89 + 1);

	// printf("k = %d\n",k);

	//	printf("Iterator = %d,\n",i);
	int b2[2];
	MPI_Bcast(b2, 2, MPI_INT, root, MPI_COMM_WORLD);
	char *hol = new char[b2[1] + 1];
	MPI_Bcast(hol, b2[1], MPI_CHAR, root, MPI_COMM_WORLD);
	// std::cout << "rank= " << rank << " b2 " << b2[0] << " " << b2[1] << std::endl;

	char **genes = new char *[b2[0]];
	hol[b2[1]] = '\0';
	// std::cout << "Len = " << b2[1] << " String = " << hol << std::endl;
	//	genes[i] = hol;
	// std::cout << "i = " << i << " Genes = " << genes[i] << std::endl;
//	MPI_Barrier(MPI_COMM_WORLD);

	int it = 0;
	char *token = strtok(hol, " ");
	// loop through the string to extract all other tokens
	while (token != NULL)
	{
		//   printf( "rank = %d, it = %d\n", rank, it); //printing each token
		int strl = strlen(token);
		// printf(" rank = %d, it = %d, strl = %d\n",rank,it, strl);
		char *tmp = new char[strl + 1];
		strcpy(tmp, token);
		tmp[strl] = '\0';
		genes[it] = tmp;
		// printf("rank = %d Token = %s, len = %d\n",rank,token,strl);
		it++;
		token = strtok(NULL, " ");
		// printf("rank = %d, it = %d, next token = %s\n",rank, it,token);
	}
	// printf("COMPLETE\n");
	// MPI_Barrier(MPI_COMM_WORLD);
	// 	printf( "rank = %d, Complete\n", rank);

	delete[] hol;

	// std::cout << "rank = " << rank << " Genes = " << genes[5] << std::endl;
	// printf("Rank = %d, complete",rank);
	// MPI_Barrier(MPI_COMM_WORLD);
	// exit(0);
	// TODO

	while (true)
	{
		int status;
		MPI_Status sta69;

		MPI_Recv(&status, 1, MPI_INT, root, 69, MPI_COMM_WORLD, &sta69);

		if (status == 1)
		{

			break;
		}

		int ints[6];
		int l1, l2, pxy, pgap, id1, id2;
		MPI_Status sta1, sta2, sta3, sta4, sta5, sta6;

		MPI_Recv(ints, 6, MPI_INT, root, 0, MPI_COMM_WORLD, &sta1);

		l1 = ints[0];
		l2 = ints[1];
		pxy = ints[2];
		pgap = ints[3];
		id1 = ints[4];
		id2 = ints[5];

		// MPI_Irecv(buf1, l1, MPI::CHAR, root, 1, MPI_COMM_WORLD, &req[0]);
		// MPI_Irecv(buf2, l2, MPI::CHAR, root, 2, MPI_COMM_WORLD, &req[1]);
		// MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[])

		char xans[l2 + l1 + 1];
		char yans[l2 + l1 + 1];

		if (len < l1 + l2 + 2)
		{
			array.re_alloc(l1 + 1, l2 + 1);
			len = l1 + l2 + 2;
		}
		else
		{
			array.set0(l1 + 1, l2 + 1);
		}

		//	MPI_Waitall(2,req, stats);

		// printf("rank = %d START, l1 = %d, l2 = %d, len = %d, id1 = %d, id2 = %d ",rank ,l1+1,l2+1,len,id1,id2);
	//	std::cout << "HELLO " << rank << std::endl;
		//std::cout << genes[id1] << " " << genes[id2] << " " << rank << std::endl;
		int ack = getMinimumPenaltyOMP(genes[id1], genes[id2], l1, l2, pxy, pgap, xans, yans, array);
	//	printf("rank = %d FINISH \n", rank);
		int l = l1 + l2;

		int id = 1;
		int a;
		for (a = l; a >= 1; a--)
		{
			if ( yans[a] == '_' &&  xans[a] == '_')
			{
				id = a + 1;
				break;
			}
		}

		std::string align1 = "";
		std::string align2 = "";
		for (a = id; a <= l; a++)
		{
			align1.append(1,  xans[a]);
		}
		for (a = id; a <= l; a++)
		{
			align2.append(1,  yans[a]);
		}
		std::string align1hash = sw::sha512::calculate(align1);
		std::string align2hash = sw::sha512::calculate(align2);
		std::string problemhash = sw::sha512::calculate(align1hash.append(align2hash));

		int data[2];
		data[0] = ack;
		data[1] = problemhash.length();
		//	std::cout << "Final Hash" << problemhash << "Length" << problemhash.length() << std::endl;
		MPI_Send(data, 2, MPI_INT, root, 3, MPI_COMM_WORLD);
		MPI_Send(problemhash.c_str(), problemhash.length(), MPI_CHAR, root, 4, MPI_COMM_WORLD);
	}

	return;
}

// function to find out the minimum penalty
// return the minimum penalty and put the aligned sequences in xans and yans
int getMinimumPenaltyOMP(char *x, char *y, int m, int n, int pxy, int pgap, char xans[], char yans[], Array2D array)
{

	int ret = -1;
	// Specialised memoery allocation function
//	Array2D array(m + 1, n + 1);
// Array2DA arrayY(n + 1, m + 1);
// Array2DA arrayXY(m+1,n+1);

// printf("STARTED\n");
// omp_set_num_threads(10);

// intialising the table
	//#pragma omp for schedule(static)
	for (int i = 0; i <= m; i++)
	{
		array(i, 0) = i * pgap;
		//	arrayY(0, i) = i * pgap;
		//	arrayXY(i, 0) = i * pgap;
	}
	// printf("INIT COMPELTE\n");

	//#pragma omp for schedule(static)
	for (int i = 0; i <= n; i++)
	{
		array(0, i) = i * pgap;
		//	arrayY(i, 0) = i * pgap;
		//	arrayXY(0, i) = i * pgap;
	}

	// printf("MALLOC COMPLETE\n");
	omp_set_num_threads(omp_get_max_threads());
	//	printf("max threads = %d\n",omp_get_max_threads());
	int ready[8 + omp_get_max_threads()] = {0};
	ready[0] = m;

	int done = 0;
	int i, j, id, yfirst, ylast, t;
	#pragma omp parallel default(none) private(i, j, id, yfirst, ylast, t) shared(x, y, array, pxy, pgap, n, m, ready, done)
	{
		t = omp_get_num_threads();
		id = omp_get_thread_num();
		// Precalcualte first and last to scrape for time?
		yfirst = floor(((double)id * n) / ((double)t)) + 1;
		ylast = floor(((double)(id + 1) * n) / ((double)t));

		for (i = 1; i <= m; i = i + 2)
		{
			// Read write conflicta
			int red;
			#pragma omp atomic read
			red = ready[id];

			while (red <= i - 1)
			{
// printf("done = %d id = %d\n",done,id);
				#pragma omp atomic read
				red = ready[id];
				// printf("id = %d read = %d\n",id,red);
			}

			// printf("i = %d id = %d, array = %d\n",i,id,array(i,yfirst-1));
			for (j = yfirst; j <= ylast; j++)
			{
				if (x[i - 1] == y[j - 1])
				{
					array(i, j) = array(i - 1, j - 1);
				}
				else
				{
					array(i, j) = min3(array(i - 1, j - 1) + pxy, array(i - 1, j) + pgap, array(i, j - 1) + pgap);
				}

				if (i == m)
				{
					continue;
				}
				if (x[i] == y[j - 1])
				{
					array(i + 1, j) = array(i, j - 1);
				}
				else
				{
					array(i + 1, j) = min3(array(i, j - 1) + pxy, array(i, j) + pgap, array(i + 1, j - 1) + pgap);
				}
			}
			#pragma omp atomic write
			ready[id + 1] = i;
			//	#pragma omp barrier
			// 	printf(" updated id = %d red = %d \n",id+1,ready[id+1]);
		}

		// printf("exit m = %d, id = %d,  red[%d] = %d,\n",m,id,id+1,ready[id+1]);
	}

	//	printf("DONE\n");

	// Reconstructing the solution
	int l = n + m; // maximum possible length

	i = m;
	j = n;

	int xpos = l;
	int ypos = l;

	// Back tracking is virtually uncahnged hence uncommented.
	while (!(i == 0 || j == 0))
	{
		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] =  x[i - 1];
			yans[ypos--] =  y[j - 1];
			i--;
			j--;
		}
		else if (array(i - 1, j - 1) + pxy == array(i, j))
		{
			xans[xpos--] =  x[i - 1];
			yans[ypos--] =  y[j - 1];
			i--;
			j--;
		}
		else if (array(i - 1, j) + pgap == array(i, j))
		{
			xans[xpos--] =  x[i - 1];
			yans[ypos--] =  '_';
			i--;
		}
		else if (array(i, j - 1) + pgap == array(i, j))
		{
			xans[xpos--] =  '_';
			yans[ypos--] =  y[j - 1];
			j--;
		}
	}
	while (xpos > 0)
	{
		if (i > 0)
			xans[xpos--] =  x[--i];
		else
			xans[xpos--] =  '_';
			
	}
	while (ypos > 0)
	{
		if (j > 0)
			yans[ypos--] =  y[--j];
		else
			yans[ypos--] =  '_';
	}
	ret = array(m, n);
	// array.free();

	return ret;

}

// mpicxx -O3 -fopenmp -o haohail-seqalkway haohail-seqalkway.cpp

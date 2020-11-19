# COMP90025-Project-2B
Sequence Alignment with MPI and Open MP

What is sequence Alignment? 
See https://en.wikipedia.org/wiki/Sequence_alignment 

To compile please ensure you have MPI installed and a OPENMP supporting C Compiler.

Then do this: '''mpicxx -O3 -fopenmp -o 2B 2B.cpp'''

Genome Sequence Data is quite freely availiable but as a demo on how to use this program. 


Define a text file with the following properties: 

3
2
3
AGGGCT
AGGCAAAA
GGGCT

Where the first number is the number of sequences we want to match
The second number is the penalty weight
The third number is the gap weight

Then we do 
'''$ mpirun ./2B < [textfile name here]'''


Expected Output

Time: 132
us602d0f604e8fb908195d53e681094f7d063c4168a33a18f32b4ca3d29f27073a486dca2ab98aab9eb47f5c407b5c59b8e6c0fa8ef4d07d131b8d6a66a37a065f
5 4 9


Final result: 3rd out of 22 students. approximately few seconds slower than the best result.

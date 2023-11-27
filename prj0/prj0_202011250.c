#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 202011250 고정현

//functions for convenience
double** allocateMemory(int m, int n); 
void releaseMemory(double** A, int m);
void printMatrix(double** A, int m, int n, char name[]);

//functions to implement in prj0
//transposeMatrix, normalizeVector, calculateLength
//scaleMatrix, multiplyTwoMatrices, addTwoMatrices
double** transposeMatrix(double** A, int m, int n);
double** normalizeVector(double** v, int n);
double calculateLength(double** v, int n);
void scaleMatrix(double** A, int m, int n, double c);
double** multiplyTwoMatrices(double** A, int m, int n, double** B, int l, int k);
double** addTwoMatrices(double** A, int m, int n, double** B, int l, int k);

int main() {
	// double** A;
	// double** B;
	// double** v; //vector
	// double** w;	//vector
	// int m = 5; 	//number of rows
	// int n = 10;	//number of columns

	// //Test transposeMatrix
	// A = allocateMemory(m,n);
	// for (int i = 0; i < m; i++)
	// 	for (int j = 0; j < n; j++)
	// 		A[i][j] = (double) i*j;
	// printMatrix(A,m,n,"A");

	// B = transposeMatrix(A,m,n);
	// printMatrix(B,n,m,"B");

	// //Test normalizeVector
	// v = allocateMemory(m,1);
	// for (int i = 0; i < m; i++)
	// 	v[i][0] = i;
	// w = normalizeVector(v,m);
	// printMatrix(v,m,1,"v");
	// printMatrix(w,m,1,"v_normalized");

	// //release all the memory allocated
	// releaseMemory(A, m);
	// releaseMemory(B, n);
	// releaseMemory(v, m);
	// releaseMemory(w, m);

	printf("202011250\n");
	
	double** A;
	double** H1;
	double** v;
	double** w;
	double** H;
	double** HT;
	double** B;
	double** C;
	int n = 2;

	A = allocateMemory(n, n);
	A[0][0] = 1;
	A[0][1] = 2;
	A[1][0] = 3;
	A[1][1] = 4;
	printMatrix(A, n, n, "A");

	H1 = allocateMemory(n, n);
	H1[0][0] = 1;
	H1[0][1] = 1;
	H1[1][0] = 1;
	H1[1][1] = -1;
	printMatrix(H1, n, n, "H~");
	
	v = allocateMemory(n, 1);
	for (int i=0 ; i<n ; i++)
		v[i][0] = H1[i][0];
	w = allocateMemory(n, 1);
	for (int i=0 ; i<n ; i++)
		w[i][0] = H1[i][1];

	v = normalizeVector(v, n);
	w = normalizeVector(w, n);

	H = allocateMemory(n, n);
	for (int i=0 ; i<n ; i++)
		H[i][0] = v[i][0];
	for (int i=0 ; i<n ; i++)
		H[i][1] = w[i][0];
	printMatrix(H, n, n, "H");

	HT = allocateMemory(n, n);
	HT = transposeMatrix(H, n, n);
	
	B = allocateMemory(n, n);
	B = multiplyTwoMatrices(HT, n, n, A, n, n);
	B = multiplyTwoMatrices(B, n, n, H, n, n);
	printMatrix(B, n, n, "B = HTAH");

	C = allocateMemory(n, n);
	C = multiplyTwoMatrices(H, n, n, B, n, n);
	C = multiplyTwoMatrices(C, n, n, HT, n, n);
	printMatrix(C, n, n, "C = HBHT");

	releaseMemory(A, n);
	releaseMemory(H1, n);
	releaseMemory(H, n);
	releaseMemory(HT, n);
	releaseMemory(B, n);
	releaseMemory(C, n);
	releaseMemory(v, n);
	releaseMemory(w, n);

	return 0;
}

//functions for convenience
double** allocateMemory(int m, int n) {
	double** A;
	A = (double**) malloc(sizeof(double*) * m);
	for (int i = 0; i < m; i++) {
		A[i] = (double*) malloc(sizeof(double) * n);
	}
	return A;
}


void releaseMemory(double** A, int m) {
	for (int i = 0; i < m; i++)
		free(A[i]);
	free(A);
}

void printMatrix(double** A, int m, int n, char name[]) {
	printf("\n%s = \n", name);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%lf ", A[i][j]);
		printf("\n");
	}
}

//functions to implement in prj0 
double** transposeMatrix(double **A, int m, int n) {
	double** B = allocateMemory(n, m);

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			B[j][i] = A[i][j];	
	
	return B;
}	

double** normalizeVector(double** v, int m) {
	double** w;
	double len = 0.0;

	for (int i = 0; i < m; i++)
		len += v[i][0]*v[i][0];	
	len = sqrt(len);

	w = allocateMemory(m,1);
	for (int i = 0; i < m; i++)
		w[i][0] = v[i][0]/len;
	
	return w;
}

double calculateLength(double** v, int n)
{
	double len = 0.0;

	for (int i=0 ; i<n ; i++) 
		len += v[i][0] * v[i][0];

	return sqrt(len);
}

void scaleMatrix(double** A, int m, int n, double c)
{
	for (int i=0 ; i<m ; i++)
		for (int j=0 ; j<n ; j++)
			A[i][j] *= c;
}

double** multiplyTwoMatrices(double** A, int m, int n, double** B, int l, int k)
{
	if (n!=l)
		return NULL;

	double** result = allocateMemory(m, k);

	for (int i=0 ; i<m ; i++)
		for (int j=0 ; j<k ; j++)
		{
			result[i][j] = 0.0;
			for (int o=0 ; o<n ; o++)
				result[i][j] += A[i][o]*B[o][j];
		}

	return result;
}

double** addTwoMatrices(double** A, int m, int n, double** B, int l, int k)
{
	if (m!=l && n!=k)
		return NULL;
		
	double** result = allocateMemory(m, n);

	for (int i=0 ; i<m ; i++)
		for (int j=0 ; j<n ; j++)
			result[i][j] = A[i][j]+B[i][j];

	return result;
}
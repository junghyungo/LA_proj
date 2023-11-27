// 202011250 고정현

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>

//비트맵 헤더를 한묶음으로
typedef struct tagBITMAPHEADER {
	BITMAPFILEHEADER bf;
	BITMAPINFOHEADER bi;
	RGBQUAD hRGB[256]; //이 코드에서는 필요없음 (8bit에만 필요)
}BITMAPHEADER;

//비트맵을 읽어와서 화소정보의 포인터를 리턴
BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int*imgHeight, char* filename);

//비트맵 파일 쓰기
void writeBitmapFile(int bytesPerPixel, BITMAPHEADER outputHeader, BYTE* output, int imgSize, char* filename);

double** allocateMemory(int m, int n);
void releaseMemory(double** A, int m);
double** constructHaarMatrixRecursive(int n);
double** applyKroneckerProduct(double** A, int n, int a, int b);
double** concatenateTwoMatrices(double** hl, double** hr, int n);
double** normalize(double** A, int n);
double** transposeMatrix(double **A, int m, int n);
double** multiplication(double **A, double**B, int m, int p, int n);
void printMatrix(double** A, int m, int n, char name[]);
double** constructIdentity(int k);

int main() {
	
	/*******************************************************************/
	/*************************** Read image  ***************************/
	/*******************************************************************/
	BITMAPHEADER originalHeader;	//비트맵의 헤더부분을 파일에서 읽어 저장할 구조체
	BITMAPHEADER outputHeader;		//변형을 가한 헤더부분을 저장할 구조체

	int imgSize, imgWidth, imgHeight;					//이미지의 크기를 저장할 변수
	int bytesPerPixel = 3;			//number of bytes per pixel (1 byte for R,G,B respectively)

	char imgLena[] = "image_lena_24bit.bmp";
	char img1[] = "tiger.bmp";
	char img2[] = "MickeyMouse.bmp";

	BYTE* image = loadBitmapFile(bytesPerPixel, &originalHeader, &imgWidth, &imgHeight, imgLena); //비트맵파일을 읽어 화소정보를 저장 (불러들이는 이미지는 .c와 같은 폴더에 저장)
	if (image == NULL) return 0;

	imgSize = imgWidth * imgHeight; // total number of pixels

	BYTE* output = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);				//결과값을 저장할 포인터 선언 및 메모리 할당
	outputHeader = originalHeader;	

	/*******************************************************************/
	/************************ Perform HWT/IHWT *************************/
	/*******************************************************************/
	// 1(a)
	//이미지 행렬 A 구성 (RGB값이 있으므로 픽셀당 값 하나씩만 읽어서 imgWidth x imgHeight 행렬 구성)
	double** A = allocateMemory(imgHeight, imgWidth);

	for (int i=0 ; i<imgHeight ; i++)
		for (int j=0 ; j<imgWidth ; j++)
			A[i][j] = image[(i*imgWidth+j)*bytesPerPixel];

	// 1(b)
	//Haar matrix H 구성 (orthonormal column을 갖도록 구성)
	int n = imgHeight; //이미지가 정사각형(Height==Width)이라고 가정; n = 2^t,t=0,1,2,...

	double** Htmp = constructHaarMatrixRecursive(n);
	double** H = normalize(Htmp, n);

	// 1(c)
	//HWT 수행: 행렬곱 B = H'*A*H
	double** HT = transposeMatrix(H, n, n);

	double** Btmp = multiplication(HT, A, n, n, n);
	double** B = multiplication(Btmp, H, n, n, n);

	// 1(d)
	//행렬 B 자르기: B의 upper left corner(subsquare matrix)를 잘라 Bhat에 저장
	int k = n; 	// n=512

	double** Bhat = allocateMemory(n, n);

	for (int i=0 ; i<n ; i++) {
        for (int j=0 ; j<n ; j++) {
			Bhat[i][j] = 0.0;
            if (i<k && j<k) {
                Bhat[i][j] = B[i][j];
            }
        }
    }
	
	// 1(e)
	//IHWT 수행: Ahat = H*Bhat*H'
	double** Ahattmp = multiplication(H, Bhat, n, n, n);
	double** Ahat = multiplication(Ahattmp, HT, n, n, n);

	// 1(f)
	/*******************************************************************/
	/******************* Write reconstructed image  ********************/
	/*******************************************************************/
	//Ahat을 이용해서 위의 image와 같은 형식이 되도록 구성 (즉, Ahat = [a b;c d]면 [a a a b b b c c c d d d]를 만들어야 함)
	BYTE* Are = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	int index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are[index] = (BYTE)Ahat[i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are, imgSize, "output_k_512.bmp");

	// 3(c)
	double** Hl = allocateMemory(n/2, n);
	double** Hh = allocateMemory(n/2, n);

	for (int i=0 ; i<n/2 ; i++) {
		for (int j=0 ; j<n ; j++) {
			Hl[i][j] = HT[i][j];
			Hh[i][j] = HT[i+n/2][j];
		}
	}

	double** HlT = transposeMatrix(Hl, n/2, n);
	double** HhT = transposeMatrix(Hh, n/2, n);

	// 1. HlT*Hl*A*HlT*Hl
	double** tmp1_1 = multiplication(HlT, Hl, n, n/2, n);
	double** tmp1_2 = multiplication(tmp1_1, A, n, n, n);
	double** tmp1_3 = multiplication(tmp1_2, HlT, n, n, n/2);
	double** HlTHlAHlTHl = multiplication(tmp1_3, Hl, n, n/2, n);
	BYTE* Are1 = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are1[index] = (BYTE)HlTHlAHlTHl[i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are1, imgSize, "output_1.bmp");

	// 2. HlT*Hl*A*HhT*Hh
	double** tmp2_1 = multiplication(HlT, Hl, n, n/2, n);
	double** tmp2_2 = multiplication(tmp2_1, A, n, n, n);
	double** tmp2_3 = multiplication(tmp2_2, HhT, n, n, n/2);
	double** HlTHlAHhTHh = multiplication(tmp2_3, Hh, n, n/2, n);
	BYTE* Are2 = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are2[index] = (BYTE)HlTHlAHhTHh[i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are2, imgSize, "output_2.bmp");

	// 3. HhT*Hh*A*HlT*Hl
	double** tmp3_1 = multiplication(HhT, Hh, n, n/2, n);
	double** tmp3_2 = multiplication(tmp3_1, A, n, n, n);
	double** tmp3_3 = multiplication(tmp3_2, HlT, n, n, n/2);
	double** HhTHhAHlTHl = multiplication(tmp3_3, Hl, n, n/2, n);
	BYTE* Are3 = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are3[index] = (BYTE)HhTHhAHlTHl[i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are3, imgSize, "output_3.bmp");

	// 4. HhT*Hh*A*HhT*Hh
	double** tmp4_1 = multiplication(HhT, Hh, n, n/2, n);
	double** tmp4_2 = multiplication(tmp4_1, A, n, n, n);
	double** tmp4_3 = multiplication(tmp4_2, HhT, n, n, n/2);
	double** HhTHhAHhTHh = multiplication(tmp4_3, Hh, n, n/2, n);
	BYTE* Are4 = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are4[index] = (BYTE)HhTHhAHhTHh[i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are4, imgSize, "output_4.bmp");

	// 3(d)
	double** Hll = allocateMemory(n/4, n);
	double** Hlh = allocateMemory(n/4, n);

	for (int i=0 ; i<n/4 ; i++) {
		for (int j=0 ; j<n ; j++) {
			Hll[i][j] = HT[i][j];
			Hlh[i][j] = HT[i+n/4][j];
		}
	}

	double** HllT = transposeMatrix(Hll, n/4, n);
	double** HlhT = transposeMatrix(Hlh, n/4, n);

	// 5. HllT*Hll*A*HllT*Hll
	double** tmp5_1 = multiplication(HllT, Hll, n, n/4, n);
	double** tmp5_2 = multiplication(tmp5_1, A, n, n, n);
	double** tmp5_3 = multiplication(tmp5_2, HllT, n, n, n/4);
	double** HllTHllAHllTHll = multiplication(tmp5_3, Hll, n, n/4, n);
	BYTE* Are5 = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are5[index] = (BYTE)HllTHllAHllTHll[i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are5, imgSize, "output_5.bmp");

	// 6. HllT*Hll*A*HlhT*Hlh
	double** tmp6_1 = multiplication(HllT, Hll, n, n/4, n);
	double** tmp6_2 = multiplication(tmp6_1, A, n, n, n);
	double** tmp6_3 = multiplication(tmp6_2, HlhT, n, n, n/4);
	double** HllTHllAHlhTHlh = multiplication(tmp6_3, Hlh, n, n/4, n);
	BYTE* Are6 = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are6[index] = (BYTE)HllTHllAHlhTHlh[i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are6, imgSize, "output_6.bmp");

	// 7. HlhT*Hlh*A*HllT*Hll
	double** tmp7_1 = multiplication(HlhT, Hlh, n, n/4, n);
	double** tmp7_2 = multiplication(tmp7_1, A, n, n, n);
	double** tmp7_3 = multiplication(tmp7_2, HllT, n, n, n/4);
	double** HlhTHlhAHllTHll = multiplication(tmp7_3, Hll, n, n/4, n);
	BYTE* Are7 = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are7[index] = (BYTE)HlhTHlhAHllTHll [i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are7, imgSize, "output_7.bmp");

	// 8. HlhT*Hlh*A*HlhT*Hlh
	double** tmp8_1 = multiplication(HlhT, Hlh, n, n/4, n);
	double** tmp8_2 = multiplication(tmp8_1, A, n, n, n);
	double** tmp8_3 = multiplication(tmp8_2, HlhT, n, n, n/4);
	double** HlhTHlhAHlhTHlh = multiplication(tmp8_3, Hlh, n, n/4, n);
	BYTE* Are8 = (BYTE*) malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	index = 0;
	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			for (int k=0 ; k<bytesPerPixel ; k++) {
				Are8[index] = (BYTE)HlhTHlhAHlhTHlh[i][j];
				index++;
			}
		}
	}
	writeBitmapFile(bytesPerPixel, outputHeader, Are8, imgSize, "output_8.bmp");


	free(image);
	free(output);
	for (int i = 0; i < imgHeight; i++)
		free(A[i]);
	free(A);
	free(Are);

	return 0;
}

BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int* imgHeight, char* filename)
{
	FILE* fp = fopen(filename, "rb");	//파일을 이진읽기모드로 열기
	if (fp == NULL)
	{
		printf("파일로딩에 실패했습니다.\n");	//fopen에 실패하면 NULL값을 리턴
		return NULL;
	}
	else
	{
		fread(&bitmapHeader->bf, sizeof(BITMAPFILEHEADER), 1, fp);	//비트맵파일헤더 읽기
		fread(&bitmapHeader->bi, sizeof(BITMAPINFOHEADER), 1, fp);	//비트맵인포헤더 읽기
		//fread(&bitmapHeader->hRGB, sizeof(RGBQUAD), 256, fp);	//색상팔렛트 읽기 (24bitmap 에서는 존재하지 않음)

		*imgWidth = bitmapHeader->bi.biWidth;
		*imgHeight = bitmapHeader->bi.biHeight;
		int imgSizeTemp = (*imgWidth) * (*imgHeight);	// 이미지 사이즈를 상위 변수에 할당

		printf("Size of image: Width %d   Height %d\n", bitmapHeader->bi.biWidth, bitmapHeader->bi.biHeight);
		BYTE* image = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSizeTemp);	//이미지크기만큼 메모리할당

		fread(image, bytesPerPixel*sizeof(BYTE), imgSizeTemp, fp);//이미지 크기만큼 파일에서 읽어오기

		fclose(fp);
		return image;
	}
}
void writeBitmapFile(int bytesPerPixel, BITMAPHEADER outputHeader, BYTE* output, int imgSize, char* filename)
{
	FILE* fp = fopen(filename, "wb");

	fwrite(&outputHeader.bf, sizeof(BITMAPFILEHEADER), 1, fp);
	fwrite(&outputHeader.bi, sizeof(BITMAPINFOHEADER), 1, fp);
	//fwrite(&outputHeader.hRGB, sizeof(RGBQUAD), 256, fp); //not needed for 24bitmap
	fwrite(output, bytesPerPixel*sizeof(BYTE), imgSize, fp);
	fclose(fp);
}

double** allocateMemory(int m, int n) {
	double** A = (double**)malloc(sizeof(double*) * m);
	for (int i = 0; i < m; i++)
		A[i] = (double*)malloc(sizeof(double) * n);
	return A;
}
void releaseMemory(double** A, int m) {
	for (int i = 0; i < m; i++)
		free(A[i]);
	free(A);
}

double** constructHaarMatrixRecursive(int n) {
	double** h;

	if (n > 2)
		h = constructHaarMatrixRecursive(n/2);
	else {
		h = allocateMemory(2,2);
		//H = [1 1; 1 -1]
		h[0][0] = 1; 
		h[0][1] = 1; 
		h[1][0] = 1; 
		h[1][1] = -1; 
		return h; 
	}

	// construct the left half (Kronecket product of h and [1;1])
	double** Hl = applyKroneckerProduct(h, n, 1, 1);
	releaseMemory(h, n/2);

	// construct the right half (Kronecker product of I and [1;-1])
	double** I = constructIdentity(n/2);
	double** Hr = applyKroneckerProduct(I, n, 1, -1); 
	releaseMemory(I, n/2);

	// merge hl and hr
	double** H = concatenateTwoMatrices(Hl, Hr, n); //H = [Hl Hr]
	releaseMemory(Hl, n);
	releaseMemory(Hr, n);

	return H;
}
double** applyKroneckerProduct(double** A, int n, int a, int b) {
	double** h = allocateMemory(n, n/2);
	
	for (int j=0 ; j<n/2 ; j++) {
		for (int i=0 ; i<n/2 ; i++) {
			h[2*i][j] = A[i][j]*a;
			h[2*i+1][j] = A[i][j]*b;
		}
	}

	return h;
}
double** concatenateTwoMatrices(double** hl, double** hr, int n) {
	double** H = allocateMemory(n, n);

	for (int i=0 ; i<n ; i++) {
		for (int j=0 ; j<n ; j++) {
			if (j < n/2)
				H[i][j] = hl[i][j];
			else
				H[i][j] = hr[i][j-n/2];
		}
	}

	return H; 
}
double** normalize(double** A, int n) {

	double** B = allocateMemory(n, n);

    for (int j=0 ; j<n ; j++) {
        double len = 0.0;
        for (int i=0 ; i<n ; i++)
            len += A[i][j]*A[i][j];
        len = sqrt(len);

        if (len != 0.0) { // 0으로 나누는 것을 방지
            for (int i=0 ; i<n ; i++)
                B[i][j] = A[i][j]/len;
        } 
		else { 
            for (int i=0 ; i<n ; i++)
                B[i][j] = A[i][j];
        }
    }

    return B;
}
double** transposeMatrix(double **A, int m, int n) {
	double** B = allocateMemory(n, m);

	for (int i=0 ; i<m ; i++)
		for (int j=0 ; j<n ; j++)
			B[j][i] = A[i][j];	
	
	return B;
}	
double** multiplication(double **A, double**B, int m, int p, int n) {
	double** C = allocateMemory(m, n);

	for (int i=0 ; i<m ; i++) {
        for (int j=0 ; j<n ; j++) {
			C[i][j] = 0.0;
            for (int k=0 ; k<p ; k++) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
	
	return C;
}

void printMatrix(double** A, int m, int n, char name[]) {
	printf("\n%s = \n", name);

	for (int i=0 ; i<m ; i++) {
		for (int j=0 ; j<n ; j++)
			printf("%lf ", A[i][j]);
		printf("\n");
	}
}
double** constructIdentity(int k) {
 	double** I = allocateMemory(k,k);

 	for (int i=0 ; i<k ; i++) {
 		for (int j=0 ; j<k ; j++) {
 			if (i != j)
 				I[i][j] = 0.0;
 			else
 				I[i][j] = 1.0;
 		}
 	}

 	return I;
}
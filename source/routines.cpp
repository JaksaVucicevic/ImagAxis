#include <cstdio>
#include <cstdlib>
#include "routines.h"
#include "GRID.h"
#include <vector>
#include "omp.h"
#include "nrutil.h"

using namespace std;

double sign(double x)
{ //returns 1 if x positive, -1 if x negative, 0 if x=0
  if (x==0) return 0;
  else
    return (x>=0) ? +1 : -1; 
}

double sqr(double x)
{ //returns square of x
  return x*x;
} 

double cub(double x)
{ //returns square of x
  return x*x*x;
} 


int pow(int base, int exp)
{
  int res = 1;
  for (int i=0; i<exp; i++)
    res *= base;
  return res;
}

complex<double> sqr(complex<double> x)
{
 return x*x;
}

/*double abs(double x)
{
  return (x>=0) ? x : -x; 
}*/

/*double abs(complex<double> x)
{
  return  sqrt(   sqr( real(x) ) 
                + sqr( imag(x) ) ); 
}*/


//---------------------------- SPLINES -----------------------------------------------------------//

double ParabolaFrom3points(double* Y, double* X)
{
  double x1 = X[0];
  double x2 = X[1];
  double x3 = X[2];
  double y1 = Y[0];
  double y2 = Y[1];
  double y3 = Y[2];

  double c = -( (-sqr(x2)*x3*y1 + x2*sqr(x3)*y1 + sqr(x1)*x3*y2 - x1*sqr(x3)*y2 - sqr(x1)*x2*y3 + x1*sqr(x2)*y3)
                            / ( (x2 - x3)*(sqr(x1) - x1*x2 - x1*x3 + x2*x3) )
              );

  return c;
}

double ParabolaFrom3points(double* Y, double* X, double x)
{
  double x1 = X[0];
  double x2 = X[1];
  double x3 = X[2];
  double y1 = Y[0];
  double y2 = Y[1];
  double y3 = Y[2];


  double a = -( (-x2*y1 + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)
                     / ( (-x1 + x2)*(x2 - x3)*(-x1 + x3) )
              ) ;
  double b = -( (sqr(x2)*y1 - sqr(x3)*y1 - sqr(x1)*y2 + sqr(x3)*y2 + sqr(x1)*y3 -  sqr(x2)*y3)
                                 / ( (x1 - x2)*(x1 - x3)*(x2 - x3) )
              ); 
  double c = -( (-sqr(x2)*x3*y1 + x2*sqr(x3)*y1 + sqr(x1)*x3*y2 - x1*sqr(x3)*y2 - sqr(x1)*x2*y3 + x1*sqr(x2)*y3)
                            / ( (x2 - x3)*(sqr(x1) - x1*x2 - x1*x3 + x2*x3) )
              );

  return a*sqr(x) + b*x + c;

}

complex<double> ParabolaFrom3points(complex<double>* Y, double* X)
{
  double ReY[3];
  double ImY[3];  
  for(int i = 0; i<3; i++)
  {
    ReY[i] = real(Y[i]);
    ImY[i] = imag(Y[i]);
  }

  return complex<double>(ParabolaFrom3points(ReY, X), ParabolaFrom3points(ImY, X));
}

complex<double> ParabolaFrom3points(complex<double>* Y, double* X, double x)
{
  double ReY[3];
  double ImY[3];  
  for(int i = 0; i<3; i++)
  {
    ReY[i] = real(Y[i]);
    ImY[i] = imag(Y[i]);
  }

  return complex<double>(ParabolaFrom3points(ReY, X, x), ParabolaFrom3points(ImY, X, x));
}

//--- cubic---//

double CubicFrom4points(double* Y, double* X)
{
  double x1 = X[0];
  double x2 = X[1];
  double x3 = X[2];
  double x4 = X[3];
  double y1 = Y[0];
  double y2 = Y[1];
  double y3 = Y[2];
  double y4 = Y[3];

  double d = (x1 *(x1 - x3) *x3 *(x1 - x4) *(x3 - x4) *x4 *y2 + 
      x2 *sqr(x4) *(-cub(x3) *y1 + sqr(x3) *x4 *y1 + sqr(x1) *(x1 - x4) *y3) + 
      sqr(x1) *x2 *sqr(x3) *(-x1 + x3) *y4 + 
      cub(x2) *(x4 *(-sqr(x3) *y1 + x3 *x4 *y1 + x1 *(x1 - x4) *y3) + 
         x1 *x3 *(-x1 + x3) *y4) + 
      sqr(x2) *(x1* x4 *(-sqr(x1) + sqr(x4)) *y3 + cub(x3) *(x4 *y1 - x1* y4) + 
         x3 *(-cub(x4) *y1 + cub(x1) *y4)))/((x1 - x2) *(x1 - x3) *(x2 - 
        x3) *(x1 - x4) *(x2 - x4) *(x3 - x4));

  return d;

}


double CubicFrom4points(double* Y, double* X, double x)
{
  double x1 = X[0];
  double x2 = X[1];
  double x3 = X[2];
  double x4 = X[3];
  double y1 = Y[0];
  double y2 = Y[1];
  double y3 = Y[2];
  double y4 = Y[3];

  double a = (x1* (x1 - x4) * x4 * (y2 - y3) + 
      sqr(x3)* (x4 *y1 + x1* y2 - x4* y2 - x1* y4) + 
      sqr(x2)* (-x4* y1 - x1* y3 + x4* y3 + x3 *(y1 - y4) + x1* y4) + 
      x2* (sqr(x4) *(y1 - y3) + sqr(x1) * (y3 - y4) + sqr(x3) *(-y1 + y4)) + 
      x3* (sqr(x4) *(-y1 + y2) + sqr(x1) *(-y2 + y4)))/((x1 - x2) *(x1 - 
        x3) *(x2 - x3) *(x1 - x4) *(x2 - x4) *(x3 - x4)); 
 
  double b = (-x1 *(x1 - x4) *x4 *(x1 + x4) *(y2 - y3) + 
      x3 *(cub(x4) *(y1 - y2) + cub(x1) *(y2 - y4)) + 
      cub(x3) *(-x4* y1 - x1 *y2 + x4 *y2 + x1 *y4) + 
      cub(x2) *(x4 *y1 + x1 *y3 - x4* y3 - x1 *y4 + x3 *(-y1 + y4)) + 
      x2 *(cub(x4) *(-y1 + y3) + cub(x3) *(y1 - y4) + cub(x1) *(-y3 + y4)))/((x1 - 
        x2) *(x1 - x3) *(x2 - x3) *(x1 - x4) *(x2 - x4) *(x3 - x4));

  double c = (sqr(x1) *(x1 - x4) * sqr(x4) *(y2 - y3) + 
      cub(x3) *(sqr(x4) *(y1 - y2) + sqr(x1) *(y2 - y4)) + 
      sqr(x2) *(cub(x4) *(y1 - y3) + cub(x1) *(y3 - y4) + cub(x3) *(-y1 + y4)) + 
      sqr(x3) *(cub(x4) *(-y1 + y2) + cub(x1) *(-y2 + y4)) + 
      cub(x2) *(sqr(x4) *(-y1 + y3) + sqr(x3) *(y1 - y4) + 
         sqr(x1) *(-y3 + y4)))/((x1 - x2) *(x1 - x3) *(x2 - x3) *(x1 - 
        x4) *(x2 - x4) *(x3 - x4));

  double d = (x1 *(x1 - x3) *x3 *(x1 - x4) *(x3 - x4) *x4 *y2 + 
      x2 *sqr(x4) *(-cub(x3) *y1 + sqr(x3) *x4 *y1 + sqr(x1) *(x1 - x4) *y3) + 
      sqr(x1) *x2 *sqr(x3) *(-x1 + x3) *y4 + 
      cub(x2) *(x4 *(-sqr(x3) *y1 + x3 *x4 *y1 + x1 *(x1 - x4) *y3) + 
         x1 *x3 *(-x1 + x3) *y4) + 
      sqr(x2) *(x1* x4 *(-sqr(x1) + sqr(x4)) *y3 + cub(x3) *(x4 *y1 - x1* y4) + 
         x3 *(-cub(x4) *y1 + cub(x1) *y4)))/((x1 - x2) *(x1 - x3) *(x2 - 
        x3) *(x1 - x4) *(x2 - x4) *(x3 - x4));

  return a*cub(x) + b*sqr(x) + c*x + d;




/*
{a -> (x1 (x1 - x4) x4 (y2 - y3) + 
      x3^2 (x4 y1 + x1 y2 - x4 y2 - x1 y4) + 
      x2^2 (-x4 y1 - x1 y3 + x4 y3 + x3 (y1 - y4) + x1 y4) + 
      x2 (x4^2 (y1 - y3) + x1^2 (y3 - y4) + x3^2 (-y1 + y4)) + 
      x3 (x4^2 (-y1 + y2) + x1^2 (-y2 + y4)))/((x1 - x2) (x1 - 
        x3) (x2 - x3) (x1 - x4) (x2 - x4) (x3 - x4)), 
  b -> (-x1 (x1 - x4) x4 (x1 + x4) (y2 - y3) + 
      x3 (x4^3 (y1 - y2) + x1^3 (y2 - y4)) + 
      x3^3 (-x4 y1 - x1 y2 + x4 y2 + x1 y4) + 
      x2^3 (x4 y1 + x1 y3 - x4 y3 - x1 y4 + x3 (-y1 + y4)) + 
      x2 (x4^3 (-y1 + y3) + x3^3 (y1 - y4) + x1^3 (-y3 + y4)))/((x1 - 
        x2) (x1 - x3) (x2 - x3) (x1 - x4) (x2 - x4) (x3 - x4)), 
  c -> (x1^2 (x1 - x4) x4^2 (y2 - y3) + 
      x3^3 (x4^2 (y1 - y2) + x1^2 (y2 - y4)) + 
      x2^2 (x4^3 (y1 - y3) + x1^3 (y3 - y4) + x3^3 (-y1 + y4)) + 
      x3^2 (x4^3 (-y1 + y2) + x1^3 (-y2 + y4)) + 
      x2^3 (x4^2 (-y1 + y3) + x3^2 (y1 - y4) + 
         x1^2 (-y3 + y4)))/((x1 - x2) (x1 - x3) (x2 - x3) (x1 - 
        x4) (x2 - x4) (x3 - x4)), 
  d -> (x1 (x1 - x3) x3 (x1 - x4) (x3 - x4) x4 y2 + 
      x2 x4^2 (-x3^3 y1 + x3^2 x4 y1 + x1^2 (x1 - x4) y3) + 
      x1^2 x2 x3^2 (-x1 + x3) y4 + 
      x2^3 (x4 (-x3^2 y1 + x3 x4 y1 + x1 (x1 - x4) y3) + 
         x1 x3 (-x1 + x3) y4) + 
      x2^2 (x1 x4 (-x1^2 + x4^2) y3 + x3^3 (x4 y1 - x1 y4) + 
         x3 (-x4^3 y1 + x1^3 y4)))/((x1 - x2) (x1 - x3) (x2 - 
        x3) (x1 - x4) (x2 - x4) (x3 - x4))}}
*/

}


complex<double> CubicFrom4points(complex<double>* Y, double* X, double x)
{
  double ReY[4];
  double ImY[4];  
  for(int i = 0; i<4; i++)
  {
    ReY[i] = real(Y[i]);
    ImY[i] = imag(Y[i]);
  }

  return complex<double>(CubicFrom4points(ReY, X, x), CubicFrom4points(ImY, X, x));
}

complex<double> CubicFrom4points(complex<double>* Y, double* X)
{
  double ReY[4];
  double ImY[4];  
  for(int i = 0; i<4; i++)
  {
    ReY[i] = real(Y[i]);
    ImY[i] = imag(Y[i]);
  }

  return complex<double>(CubicFrom4points(ReY, X), CubicFrom4points(ImY, X));
}



//-----------------------------------------------------------------------------------------------------//




//--------------------------------------------- IO -----------------------------------------------------//


void PrintFunc3D(const char* FileName, int N, complex<double>** Y, double* X)          
{ 
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i+=20)
    for (int j=0; j<N; j+=20)
      fprintf(f,"%.15le %.15le %.15le %.15le\n", X[i], X[j], real(Y[i][j]), imag(Y[i][j]));
  fclose(f);
}

void PrintFunc(const char* FileName, int N, int M, double** Y, double* X)
{ 
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)   
  { 
    fprintf(f, "%.15le", X[i]);
    for (int j=0; j<M; j++)
    { 
       // loop through and store the numbers into the file
       fprintf(f, "%.15le", Y[i][j] );
    }
   fprintf(f, "\n");
  }
  fclose(f);
}

void PrintFunc(const char* FileName, int N, complex<double>* Y, double* X)          
{ 
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%.15le %.15le %.15le\n", X[i], real(Y[i]), imag(Y[i]));
  fclose(f);
}

void PrintFunc(const char* FileName, int N, double* Y, double* X)
{ 
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)  
    fprintf(f,"%.15le %.15le\n", X[i], Y[i]);
  fclose(f);
}

void PrintFunc(const char* FileName, int N, double* Y)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%d %.15le \n", i, Y[i]);
  fclose(f);
}

void PrintFunc(const char* FileName, std::vector< complex<double> > Y, std::vector<double> X)          
{ 
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<X.size(); i++)
    fprintf(f,"%.15le %.15le %.15le\n", X[i], real(Y[i]), imag(Y[i]));
  fclose(f);
}

void GetDimensions(const char* FileName, int &N, int &M)
{
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");
  
  //----count rows----//
  int i=0;
  char str[1000];  
  while (!feof(f))
  {  
    fgets ( str, 1000, f ); 
    if (str[0]=='\n') continue;
    //printf(str);

    i++;
  }
  N=i-1;
  fclose(f);
  //----count columns---//
  f = fopen(FileName, "r");
  i=1;
  int j=0;
  fgets ( str, 1000, f ); 
  while (str[i+1] != '\n')
  {  if ((str[i]!=' ')and(str[i-1]==' ')) j++;   
     i++;
  }
  M=j+1;

  //---close file-----//
  fclose(f);
}

void ReadFunc(const char* FileName, int &N, int &M, double** &X)
{ 
  GetDimensions(FileName, N, M);
  printf("N: %d M: %d\n", N, M);
 
  X = new double*[N];
  for (int i=0; i<N; i++)
    X[i] = new double[M];

  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) printf("FILE NOT FOUND!");
  
  for (int i=0; i<N; i++)
    for (int j=0; j<M; j++)   		
    { double Dummy;
      fscanf(f, "%le", &Dummy);
      X[i][j]=Dummy;
    }
  fclose(f);
}

void ReadFunc(const char* FileName, int &N, double* &Y, double* &X)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");
  
  //----count rows----//
  int i=0;
  char str[1000];  
  int prelines = 0;
  while (!feof(f))
  {  
    fgets ( str, 1000, f ); 
    if (str[0]=='\n') continue;
    if (str[0]=='#') { prelines++; continue; };
    i++;
  }
  N=i-1;
  printf("N: %d, prelines: %d \n", N, prelines);
  fclose(f);
 
  X = new double[N];
  Y = new double[N];

  f = fopen(FileName, "r");

  for (int i=0; i<prelines; i++)
    fgets ( str, 1000, f ); 
   
  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2,Dummy3;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);
    X[i]=Dummy1;
    Y[i]=Dummy2;
  }
  fclose(f);
}


void ReadFunc(const char* FileName, int &N, complex<double>* &Y, double* &X, bool PurelyReal=false)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");
  
  //----count rows----//
  int i=0;
  char str[1000];  
  int prelines = 0;
  while (!feof(f))
  {  
    fgets ( str, 1000, f ); 
    if (str[0]=='\n') continue;
    if (str[0]=='#') { prelines++; continue; };
    i++;
  }
  N=i-1;
  printf("N: %d, prelines: %d \n", N, prelines);
  fclose(f);
 
  X = new double[N];
  Y = new complex<double>[N];

  f = fopen(FileName, "r");

  for (int i=0; i<prelines; i++)
    fgets ( str, 1000, f ); 
   
  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2,Dummy3;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);
    if (not PurelyReal) fscanf(f, "%le", &Dummy3);
    X[i]=Dummy1;
    Y[i]=complex<double>(Dummy2, (PurelyReal) ? 0.0 : 0.0/*Dummy3*/);
  }
  fclose(f);
}
//---------------- vectors and matrices--------------------//

void Copy( int N, double** A, double** B, int offset )
{
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    B[i+offset][j+offset] = A[i+offset][j+offset];
}

void Multiply(int N, double** A,  int A_offset_i,  int A_offset_j,  
                     double** B,  int B_offset_i,  int B_offset_j,
                     double** AB, int AB_offset_i, int AB_offset_j)
{ 
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
  { AB[ AB_offset_i + i ][ AB_offset_j + j ] = 0.0;
    for(int k=0; k<N; k++)
      AB[ AB_offset_i + i ][ AB_offset_j + j ] += A[ A_offset_i + i ][ A_offset_j + k ] * B[ B_offset_i + k ][ B_offset_j + j ];
  }
}

void Multiply(int N, double** A, double** B, double** AB, int offset, int offset_j)
{
  int offset_i; 
  if (offset_j<0) offset_j = offset;
  if (offset<0) { offset_i = 0; offset_j=0; }
  else offset_i = offset;
  
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
  { AB[ offset_i + i ][ offset_j + j ] = 0.0;
    for(int k=0; k<N; k++)
      AB[ offset_i + i ][ offset_j + j ] += A[ offset_i + i ][ offset_j + k ] * B[ offset_i + k ][ offset_j + j ];
  }
}

void MultiplyByMatrix(int N, double* v, double** m)
{
  double* res = new double[N];
  for (int j=0; j<N; j++)
  { res[j] = 0;
    for (int i=0; i<N; i++)
      res[j] += m[i][j]*v[i]; 
  }
  for (int i=0; i<N; i++) v[i] = res[i];
  delete [] res;
}

void CreateRotationMatrix(double** m, int N, double angle, int* plane)
{
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) 
      if (i==j) m[i][j] = 1.0;
      else m[i][j] = 0.0;
  m[plane[0]][plane[0]] = cos(angle);
  m[plane[1]][plane[1]] = cos(angle);
  m[plane[0]][plane[1]] = sin(angle);
  m[plane[1]][plane[0]] = -sin(angle);
}

void RotateVector(int N, double* v, double angle, int* plane)
{ 
  double** RotationMatrix = new double*[N];  
  for (int i=0; i<N; i++)
    RotationMatrix[i] = new double[N];
  CreateRotationMatrix(RotationMatrix, N, angle, plane);
  MultiplyByMatrix(N, v, RotationMatrix);
  delete [] RotationMatrix;
}

void RotateVector2D(double* v, double angle)
{
  RotateVector(2, v, angle, (int []){0,1} );
}


#define NRANSI
#define TINY 1.0e-20

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=new double[n+1];
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) printf("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
        delete [] vv;
}

void InvertMatrix(int N, double** A, double** invA, double &det, int offset, int offset_j)
{
  int offset_i; 
  if (offset_j<0) offset_j = offset;
  if (offset<0) { offset_i = 0; offset_j=0; }
  else offset_i = offset;

  double** a = new double*[N+1];
  for(int i=0; i<=N; i++) 
  { a[i] = new double[N+1];
    for(int j=0; j<=N; j++) 
      a[i][j] = ( (i>0) and (j>0) ) ? A[i-1+offset_i][j-1+offset_j] : 0.0; 
  }

  //PrintMatrixOnScreen(N, A);
  //PrintMatrixOnScreen(N+1, a);


  double *col = new double[N+1];
  double d;
  int i,j,*indx = new int[N+1];

  ludcmp(a,N,indx,&d); //Decompose the matrix just once.
  //printf("ludcmp done\n");
  //PrintMatrixOnScreen(N+1, a);
  det = d;
  for(j=1;j<=N;j++) det *= a[j][j];

  
  for(j=1;j<=N;j++) 
  { //Find inverse by columns.
    for(i=1;i<=N;i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(a,N,indx,col);
    for(i=1;i<=N;i++) invA[i-1+offset_i][j-1+offset_j]=col[i];
  }

  for(j=0;j<=N;j++)
    delete [] a[j];
  delete [] a;
  
  delete [] indx;
  delete [] col;
}

double Determinant(int N, double** A, int offset, int offset_j)
{ 
  int offset_i; 
  if (offset_j<0) offset_j = offset;
  if (offset<0) { offset_i = 0; offset_j=0; }
  else offset_i = offset;

  double** a = new double*[N+1];
  for(int i=0; i<=N; i++) 
  { a[i] = new double[N+1];
    for(int j=0; j<=N; j++) 
      a[i][j] = ( (i>0) and (j>0) ) ? A[i-1+offset_i][j-1+offset_j] : 0.0; 
  }

  //PrintMatrixOnScreen(N, A);
  //PrintMatrixOnScreen(N+1, a);


  double *col = new double[N+1];
  double d;
  int i,j,*indx = new int[N+1];

  ludcmp(a,N,indx,&d); //Decompose the matrix just once.
  //printf("ludcmp done\n");
  //PrintMatrixOnScreen(N+1, a);
  double det = d;
  for(j=1;j<=N;j++) det *= a[j][j];
  
  delete [] indx;
  delete [] col;
  for(j=0;j<=N;j++)
    delete [] a[j];
  delete [] a;

  return det;
}

#undef TINY
#undef NRANSI



void PrintMatrixForGnuplot(const char* FN, int N, double** A)
{
  FILE* f = fopen(FN, "w");
  for(int i=0; i<N; i++)
  { for(int j=0; j<N; j++) 
      fprintf(f,"%d %d %.15le\n", i, j, A[i][j]);
    fprintf(f,"\n");
  }
  fclose(f);
}


void PrintMatrixInFile(const char* FN, const char* tag, int N, double** A)
{
  FILE* f = fopen(FN, tag);
  for(int j=0; j<N; j++)
  { for(int i=0; i<N; i++) 
       if (A[i][j]==0.0)
         fprintf(f,".\t\t");
       else
         fprintf(f,"%.3f\t\t", A[i][j]);
     fprintf(f,"\n");
  }
  fprintf(f,"\n");
  fclose(f);
}

void PrintMatrixOnScreen(int N, double** A)
{
  for(int i=0; i<N; i++) 
  { for(int j=0; j<N; j++) 
       if (A[i][j]==0.0)
         printf(".\t");
       else
         printf("%4.3f\t", A[i][j]);
     printf("\n");
  }
  printf("\n");
}

void PrintMatrixOnScreen(int N, double** A, int offset_i, int offset_j)
{
  for(int i=0; i<N; i++) 
  { for(int j=0; j<N; j++) 
       if (A[i+offset_i][j+offset_j]==0.0)
         printf(".\t");
       else
         printf("%4.3f\t", A[i+offset_i][j+offset_j]);
     printf("\n");
  }
  printf("\n");
}

void PrintMatrixOnScreen(int N, int M, double** A)
{
  for(int j=0; j<M; j++)
  { for(int i=0; i<N; i++) 
       if (A[i][j]==0.0)
         printf(".\t");
       else
         printf("%.3f\t", A[i][j]);
     printf("\n");
  }
  printf("\n");
}

void PrintMatrixOnScreen(int N, int M, int** A)
{
  for(int j=0; j<M; j++)
  { for(int i=0; i<N; i++) 
       if (A[i][j]==0)
         printf(".\t");
       else
         printf("%d\t", A[i][j]);
     printf("\n");
  }
  printf("\n");
}


                              
//------------------ integral routine ---------------------//
double TrapezIntegral(int N, double Y[], double X[])
{
  
  double sum = Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*(X[i+1]-X[i-1]);
  return sum*0.5;
}

complex<double> TrapezIntegral(int N, complex<double> Y[], double X[])
{
  complex<double> sum = Y[0]*complex<double>(X[1]-X[0]) +
                         Y[N-1]*complex<double>(X[N-1]-X[N-2]);
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*complex<double>(X[i+1]-X[i-1]);

  return sum*0.5;
}

//------------------ integral routine ---------------------//

double TrapezIntegralMP(int N, double Y[], double X[])
{
#ifdef _OMP

  double sum = 0.0;
  double sum0 = Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);
  int Nt;
  double* psum = new double[8];
  #pragma omp parallel shared(psum)
  { Nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    psum[tid] = 0;
    for (int i=tid+1; i<N-1; i+=Nt)
      psum[tid]+=Y[i]*(X[i+1]-X[i-1]);
    //printf("inside proc %d, psum = %le\n",tid,psum[tid]);
  }
  for (int i = 0; i<Nt; i++)
  {  sum += psum[i];
     //printf("outside proc %d, psum = %le\n",i,psum[i]);
  }
  delete [] psum;
  return (sum+sum0)*0.5;

#else

  return TrapezIntegral(N,Y,X);

#endif
}

complex<double> TrapezIntegralMP(int N, complex<double> Y[], double X[])
{

#ifdef _OMP
  complex<double> sum0 = Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);

  complex<double> sum = 0.0;
  int Nt;
  complex<double>* psum = new complex<double>[8];
  #pragma omp parallel shared(psum)
  { Nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    psum[tid] = 0;
    for (int i=tid+1; i<N-1; i+=Nt)
      psum[tid]+=Y[i]*(X[i+1]-X[i-1]);
  }
  for (int i = 0; i<Nt; i++) sum += psum[i];
  delete [] psum;
  return (sum + sum0)*0.5;

#else

  return TrapezIntegral(N,Y,X);

#endif

}

/*
double TrapezIntegral(std::vector<double> Y, std::vector<double> X)
{
  double sum = 0.0;
  int N = Y.size();
  sum += Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*(X[i+1]-X[i-1]);
  return sum*0.5;
}
*/
complex<double> TrapezIntegral(std::vector< complex<double> > Y, std::vector<double> X)
{
  complex<double> sum = 0.0;
  int N = Y.size();
  sum += Y[0]*complex<double>(X[1]-X[0]) + 
         Y[N-1]*complex<double>(X[N-1]-X[N-2]);
  
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*complex<double>(X[i+1]-X[i-1]);
  
  return sum*complex<double>(0.5);
}



double EllipticIntegralFirstKind(double x)
{ //integration goes from 0 to 1
  //split integration interval
  int Nlog=1000;
  int Nlin=1000;
  int N;
  double* omega;
  GRID grid;
  grid.InitGrid(Nlog,Nlin,1.0,1.0e-1,1.0e-6);
  grid.GetGrid(N,omega);

  double* integrand = new double[N];  
  for (int i=N/2; i<N; i++) integrand[i] = 0;
  for (int i=0; i<N/2; i++) integrand[i] = 1.0 / sqrt( (1-sqr(omega[i]+1.0))*(1-sqr((omega[i]+1.0)*x))  );
  
  double res = TrapezIntegral(N,integrand,omega);
  delete [] integrand;
  return res;
}

complex<double> EllipticIntegralFirstKind(complex<double> x)
{ //integration goes from 0 to 1
  //split integration interval
  int Nlog=3000;
  int Nlin=3000;
  int N;
  double* omega;
  GRID grid;
  grid.InitGrid(Nlog,Nlin,1.0,1.0e-1,1.0e-14);
  grid.GetGrid(N,omega);

  complex<double>* integrand = new complex<double>[N];  
  for (int i=N/2; i<N; i++) integrand[i] = 0;
  for (int i=0; i<N/2; i++) integrand[i] = 1.0 / (sqrt( complex<double> ( (1.0 - sqr(omega[i]+1.0))
                                                                              * ( 1.0 - sqr((omega[i]+1.0)*x) )
                                                                         )
                                                       )
                                                  );
  /*char FN[50];
  sprintf(FN,"elliptic.x%.3f",real(x));
  PrintFunc(FN,N, integrand, omega);*/
  complex<double> res = TrapezIntegral(N,integrand,omega);
  delete [] integrand;
  delete [] omega;
  grid.~GRID();
  return res;
}

//------Sine Integral, from Num Recipes -----//

double SI(double x)
{ 
  //--------consts----------//
  double EPS = 3.0e-16;
  double EULER = 0.577215664901533;
  int MAXIT = 100;
  double PIBY2 = 1.570796326794897;
  double FPMIN = 1.0e-30;
  double TMIN = 2.0;
  complex<double> ONE =  complex<double>(1.0,0.0);
  //--------vars-------------//
  double si;
  int i,k,odd;
  double a,err,fact,sign,sum,sumc,sums,t,term;
  complex<double> h,b,c,d,del;
  t=fabs(x);
                                           // Special case.
  if (t == 0.0) {
      si=0.0;
//      *ci = -1.0/FPMIN;
//      return si;
  }
                                            //Evaluate continued fraction by modified
  if (t > TMIN) {
                                               // Lentz’s method (§5.2).
      b=complex<double>(1.0,t);
      c=complex<double>(1.0/FPMIN,0.0);
      d=h=ONE/b;
      for (i=2;i<=MAXIT;i++) {
          a = -(i-1)*(i-1);
          b+=complex<double>(2.0,0.0);
                                            //Denominators cannot be zero.
          d=ONE/( complex<double>(a,0)*d +b );
          c=b+complex<double>(a,0.0)/c;
          del=c*d;
          h=h*del;
          if (fabs(real(del)-1.0)+fabs(imag(del)) < EPS) break;
      }
      if (i > MAXIT) printf("cf failed in cisi\n");
      h=complex<double>(cos(t),-sin(t))*h;
//      *ci = -real(h);
      si=PIBY2+imag(h);
                                            //Evaluate both series simultaneously.
  } else {
                                            //Special case: avoid failure of convergence
      if (t < sqrt(FPMIN)) {
                                               // test because of underflow.
          sumc=0.0;
          sums=t;
      } else {
          sum=sums=sumc=0.0;
          sign=fact=1.0;
          odd=1;
          for (k=1;k<=MAXIT;k++) {
              fact *= t/k;
              term=fact/k;
              sum += sign*term;
              err=term/fabs(sum);
              if (odd) {
                  sign = -sign;
                  sums=sum;
                  sum=sumc;
              } else {
                  sumc=sum;
                  sum=sums;
              }
              if (err < EPS) break;
              odd=!odd;
          }
          if (k > MAXIT) printf("maxits exceeded in cisi\n");
      }
      si=sums;
//      *ci=sumc+log(t)+EULER;
  }
  if (x < 0.0) si = -(si);
  return si;
}


//-------------------- DOSes and Init Deltas -----------------//

double DOS(int DOStype, double t, double om)
{ 
  switch (DOStype)
  {
    case DOStypes::SemiCircle: 
            if ( sqr(2*t)<sqr(om) )
              return 0.0;
            else 
              return sqrt( sqr(2*t)-sqr(om) ) 
                     / ( 2*pi*sqr(t) ); 
            break;
    case DOStypes::Gaussian: 
            { 
              double d = 1/(2*t*sqrt(pi))*exp( -sqr( om )/( 4*sqr(t) ) ); 
              return (d>0.001) ? d : 0;
            }
            break;
    case DOStypes::Insulator:
           {
            double a=-2.5, b =2.5;
            if ( sqr(2*t)<sqr(om-a) and sqr(2*t)<sqr(om-b))
                return 0.0;
            else 
              if (sqr(2*t)>sqr(om-a))
                return sqrt( sqr(2*t)-sqr(om-a) )
                       / ( 4*pi*sqr(t) );
              else
                return sqrt( sqr(2*t)-sqr(om-b) )
                       / ( 4*pi*sqr(t) );                                                                                
           } break;     
    case DOStypes::SquareLattice:
           { if (om==0.0) return 100.00;
             if (abs(om)>=4*t) return 0.0;
             else return 1.0 / (sqr(pi)*2*t) 
                         * EllipticIntegralFirstKind( sqrt( 1.0-sqr(om/(4*t))) );
           }  break;  
    case DOStypes::CubicLattice:           
           {  if (abs(om)>=6*t) return 0.0;
              int N = 6000;
              double* omega;
              GRID grid;
              grid.InitGrid(0, N, 3.141, 0.0, 0.0);
              grid.GetGrid(N,omega);
              double* d = new double[N];
              for (int i=0; i<N/2; i++) d[i]=0;
              for (int i=N/2; i<N; i++)
              {  complex<double> f = 4.0 * t / ( om - 2 * t * cos(omega[i]) + complex<double>(0.0, 0.001) );
                 //if (f>1.0)
                   d[i] = imag( f * EllipticIntegralFirstKind( f ) );
                 //else d[i]=0;
                 //if (om<1.0) printf("omega %.3f: f=%.3f\n",om,f);
              }
              /*if (om<1.0)
              {  char FN[50];
                 sprintf(FN,"integrand.w%.3f",om);              
                 PrintFunc(FN,N,d,omega);
              }*/

              double Res = -1.0/pi * TrapezIntegral(N,d,omega);
              delete [] d;
              delete [] omega;
              grid.~GRID();
              return 1.0 / (2 * sqr(pi) * t) 
                     * Res;          
           }  break;
    //Add here different DOS types!!! 
    default: { printf("DOStype %d not implemented!", DOStype); exit(1); }
  }
}

void WriteCubicDosToFile()
{ 
  int N = 3000;
  double* omega;
  double* dos = new double[N];
  GRID grid;
  grid.InitGrid(0, N, 3.5, 0.0, 0.0);
  grid.GetGrid(N,omega);
  for(int i=N/2; i<N; i++)
  { dos[i] = DOS(DOStypes::CubicLattice, 0.5, omega[i]);
    dos[N-1-i]=dos[i];  
    printf("Write2File: omega %.3f: %.3f\n", omega[i], dos[i]); 
  }
  PrintFunc("CubicLatticeDOS",N,dos,omega);
  delete [] dos;
  delete [] omega;
}

void ReadDosFromFile(const char* FN, int N, double* omega, double* DOS)
{
  int n,m;
  double** dos;
  ReadFunc(FN, n, m, dos);
  for (int i = 0; i<N; i++)
    DOS[i] = interpl(n, dos[1], dos[0], omega[i]);
  for (int i=0; i<m; i++)
    delete [] dos[i];
  delete [] dos;
}

double interpl(int N, double* Y, double* X, double x)
{
  if ((x < X[0])or(x > X[N-1])) return 0;
  else
  { int i;
    for (i = 0; i < N; i++)
      if (X[i]>x) break;
    return Y[i-1] + (Y[i]-Y[i-1]) / (X[i]-X[i-1])
                                  * (x-X[i-1]);
  }
}

complex<double> interpl(int N, complex<double>* Y, double* X, double x)
{
  double* ReY = new double[N];
  double* ImY = new double[N];
  for(int i=0; i<N; i++)
  {
    ReY[i] = real(Y[i]);
    ImY[i] = imag(Y[i]);
  }

  double ReYx = interpl(N, ReY, X, x);
  double ImYx = interpl(N, ImY, X, x);

  delete [] ReY;
  delete [] ImY;

  return complex<double>(ReYx,ImYx);
}

void InitDOS(int DOStype, double t, int N, double* omega, double* dos)
{
  for (int i=0; i<N; i++) dos[i] = DOS(DOStype, t, omega[i]);
}

void InitDelta(int DOStype, 
               int N, 
               double V, 
               double mu, 
               double eta, 
               double t,
               double* omega,  
               complex<double>* Delta,
               const char* FNDos)
{  
  printf("Making  Delta: t=%.3f\n",t);
  // prepare dos
  
//------------------Get DOS----------------------//
  double* dos = new double[N];
  if (FNDos=="")
    for(int i=0; i<N; i++)
      dos[i] = DOS(DOStype, t, omega[i]);
  else
    ReadDosFromFile(FNDos, N, omega, dos);

  PrintFunc("InitDelta.DOS",N,dos,omega);
  
//-----------------------------------------------//
  for (int i=0; i<N; i++)
  {  
    //treat integrand carefully
    double D = 0.0;
    complex<double> LogTerm = 0.0;
    D = dos[i];

    LogTerm =  ((i==0)or(i==N-1)) ? 0.0 
                                  : complex<double>(D, 0.0) * log( complex<double> ( (mu + omega[i] + omega[N-1])
                                                                                    /(mu + omega[i] - omega[N-1]) ) );
    //printf("LogTerm(%.3f) = %.3f, %.3f\n",omega[i], real(LogTerm),imag(LogTerm));
    //create integrand array
    complex<double>* d = new complex<double>[N];
    for (int j=0; j<N; j++)
    {  d[j] = (i==j) ? 0.0 : complex<double>(dos[j] - D, 0.0) 
                            / complex<double>( mu + omega[i] - omega[j] ); 
       //if (i==j) printf("d(i=j) = %.3f, %.3f\n", real(d[j]), imag(d[j]) );
    }
      
    //integrate to get G 
    Delta[i] = - sqr(V)*(TrapezIntegral(N, d,omega) + LogTerm) ; 
    delete [] d; 
  }
  PrintFunc("DeltaMade",N,Delta,omega);

//---------------------------------------------------------------------------//

  delete [] dos;
  printf("Intgral Delta: %.6f", imag(TrapezIntegral(N,Delta,omega)));
}

void InitDeltaFromSIAMOutputFile(const char* FN, int N, double* omega, complex<double>* Delta)
{
  int n,m;
  double** output;
  ReadFunc(FN, n, m, output);
  for (int i = 0; i<N; i++)
    Delta[i] = complex<double>( interpl(n, output[2], output[0], omega[i]), 
                              interpl(n, output[3], output[0], omega[i]) );
  for (int i=0; i<m; i++)
    delete [] output[i];
  delete [] output;
}

void InitInsulatorDelta(double U,
                        int N, 
                        double t,  
                        double* omega,
                        complex<double>* Delta)
{ printf("Making insulator Delta from Hubbar I approx: U=%.3f, t=%.3f\n",U,t);
  for (int i=0; i<N; i++)
  { 
    double Gat = 0.5 * ( 1.0/(omega[i] + U/2.0) + 1.0/(omega[i] - U/2.0) );
    complex<double> Sqr = sqrt( complex<double>(1.0)/sqr(Gat)  - 4.0*sqr(t) );
    complex<double> G = (1.0/Gat + real(Sqr) < 0) ? (1.0/Gat + conj(Sqr) ) / (2*sqr(t))
                                                  : (1.0/Gat - Sqr ) / (2*sqr(t));
    Delta[i] = sqr(t)*G;
  }
  PrintFunc("initdelta",N,Delta,omega);
   printf("Intgral Delta: %.6f", imag(TrapezIntegral(N,Delta,omega)));
}

void InitArrays(int N, double*** a, const int* Ns)
{
  for (int i=0; i<N; i++)
    a[i][0] = new double[Ns[i]];  
}

void InitArrays(int N, double*** a, const int* Ns, double (**func)(int))
{
  for (int i=0; i<N; i++)
  {   a[i][0] = new double[Ns[i]];  
      if (func[i]!=NULL)
        for (int j=0; j<Ns[i]; j++) 
        {  a[i][0][j]=func[i](j);
           printf("j: %d a: %f\n",j,a[i][0][j]);
        }
  }
}
      
void ReleaseMemory(int N, double*** a)
{
  for (int i=0; i<N; i++)
    delete [] a[i][0];
     
}





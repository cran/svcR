# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <ctype.h>
# include "R.h"
# include "Rinternals.h"
# include "Rdefines.h"

FILE * f = 0;


int  fibiter(int a, int b, int c);
void fib(int*n);
void arsim(double *x, long *n,double *phi);

void SizeMat_C(	char  ** NomFile,
				double * ListMax);
void LireMat_C(	char  ** NomFile, 
				int    * MaxLin, 
				int    * MaxCol, 
				double * OutMat);
void MinMaxMat_C (	double	* matrix, 
					int		* NRows , 
					int		* NCols , 
					int		* Cx, 
					int		* Cy , 
					double	* MinMaxXY);

int	Compare( const void *a, const void *b);
int cmpint (const void * it1, const void * it2);

void     calcKernelMatrix_C (   double	* matrix, 
								int		* MatDist, 
								int		* q, 
								int		* Ncols, 
								int		* NRows, 
								int		* KernChoice, 
								double	* MatKernel);
double   Kernel_C ( int		 q, 
					double * Vec1, 
					double * Vec2, 
					int		 Ncols, 
					int		 Choice );
double   KernelLinear_C (   double	* Vec1, 
							double	* Vec2, 
							int		  Ncols );
double   KernelGaussian_C ( int		 q, 
							double * Vec1, 
							double * Vec2, 
							int		 Ncols );
double   KernelGaussianDist_C ( int		 q, 
								double * Vec1 );

void    CalcWcluster_C (	double * MatKern, 
							int    * MxIter, 
							int    * nlin, 
							double * MaxValA, 
							double * niu, 
							int    * MnW, 
							int    * MaxW, 
							double * AroundNull, 
							double * VectorsYA );
double  CritereWcluster ( double * VecteurA, 
						  double * MatrixK, 
						  int      nlin );
int     ConstraintCluster1 ( double   BoundS, 
							 double * VecteurA, 
							 int      nlin );
int     ConstraintCluster2 (   double * VecteurA, 
							   int      nlin, 
							   double   AroundNull );
double * MakeA ( int    nlin, 
				 double niu, 
				 double AroundNull, 
				 double MaxValA );

void     ListSVPoints_C ( double	* VectorsYA , 
						int		* nlin , 
						double	* MxValA, 
						double	* ArouNullVA , 
						double	* ListPoints);
void     vectorWcluster (	double	* DataMatrix, 
							int		* nlin, 
							int		* ncol , 
							double	* VectorsYA, 
							double	* VecW );

void     scalarRO ( double	* DataMatrix, 
					int		* q, 
					int		* nlin, 
					int		* ncol, 
					double	* VecW, 
					int		* Kchoice, 
					double	* VectorsYA, 
					double	* ro );
void     RadiusCluster ( double	* VectorsYA, 
						 int	* nlin, 
						 double	* MxValA, 
						 double	* ArouNullVA, 
						 double	* MatrixK , 
						 double	* R);
void     RadiusPoint (  double	* Vec, 
						int		* q, 
						int		* nlin, 
						int		* ncol , 
						double	* VectorsYA, 
						int		* Kchoice, 
						double	* DataMatrix, 
						double	* MatrixKern, 
						double	* Rad );
void     RadiusData ( int		* IndicePoint, 
					  int		* nlin , 
					  double	* VectorsYA, 
					  double	* MatrixKern, 
					  double	* Rad );
void     SmallR  (  int		* nlin, 
					double	* RadiusC, 
					double	* niu, 
					double	* VectorsYA , 
					double	* MatrixKern, 
					double	* resu );

void    ClusterLabeling_C (	double	* DataMatrix, 
							double	* MatrixKern,
							int		* Ngrid, 
							int		* ncol, 
							int		* nlin, 
							int		* q, 
							double	* MxX,
							double	* MnX, 
							double	* MxY, 
							double	* MnY,
							double	* RC, 
							double	* smallr,
							int		* KChoice, 
							double	* VectorsYA,
							double	* NumPoints );
void    NbCluster_C (	double	* NumPoints, 
						int		* Ngrid, 
						int		* NCluster ) ;
void    MatchGridPoint_C (	double	* DataMatrix,
							double	* MatrixKern,
							int		* Ngrid,
							int		* ncol,
							int		* nlin,
							double	* MxX,
							double	* MnX,
							double	* MxY,
							double	* MnY,
							int		* NumCluster,
							int		* Knn,
							int		* Cx,
							int		* Cy,
							double	* NumPoints,
							double	* ClassPoints );
void    Evaluation_C (	double	* DataMatrix, 
						int		* nlin, 
						int		* ncol, 
						int		* NCluster,
						double		* ClassPoints, 
						double	* ListMis );

# include"svcR.h"

int fibiter(int a, int b, int c){
     if( c <= 0) return b;
     else return fibiter (a+b,a,c-1);
}
void fib(int*n){
     *n=fibiter(1,0,*n);
}

//---------------------------------------------------------------------------

int Compare( const void *a, const void *b)
{
   return ( *(int *) a - *(int *)b );
}

int cmpint (const void * it1, const void * it2)
{
   int i1 = *(int *) it1;
   int i2 = *(int *) it2;
   if (i1 > i2) return 1;
   else if (i1 < i2) return -1;
   else return 0;
}

//---------------------------------------------------------------------------

void SizeMat_C( char **NomFile,	double   * ListMax ) {

FILE    *       I           = NULL;
char            Chaine[500] = "aa";
char            pszVal[256] = "2";
int             ComptVal    = 0;
int             Lin = 0, Col = 0, NbCol = 3, NbLin = 0;
double  *       InMat = NULL;
int             MaxLin, MaxCol;
int             MinL = 1000000, MinC = 1000000;

MaxLin = 0; MaxCol = 0;

if( (I= fopen(NomFile[0], "r")) == NULL ){
		fclose(I);
        exit(1);
} //finif

fseek(I, 0, SEEK_SET);
while( fgets(Chaine, 500, I) != NULL ) {
        ComptVal++;
} //finwhile

NbLin = ComptVal;
InMat = calloc( (NbCol+1) * (NbLin+1) , sizeof (double) );

ComptVal = 0;
fseek(I, 0, SEEK_SET);
while( fgets(Chaine, 500, I) != NULL ) {
        tolower(Chaine);
        sscanf(Chaine,"%d %d %s",&Lin,&Col,pszVal);
        *(InMat + ((ComptVal * NbCol) + 0)) = Lin ;
        *(InMat + ((ComptVal * NbCol) + 1)) = Col ;
        *(InMat + ((ComptVal * NbCol) + 2)) = atof(pszVal) ;
        ( MaxLin < Lin ) ? (MaxLin=Lin):(MaxLin=MaxLin);
        ( MaxCol < Col ) ? (MaxCol=Col):(MaxCol=MaxCol);
        (  MinL  > Lin ) ? (MinL=Lin):(MinL=MinL);
        (  MinC  > Col  ) ? (MinC=Col):(MinC=MinC);
        ComptVal++ ;
} //while

MaxCol = MaxCol - MinC + 1;
MaxLin = MaxLin - MinL + 1;

*(ListMax)   = *&MaxLin;
*(ListMax+1) = *&MaxCol;

free(InMat);
fclose(I);
I = 0;

return ;
}

//---------------------------------------------------------------------------

void LireMat_C(char **NomFile, int *MaxLin, int *MaxCol, double *OutMat)
{
FILE    *       I           = NULL;
char            Chaine[500] = "aa";
char            pszVal[256] = "2";
int             ComptVal    = 0;
int             Lin = 0, Col = 0, NbCol = 3, NbLin = 0;
int             iCtr, jCtr;
double  *       InMat = NULL;

if( (I= fopen(NomFile[0], "r")) == NULL ){
		fclose(I);
        exit(1);
}

fseek(I, 0, SEEK_SET);
while( fgets(Chaine, 500, I) != NULL ) {
        ComptVal++;
} //finwhile

NbLin = ComptVal;
InMat = calloc( NbCol * NbLin , sizeof (double) );

ComptVal = 0;
fseek(I, 0, SEEK_SET);
while( fgets(Chaine, 500, I) != NULL ) {
        tolower(Chaine);
        sscanf(Chaine,"%d %d %s",&Lin,&Col,pszVal);
        *(InMat + ((ComptVal * NbCol) + 0)) = Lin ;
        *(InMat + ((ComptVal * NbCol) + 1)) = Col ;
        *(InMat + ((ComptVal * NbCol) + 2)) = atof(pszVal) ;
        ComptVal++ ;
} //while


for (iCtr = 0; iCtr < NbLin; iCtr++) {
        for (jCtr = 0; jCtr < NbCol; jCtr++) {
                Lin = *(InMat + ((iCtr * NbCol) + 0)) -1;
                Col = *(InMat + ((iCtr * NbCol) + 1)) -1;
                *((OutMat) + ((Lin * *MaxCol) + Col)) = *(InMat + ((iCtr * NbCol) + 2));
        }
}

free(InMat);
fclose(I);
I = 0;

return ;
}

//########################################################################
//# Compute min max of matrix
//#
//########################################################################

//# Usage:
//#   MinMaxMat(Mat=matrice, Cx=1, Cy=2);
//#

void MinMaxMat_C ( double * matrix, int * NRows , int * NCols, int * Cx, int * Cy , double * MinMaxXY) {

int             i = 0;
int             N = 0;
int             M = 0;
double *        MaxX = 0;
double *        MinX = 0;
double *        MaxY = 0;
double *        MinY = 0;
double          Cell = 0;
int             DimX = 0;
int             DimY = 0;

N    = * NRows;
M    = * NCols;
MaxX = (MinMaxXY);
MinX = (MinMaxXY+1);
MaxY = (MinMaxXY+2);
MinY = (MinMaxXY+3);
DimX = * Cx -1;
DimY = * Cy -1;

*MaxX = -1000000;
*MinX = +1000000;
*MaxY = -1000000;
*MinY = +1000000;

for ( i = 0; i < N; i++) {

        Cell = *(matrix + ((i * M) + DimX));
        (*MaxX < Cell)? (*MaxX=Cell) : (*MaxX=*MaxX);
        (*MinX > Cell)? (*MinX=Cell) : (*MinX=*MinX);

        Cell = *(matrix + ((i * M) + DimY));
        (*MaxY < Cell)? (*MaxY=Cell) : (*MaxY=*MaxY);
        (*MinY > Cell)? (*MinY=Cell) : (*MinY=*MinY);

} //endfori

//fprintf(f, "MinMax %6.2lf \t %6.2lf \t %6.2lf \t %6.2lf \n",*MaxX,*MinX,*MaxY,*MinY  );


return;

}

//########################################################################
//# Calc Kernel Matrix
//#
//########################################################################

//# Usage:
//#   calcKernelMatrix(matrix=M)
//#

void calcKernelMatrix_C ( double * matrix, int *MatDist , int * q, int * Ncols, int * NRows, int * KernChoice, double * MatKernel) {

int			ncol	= *Ncols;
int			nrow	= *NRows;
int			i	= 0;
int			j	= 0;
double             *    Vec1    = 0;
double             *    Vec2    = 0;
double             *    Cell    = 0;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

for(i = 0; i < nrow; i++) {

	j = i ;

	while( j < nrow ) {

		*((MatKernel) + ((i * nrow) + j)) = 0;

		Vec1 = ( matrix + i*ncol );
                Vec2 = ( matrix + j*ncol );
                Cell = ( matrix + i*ncol + j );

                if( *MatDist == 0 ) {
                        *((MatKernel) + ((i * nrow) + j)) = Kernel_C ( *q, Vec1, Vec2, ncol, *KernChoice ) ;
                } else {
                        *((MatKernel) + ((i * nrow) + j)) = Kernel_C ( *q, Cell, Cell, ncol, *KernChoice ) ;
                } //finif

		*((MatKernel) + ((j * nrow) + i)) = *((MatKernel) + ((i * nrow) + j));

		j=j+1;

	} //endwhile

} //endfori

     /*   for ( i = 0; i < nrow; i++) {
                for ( j = 0; j < ncol; j++) {
                        fprintf(f, "%6.2lf ",*(matrix + ((i * ncol) + j)) );
                }
                fprintf(f, "\n");
        }
        fprintf(f, "\n");
        for ( i = 0; i < nrow; i++) {
                for ( j = 0; j < nrow; j++) {
                        fprintf(f, "%1.20lf ",*(MatKernel + ((i * nrow) + j)) );
                }
                fprintf(f, "\n");
        }
*/
//fclose(f); f=0;

return ;
}

//########################################################################
//# Define Kernel
//#
//########################################################################

//# Usage:
//#   Kernel(Vec1=v1, Vec2=v2, Choice=0)
//#

double Kernel_C ( int q, double * Vec1, double * Vec2, int Ncols, int Choice ) {

double Res = 0;

if( Choice == 0 )
		Res = KernelLinear_C(Vec1, Vec2, Ncols);
else if( Choice == 1)
		Res = KernelGaussian_C(q, Vec1, Vec2, Ncols);
else if( Choice == 2)
		Res = KernelGaussianDist_C(q, Vec1);

return (Res);
}

//########################################################################
//# Define KernelLinear
//#
//########################################################################

//# Usage:
//#   KernelLinear(V1=v1, V2=v2)
//#

double KernelLinear_C (double * Vec1, double * Vec2, int Ncols ) {

double	Res     = 0;
double	ValVec1 = 0;
double	ValVec2 = 0;
int		j       = 0;

for(j = 0; j < Ncols-1; j++) {

		ValVec1 = *(Vec1  + j);
		ValVec2 = *(Vec2  + j);
		Res = Res + ( ValVec1 * ValVec2 );
		//fprintf(f, "%d \t %d \t %d\n",ValVec1, j, Res );
}//endfori

return(Res);
}

//########################################################################
//# Define KernelGaussian
//#
//########################################################################

//# Usage:
//#   KernelGaussian(V1=v1, V2=v2)
//#

double KernelGaussian_C (int q, double * Vec1, double * Vec2, int Ncols ) {

double	Res     = 0;
double	ValVec1 = 0;
double	ValVec2 = 0;
int	j       = 0;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

for(j = 0; j < Ncols-1; j++) {

		ValVec1 = *(Vec1  + j);
		ValVec2 = *(Vec2  + j);
		Res = Res + (( ValVec1 - ValVec2 ) * ( ValVec1 - ValVec2 ));
}//endfori

//fprintf(f, "res %1.2lf \n", Res );
Res = exp( - (q * Res) );

//fprintf(f, "%d \t %1.2lf \t %1.2lf \t %d \t %1.20lf\n", q, ValVec1, ValVec2, j, Res );

// fclose(f);
// f=0;

return(Res);
}

//########################################################################
//# Define KernelGaussian with element computing a distance
//#
//########################################################################

//# Usage:
//#   KernelGaussianDist(V1=v1)
//#

double KernelGaussianDist_C (int q, double * Vec1 ) {

double  Res     = 0;
double  ValVec1 = 0;

ValVec1 = *(Vec1) ;

Res = exp( - (q * sqrt( ValVec1*ValVec1 )) );
//fprintf(f, "%1.2lf \t %1.2lf \n",ValVec1, Res );

return(Res);
}

//########################################################################
//# Define Criterium for Clustering
//#
//########################################################################

//# Usage:
//#   CritereWcluster(VecteurA=Va, MatrixK=M)
//#

double CritereWcluster (double * VecteurA , double * MatrixK, int nlin ) {

int    i     = 0;
int    j     = 0;
int    N     = 0;
double W     = 0;
double SomIJ = 0;

// W est le critere a maximiser
N = nlin;

for( i = 0 ; i < N ; i++ ) {
        SomIJ = SomIJ + *(VecteurA+i) * *(MatrixK + N*i + i);
	for( j = 0 ; j < N ; j++ ) {
                SomIJ = SomIJ - *(VecteurA+i) * *(VecteurA+j) * *(MatrixK + N*i + j);
	} //fin for j
	//fprintf(f, "SomIJ %1.6lf %1.6lf \n",*(VecteurA+i),*(MatrixK + N*i + j) );
} //fin for i

W = SomIJ ;

return (W);
}


//########################################################################
//# Define Constraint Cluster 1
//#
//########################################################################

//# Usage:
//#   ConstraintCluster1(BoundS=bs, VecteurA=Va)
//#

int ConstraintCluster1 (double BoundS, double *VecteurA, int nlin ) {

int N = 0;
int i = 0;

N = nlin;

for( i = 0 ; i < N ; i++ ) {
	if( *(VecteurA + i) < 0 || *(VecteurA + i)  > BoundS ){
		return( 0 );
	} //endif
} //fin for i

return( 1 );

}

//########################################################################
//# Define Constraint Cluster 2
//#
//########################################################################

//# Usage:
//#   ConstraintCluster2(VecteurA=Va)
//#

int ConstraintCluster2 ( double * VecteurA, int nlin, double AroundNull ) {

int    i   = 0;
double Sum = 0;
int    N   = 0;

N = nlin;

for( i = 0 ; i < N ; i++ ) {
        Sum = Sum + *(VecteurA + i);
} //fin for i

if( (Sum >  (1 + AroundNull)) || (Sum < (1 - AroundNull)) )  //en pratique
	return( 0 );

return( 1 );

}

//########################################################################
//# Define Make Vector A
//#
//########################################################################

//# Usage:
//#   MakeA(Dim=N)
//#

double * MakeA (int nlin, double niu, double AroundNull, double MaxValA ) {

double   nu     = 0;
int      i      = 0;
int      isum   = 0;
int      N      = 0;
int      K      = 0;
double * COEFFS = 0;
double   Z      = 0;
double   Sum    = 0;
time_t   t;

nu     = niu;
N      = nlin;
COEFFS = calloc( 2*N+1 , sizeof(double) );

srand((unsigned) time(&t));
for( i = 0 ; i < N ; i++ )
        *(COEFFS+i) = (rand()%100)/100.0;
Z       = (rand()%100)/110.0;

for( i = 0 ; i < N ; i++ ){
	K = 1;
	if( *(COEFFS+i) <= Z )
		*(COEFFS+i)     = 0;
	if( *(COEFFS+i) > 0.9 && K <= nu*N) {
		*(COEFFS+i)	= 1/(N*nu);
		K               = K+1;
	}
	else if( *(COEFFS+i) > Z )
		*(COEFFS+i) = *(COEFFS+i)/(nu*N);
} //#finfor

for( isum = 0 ; isum < N ; isum++ )
        Sum = Sum + *(COEFFS+isum);
        
i = 1;
while( Sum > 1 && i <= N ) {
	*(COEFFS+i)	= 0;
	i		= i+1;

        Sum = 0;
        for( isum = 0 ; isum < N ; isum++ )
                Sum = Sum + *(COEFFS+isum);
} //finwhile

while( Sum < (1 - AroundNull) && i <= N ) {
	*(COEFFS+i)	= MaxValA;
	i		= i+1;

        Sum = 0;
        for( isum = 0 ; isum < N ; isum++ )
                Sum = Sum + *(COEFFS+isum);
} //finwhile

return(COEFFS);
}

//########################################################################
//# Define Alpha factors
//#
//########################################################################

//# Usage:
//#   CalcWcluster_C(V1=v1)
//#

void CalcWcluster_C (double * MatKern, int * MxIter, int *nlin, double* MaxValA, double* niu, int* MnW, int* MaxW, double *AroundNull, double * VectorsYA ) {

double  W         = *MaxW;
double  ValW      = 0 ;
double  TabWPrec  = 0 ;
double  PreviousW = *MaxW ;
double  MinW      = *MnW ;
double  BoundSup  = 1E+10;
double  nu        = *niu;
double  MxValA	  = *MaxValA;
double  AroundNul = *AroundNull;
int     N	      = *nlin;
int     Iter      = 1;
int     MaxIter   = *MxIter;
int     i         = 0;
double * VWA      = 0;

//WYA	<- list(W="",Y="",A="")

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "w");

//creation of a random vector (randomization uniforme runif) / (randomization gaussienne rnorm)
VWA       = (double*)MakeA( N , nu, AroundNul, MxValA );         // WYA$A

// computation of WYA$A bound sup
if( nu && N )
	BoundSup = MxValA ; // 1 / ( nu * N );

TabWPrec = MinW;

while( (W > TabWPrec &&  Iter <= MaxIter) ||  (Iter <= MaxIter) ) {
	
	PreviousW = W;
	W         = CritereWcluster(VWA, MatKern, N);

	if( Iter > 1 ){
		TabWPrec = ValW;
		}

	if( W > TabWPrec ){
		ValW = W;
	}
	else {
		ValW = TabWPrec;
		free(VWA); 
		VWA =  MakeA( 2*N+1 , nu, AroundNul, MxValA);
		while(
			(ConstraintCluster1(BoundSup, VWA, N) == 0)
			||
			(ConstraintCluster2(VWA, N, AroundNul) == 0)
		) {
								free(VWA); 
                                VWA =  MakeA( N , nu, AroundNul, MxValA);
		} //#fin while
	} //#finif

	PreviousW = W;
	W         = CritereWcluster(VWA, MatKern, N);
	Iter      = Iter + 1;
} //finwhile

for ( i = 0; i < 2*N+1; i++){
        *(VectorsYA+i) = *(VWA+i);
        if( i == N )
                *(VectorsYA+N) = W;
        //fprintf(f, "%1.6lf \n",*(VectorsYA+i) );
} //endfor

//fprintf(f, "\n");
//fprintf(f, "%2.20lf \t %2.20lf \n",*(VectorsYA+N), W );

//fclose(f);
//f = 0;

return;
}

//########################################################################
//# Define Support Vector List
//#
//########################################################################

//# Usage:
//#   ListSVPoints(VA=va)
//#

void ListSVPoints_C ( double * VectorsYA , int * nlin , double * MxValA, double * ArouNullVA , double * ListPoints) {

int             N               = 0;
int    *        ListP           = 0;
double          MaxValA         = 0;
double          AroundNullVA    = 0;
int             Signe           = 0;
int             i               = 0;
int             ValVA           = 0;
int             IndVal          = 0;

N		= * nlin;
AroundNullVA    = * ArouNullVA;
MaxValA         = * MxValA;
ListP           = calloc( N , sizeof(int) );

Signe = 0;
for( i = 0; i < N ; i++ ) {
	if( *(VectorsYA+i) > AroundNullVA && *(VectorsYA+i) < MaxValA ){
	        *(ListP+i) = i;	Signe = 1;
	} //finif
} //#fin for i

if( Signe == 0){
	ValVA = 1; IndVal = 1;
	for( i = 0; i < N; i++ ) {
		if( *(VectorsYA+i) < ValVA ){
			ValVA  = *(VectorsYA+i);
			IndVal = i;
		} //#endif
	} //#fin for i
	*(ListP+IndVal) = IndVal;
} //#endif

for ( i = 0; i < N; i++){
        *(ListPoints+i) = *(ListP+i);
       // fprintf(f, "%d \n",*(ListPoints+i) );
} //endfor

return;
}

//########################################################################
//# Calcule Vecteur W for SVC
//#
//########################################################################

//# Usage:
//#   vectorWcluster(Mat=Matrice,WYA=wya);
//#

void vectorWcluster ( double * DataMatrix, int * nlin, int * ncol , double * VectorsYA, double * VecW ) {

double    *        iVecW           = 0;
int                IndLig          = 0;
int                IndCol          = 0;
int                i               = 0;
int                N               = 0;
int                M               = 0;

N		= * nlin;
M		= * ncol;
iVecW           = calloc( N , sizeof(double) );

for(IndCol = 0; IndCol < M-1; IndCol++ ) {
        for(IndLig = 0; IndLig < N; IndLig++ ) {
	        *(iVecW+IndCol) = *(iVecW+IndCol) + *(VectorsYA+IndLig) * *(DataMatrix + IndLig*M + IndCol);
        }  //#fin for IndLig
} //#fin for IndCol

/*for ( i = 0; i < M-1; i++){
        *(VecW+i) = *(iVecW+i);
        fprintf(f, "%1.6lf \n",*(VecW+i) );
} //endfor
*/

return ;
}

//########################################################################
//# Calcule Scalar RO for SVC
//#
//########################################################################

//# Usage:
//#   scalarRO(Mat=Matrice,VecW=w,WYA=wya);
//#

void scalarRO ( double * DataMatrix, int * q, int * nlin, int * ncol, double * VecW, int * Kchoice, double * VectorsYA, double * ro ) {

int      KC   = 0;
int      N    = 0;
int      M    = 0;
int      i    = 0;
double * Vec1 = 0;
double * Vec2 = 0;

KC = * Kchoice;
N  = * nlin;
M  = * ncol;

for( i = 0; i < N; i++ ) {

        Vec1 = (VecW);
        Vec2 = (DataMatrix + i*M);
	*ro  = *ro + *(VectorsYA+i) * Kernel_C( *q, Vec1, Vec2, *ncol, KC);

} //#fin for i

//fprintf(f, "\n%1.6lf \n", *ro );

return ;
}

//########################################################################
//# Define Radius of the cluster
//#
//########################################################################

//# Usage:
//#   RadiusCluster(VA=va, MatrixK=Mat)
//#

void RadiusCluster ( double * VectorsYA, int * nlin, double * MxValA, double * ArouNullVA, double * MatrixK , double * R) {

int             IndiceSV        = 0 ;
int             N               = 0;
int             i               = 1;
double          MaxValA         = 0;
double          AroundNullVA    = 0;
double          Min             = 0;

MaxValA         = *MxValA;
AroundNullVA    = *ArouNullVA;
N               = * nlin;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

while( IndiceSV == 0 && i <= N ) {
        if( *(VectorsYA+i) > AroundNullVA && *(VectorsYA+i) < MaxValA )
		IndiceSV = i ;
	i = i + 1;
} //#fin while i

Min = 1000;
if( IndiceSV == 0) for(i = 0; i < N; i++ ) {
        if( *(VectorsYA+i) < Min ){
		IndiceSV = i ;
		Min = *(VectorsYA+i);
	} //#endif
} //#finfori

//fprintf(f, "\n%d \n", IndiceSV );

if( IndiceSV != 0 ) {
        RadiusData(&IndiceSV, &N, VectorsYA, MatrixK, R);
}
else {	
        *R = 0;
} //endif

//fprintf(f, "\n R %1.6lf \n", *R );

//fclose(f); f = 0;

return ;
}

//########################################################################
//# Define Radius of a vector
//#
//########################################################################

//# Usage:
//#   RadiusPoint(Vec=v, VA=wva, Mat=Mat, MatK=MatK)
//#

void RadiusPoint ( double * Vec, int * q, int * nlin, int * ncol , double * VectorsYA, int * Kchoice, double * DataMatrix, double * MatrixKern, double *Rad ) {

int             M       = 0;
int             KC      = 0;
int             N       = 0;
int             i       = 0;
int             j       = 0;
double *        Vec1    = 0;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

N               = * nlin;
M               = * ncol;
KC              = * Kchoice;

*Rad = 0;

for( i = 0; i < N; i++  ) {

	for( j = 0; j < N; j++ ) {

		 *Rad = *Rad + ( *(VectorsYA+i) * ( *(VectorsYA+j) )  * ( *(MatrixKern + i*N + j) ) );

	} //#fin for j
			//fprintf(f, "\nin %1.6lf \n", *Rad );

	Vec1 = (DataMatrix+i*M);

	*Rad = *Rad - 2 * *(VectorsYA+i) * Kernel_C( *q, Vec1, Vec, M , KC );

} //#fin for i

*Rad = *Rad + Kernel_C( *q, Vec, Vec, M , KC );

//fprintf(f, "\n%d \n", *q );

if( *Rad > 0 ){
	*Rad = sqrt( *Rad );
}
else {
	*Rad = 0;
} //#finif

//fclose(f);
//f = 0;

return;
}

//########################################################################
//# Define Radius of a point
//#
//########################################################################

//# Usage:
//#   RadiusData(IndicePoint=1, VA=wva, MatK=Mat)
//#

void RadiusData (int * IndicePoint, int * nlin , double * VectorsYA, double * MatrixKern, double *Rad ) {

int             N       = 0;
int             i       = 0;
int             j       = 0;

*Rad    = 0 ;
N       = * nlin;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

for( i = 0; i < N; i++ ) {

	for( j = 0; j < N; j++ ) {
	        *Rad = *Rad + *(VectorsYA+i) * *(VectorsYA+j) * *(MatrixKern + i*N + j);
	} //#fin for j

        *Rad = *Rad - 2 * *(VectorsYA+i) * *(MatrixKern + i*N + *IndicePoint);

		//fprintf(f, "\n rad %1.20lf vec %1.20lf mat %1.20lf\n", *Rad , *(VectorsYA+i), *(MatrixKern + i*N + *IndicePoint) );

} //#fin for i

*Rad = *Rad + *(MatrixKern + *IndicePoint*N + *IndicePoint) ;


//fprintf(f, "\n r1 %1.6lf \n", *Rad );

if( *Rad > 0 ){
	*Rad = sqrt( *Rad );
}
else {
	*Rad = 0;
} //#finif

//fprintf(f, "\n r2 %1.6lf \n", *Rad );

//fclose(f); f = 0;

return ;
}

//########################################################################
//# Define small variation of r around R
//#
//########################################################################

//# Usage:
//#   SmallR(VA=wva, MatK=Mat)
//#

void SmallR  (  int * nlin, double * RadiusC, double * niu, double * VectorsYA , double * MatrixKern, double *resu ) {

double  C       = 0;
int     N       = 0;
double  nu      = 0;
int     i       = 0;
double  R       = 0;
double *r       = 0;
double  rmean   = 0;
double  rmax    = 0;
double  rmin    = 0;
int     compt   = 0;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

nu      = * niu;
N       = * nlin;
C       = 1 / ( nu * N );

r = calloc( N , sizeof(double) );

for( i = 0; i < N; i++ ) {
        
	R      = 0;
	*(r+i) = 0;

	if( *(VectorsYA+i) > C * 0.99 ) {
		RadiusData(&i, &N, VectorsYA, MatrixKern, &R);
		if( R > *RadiusC )
		        *(r+i) =  R*R - *RadiusC * *RadiusC;
	} //#finif
	
} //#finfor

//# calculation of the mean and the max
for( i = 0; i < N; i++ ) {

	if( *(r+i) != 0 ){
		rmean = rmean + *(r+i) ;
		compt = compt + 1;
	} //#finif
	if( fabs( *(r+i) ) > rmax )
		rmax = fabs( *(r+i) );

	if( fabs( *(r+i) ) < rmin )
		rmin = fabs( *(r+i) );

} //#finfor
if( compt != 0) 
	rmean = rmean / compt;

//*resu = abs(rmean);
*resu = rmax;

//fprintf(f, "\nresu %1.6lf %1.6lf \n", *resu,rmax  );
//fclose(f);f=0;

return  ;
}

//########################################################################
//# Compute labeling clusters
//#
//########################################################################

//# Usage:
//#   ClusterLabeling(Mat=matrice, MatK= matriceK, Cx=1,Cy=2, WYA=va);
//#

void ClusterLabeling_C ( double * DataMatrix, double * MatrixKern,
                int * Ngrid, int * ncol, int * nlin, int * q, double * MxX,
                double * MnX, double * MxY, double * MnY,
                double * RC, double * smallr,
                int * KChoice, double * VectorsYA,
                double * NumPoints ) {

int             i               = 0;
int             j               = 0;
int             N               = 0;
int             Col             = 0;
int             Lin             = 0;
int             InP             = 0;
int             KC              = 0;
double          MaxX	        = 0;
double          MinX	        = 0;
double          MaxY	        = 0;
double          MinY	        = 0;
double *        ListVec         = 0;
short  *        PointGrid       = 0;
double          x               = 0;
double          y               = 0;
double          RadiusC         = 0;
double          Rad             = 0;
double          r               = 0;
double          RR	            = 0;
int             IndLC	        = 0;
int             Signal	        = 0;
int             ii              = 0;
int             jj              = 0;
int    *        ListClusters    = 0;
int             IndListCluster  = 0;
int             NbElemCluster   = 0;
int             NbElemClusterCur= 0;
int             IndMembListClust= 0;
int             ActualPoint     = 0;
int             num             = 0;
int             k               = 0;
int             l               = 0;
int             NumberCluster   = 0;
short   *       ClassConnex     = 0;
int             IndClassFusIn   = 0;
int             IndClassFusOut  = 0;
int             indColMatConnex = 0;
int             indListClusters = 0;
int             IndElem         = 0;
int             ComptRecii      = 0;
int             ComptRecjj      = 0;
int             NumClus         = 0;

InP     = 1;
N       = *Ngrid;
MaxX	= *MxX;
MinX	= *MnX;
MaxY	= *MxY;
MinY	= *MnY;
Lin     = *nlin;
Col     = 3;
ListVec =  calloc( Col , sizeof(double) );
KC      = *KChoice;
RadiusC = *RC;
r       = *smallr;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

PointGrid    = calloc( N*N , sizeof(short) );
RR			 = (RadiusC*RadiusC + r);

for(i = 0; i < N; i++ ){
        for(j = 0; j < N; j++ ){

                x  = ( MaxX - MinX )*( (i) / (double)(N-1) ) + MinX ; ListVec[0] = x;
                y  = ( MaxY - MinY )*( (j) / (double)(N-1) ) + MinY ; ListVec[1] = y;  ListVec[2] = 1;

                RadiusPoint ( ListVec, q, &Lin, &Col , VectorsYA, &KC, DataMatrix, MatrixKern, &Rad );

                *(PointGrid +i*N +j) = (short) InP;

                *(NumPoints +InP*6 + 1) = (double) InP;
                *(NumPoints +InP*6 + 2) = (double) i;
                *(NumPoints +InP*6 + 3) = (double) j;
                *(NumPoints +InP*6 + 5) = (double) 0;

		if(  RR >=  Rad*Rad ) {
                        *(NumPoints +InP*6 + 4) = (double)1;
		}
		else {
			*(NumPoints +InP*6 + 4) = (double)0;
		} //#finif

		InP = InP +1;

	} //#finforj

} //#finfori

//#sortie test
/*
fprintf(f, "Grid \n" );
for ( j = (N-1); j != -1; j-- ){
        for( i = 0; i < N ; i++ ) {
                InP = (int) *(PointGrid +i*N + j);
                //fprintf(f, " %2.0lf ", *(NumPoints +InP*6 + 4) );
	} //endforj
        //fprintf(f, "\n" );
}  //endfori
*/

InP             = 0;
IndLC		= 0;
Signal		= 0;

for(i = 0; i < N; i++ ){
	for(j = 0; j < N; j++ ){

                InP = (short) *(PointGrid +i*N +j);

                if( *(NumPoints +InP*6 + 4) == 1) {
                        Signal	= 0;
			if( IndLC == 0 ) {	//# create a new cluster

                                ListClusters    = (int*) calloc( 2  , sizeof(int*) );

                                ListClusters[0] = (int) calloc( 1  ,  sizeof(int) );
                                ListClusters[0] = (int) 1;

                                ListClusters[1] = (int) calloc( 2  , sizeof(int) );
                                *( (int*)ListClusters[1]+0) = (int) 1;
                                *( (int*)ListClusters[1]+1) = (int) InP;

                                *(NumPoints +InP*6 + 5) = (double) 1;

				IndLC   = IndLC + 1;
				Signal  = 1;
			}
			else {				       //# look for assign to an existing cluster
				ComptRecii        = 1;
                                for( ii = (i-1); ComptRecii <= 3; ii++ ){
                                        ComptRecjj        = 1;
					for( jj = (j-1); ComptRecjj <= 3; jj++ ){


                                                if( ii <  0 ) ii = 0;   if( jj <  0 ) jj = 0;
						if( ii >= N ) ii = N-1; if( jj >= N ) jj = N-1 ;

						if( ii == i && jj == j ){ }
						else {
							if( Signal == 0)
                                                        for( IndListCluster = 1; IndListCluster <= ListClusters[0] ; IndListCluster++ ) {

                                                                NbElemCluster = *( (int*)(ListClusters[IndListCluster])+0);

								if( Signal == 0) for( IndMembListClust = 1; IndMembListClust <= NbElemCluster ; IndMembListClust++) {

                                                                        ActualPoint = *( (int*)(ListClusters[IndListCluster])+ IndMembListClust);

									if( *(NumPoints +ActualPoint*6 + 2) == ii && *(NumPoints +ActualPoint*6 + 3) == jj ){

                                                                                *( (int*)(ListClusters[IndListCluster])+0)      =   *( (int*)(ListClusters[IndListCluster])+0) +1;
                                                                                IndElem =  *( (int*)(ListClusters[IndListCluster])+0);
                                                                                ListClusters[IndListCluster]                    =  (int) realloc( (int*)ListClusters[IndListCluster]  , (IndElem+2) * sizeof(int) );
                                                                                *( (int*)(ListClusters[IndListCluster])+IndElem)=  (int) InP;

                                                                                *(NumPoints +InP*6 + 5)         = (double) IndListCluster;
										Signal                          = 1;

									} //#endif

								} //#finforIndMembListClust

							} //#finforIndListCluster
						} //#finifiijj
                                        ComptRecjj++;
					} //#endforjj
                                ComptRecii++ ;
                                } //#endforii
			} //#endif

			if( Signal == 0 ){ //#on cree un nouveau cluster

                                ListClusters[0]                 = (int) ListClusters[0] + 1;
                                num                             = (int )ListClusters[0];
                                ListClusters                    = (int*) realloc( ListClusters , (num+1) * sizeof(int*) );
                                ListClusters[num]               = (int) calloc( 2  , sizeof(int) );

                                *( (int*)(ListClusters[num])+0) = (int) 1;
                                *( (int*)(ListClusters[num])+1) = (int) InP;
                                *(NumPoints +InP*6 + 5)         = (double) num;

			} //#endif

		} //#endif

	} //#finforj
} //#finfori

//#sortie test
/*fprintf(f, "\n" );
for ( j = (N-1); j != -1; j-- ){
        for( i = 0; i < N ; i++ ) {

                InP = (short) *(PointGrid +i*N +j);
                fprintf(f, " %2.0lf ", *(NumPoints +InP*6 + 5) );
	} //endforj
        fprintf(f, "\n" );
}  //endfori
*/

if( ListClusters == 0 || ListClusters[0] == 0 ){
	free(PointGrid);	PointGrid=0;
	free(ListVec);		ListVec=0;
	free(ListClusters);	ListClusters=0;
	//fprintf(f, "fini\n" );
	//fclose(f);
	//f = 0;
	return;
}

//#calcul de la matrice de connexite entre classes
NumberCluster = (short) ListClusters[0];

ClassConnex = calloc( (NumberCluster+2)*(NumberCluster+2) , sizeof( short ) );

if( NumberCluster > 1 ) {

for(i = 0; i <= NumberCluster; i++ )
        for(j = 0; j <= NumberCluster; j++ )
	        *(ClassConnex +i*(NumberCluster+1) +j) = (short)0;

for(i = 0; i < N; i++ ) {
	for(j = 0; j < N; j++ ){

                InP = (short) *(PointGrid +i*N +j);

                ComptRecii      = 1;
                for( ii = (i-1); ComptRecii <= 3; ii++ ){
                        ComptRecjj        = 1;
                        for( jj = (j-1); ComptRecjj <= 3; jj++ ){

                                if( ii <  0  ) ii = 0;   if( jj <  0 ) jj = 0;
				if( ii >=  N ) ii = N-1; if( jj >= N ) jj = N-1 ;

				ActualPoint = (short) *(PointGrid +ii*N +jj);
				if( *(NumPoints +InP*6 + 5) != 0 &&  *(NumPoints +ActualPoint*6 + 5) != 0 ) {
					k = (int) *(NumPoints +InP*6 + 5);
					l = (int) *(NumPoints +ActualPoint*6 + 5);
					*(ClassConnex +k*(NumberCluster+1) +l) = (short) 1;
				} //#endif

                                ComptRecjj++;
			} //#endforjj
                ComptRecii++ ;
		} //#endforii
	} //#finforj
} //#finfori
       
} //#endifNumberCluster
 


//#deletion of bad classes
IndClassFusIn  = 0;
IndClassFusOut = NumberCluster;
while( IndClassFusOut > 1) {

	for( indColMatConnex = 1; indColMatConnex <=(IndClassFusOut-1) ; indColMatConnex++  )
        if( *(ClassConnex +indColMatConnex*(NumberCluster+1) +IndClassFusOut) == 1 ) {

                NbElemClusterCur                                = *( (int*)(ListClusters[indColMatConnex])+0);
             
                *( (int*)(ListClusters[indColMatConnex])+0)     = *( (int*)(ListClusters[indColMatConnex])+0)  + *( (int*)(ListClusters[IndClassFusOut])+0) ;
                NbElemCluster                                   = *( (int*)(ListClusters[indColMatConnex])+0) ;
                ListClusters[indColMatConnex]                   = (int) realloc( (int*)ListClusters[indColMatConnex]  , (NbElemCluster+2) * sizeof(int) );

                for( IndMembListClust = 1; IndMembListClust <= *( (int*)(ListClusters[IndClassFusOut])+0) ; IndMembListClust++)
                        *( (int*)(ListClusters[indColMatConnex]) + NbElemClusterCur + IndMembListClust) = *( (int*)(ListClusters[IndClassFusOut])+IndMembListClust);

                *( (short*)(ListClusters[IndClassFusOut])+0)    = (int) 0 ;
                ListClusters[0]                                 = (int) ListClusters[0] - 1;

	} //#finforindColMatConnex

	IndClassFusOut = IndClassFusOut - 1;

} //#finwhile

if( ListClusters[0] > 0 )
for( indListClusters = 1; indListClusters <= NumberCluster ; indListClusters++ ) {    //	# deletion of copied clusters (list)
	NbElemCluster   = *( (int*)(ListClusters[indListClusters])+0);
	if( NbElemCluster )
                NumClus++;
        for( IndMembListClust = 1; IndMembListClust <= NbElemCluster ; IndMembListClust++ ) {
		ActualPoint                     = *( (int*)(ListClusters[indListClusters])+IndMembListClust);
        *(NumPoints +ActualPoint*6 + 5) = (double) NumClus;        
	} //#finforIndMembListClust

} //#finforindListClusters

//#sortie test
/*
fprintf(f, "NumberCluster %d\n",NumClus );
for ( j = (N-1); j != -1; j-- ){
        for( i = 0; i < N ; i++ ) {

                InP = (int) *(PointGrid +i*N + j);
                fprintf(f, " %2.0lf ", *(NumPoints +InP*6 + 5) );

	} //endforj
        fprintf(f, "\n" );
}  //endfori
*/

//##### end of computation

free(PointGrid);PointGrid=0;
free(ClassConnex);  ClassConnex=0;
free(ListVec);   ListVec=0;
free(ListClusters); ListClusters=0;

/*fprintf(f, "fini\n" );
fclose(f);
f = 0;*/

return ;
}

//########################################################################
//# Compute Matching Points to Grid
//#
//########################################################################

//# Usage:
//#   MatchGridPoint(Mat=Mat, NumCluster=length(ListClusters), Cx= cx, Cy=cy, Knn=1);
//#

void NbCluster_C ( double * NumPoints, int * Ngrid, int * NCluster ) {

int             i               = 0;
int             j               = 0;
int             InP             = 0;
int             NG              = 0;
int    *        PointGrid       = 0;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

*NCluster       = 0;
NG              = *Ngrid;
PointGrid       = calloc( NG*NG , sizeof(int) );

for(i = 0; i < NG; i++ ){
        for(j = 0; j < NG; j++ ){
                *(PointGrid +i*NG + j) = InP;
		InP = InP +1;
        } //#endforj
} //#endforj

for ( i = 0; i < NG; i++ ){

        for(j = 0; j < NG; j++ ){

                InP = *(PointGrid +i*NG + j);
                ( *NCluster < *(NumPoints +InP*6 + 5) )? (*NCluster = *(NumPoints +InP*6 + 5)):(*NCluster=*NCluster);

        } //#endforj

}  //endfori

free(PointGrid); PointGrid=0;
//fclose(f);
//f = 0;


}

//########################################################################
//# Compute Matching Points to Grid
//#
//########################################################################

//# Usage:
//#   MatchGridPoint(Mat=Mat, NumCluster=length(ListClusters), Cx= cx, Cy=cy, Knn=1);
//#

void MatchGridPoint_C ( double * DataMatrix, double * MatrixKern,
                int * Ngrid, int * ncol, int * nlin, double * MxX,
                double * MnX, double * MxY, double * MnY, int *NCluster,
                int * Knn, int * Cx, int * Cy, double * NumPoints, double * ClassPoints ) {

int             i               = 0;
int             j               = 0;
int             NG              = 0;
int             Col             = 0;
int             Lin             = 0;
int             InP             = 0;
double          MaxX	        = 0;
double          MinX	        = 0;
double          MaxY	        = 0;
double          MinY	        = 0;
int             ii              = 0;
int             jj              = 0;
int    *        ScoreNeighbours = 0;
int             ActualPoint     = 0;
int             NumClass        = 0;
int             NPoint          = 0;
int             IndNumClass     = 0;
int				Xi              = 0;
int				Yi              = 0;
int    *        PointGrid       = 0;
int             ComptRecii      = 0;
int             ComptRecjj      = 0;
int             NumCluster      = 0;

InP     = 1;
NG      = *Ngrid;
MaxX	= *MxX;
MinX	= *MnX;
MaxY	= *MxY;
MinY	= *MnY;
Lin     = *nlin;
Col     = *ncol;
NPoint  = *nlin;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

PointGrid    = (int*) calloc( NG*NG , sizeof(int) );

for(i = 0; i < NG; i++ ){
        for(j = 0; j < NG; j++ ){
                *(PointGrid +i*NG + j) = InP;
		InP = InP +1;
        } //#endforj
} //#endforj

NumCluster = *NCluster;

if( NumCluster > 0)
for( i = 0; i < NPoint; i++ ){

	Xi  =  (int) ( fabs( *(DataMatrix + i*Col  + *Cx) - MinX )*(NG-1) / (double) ( MaxX - MinX )  ) ;
	Yi  =  (int) ( fabs( *(DataMatrix + i*Col  + *Cy) - MinY )*(NG-1) / (double) ( MaxY - MinY )  ) ;
	if( Xi >= NG) Xi = NG-1; if(Yi >= NG) Yi = NG-1;
	if( Xi <  0  ) Xi = 0;   if( Yi <  0 ) Yi = 0;

	ScoreNeighbours = (int*) calloc( (NumCluster+1) , sizeof(int) );

	ComptRecii        = 1;
	for( ii = (Xi-1); ComptRecii <= 3; ii++ ){

		ComptRecjj        = 1;  
		for( jj = (Yi-1); ComptRecjj <= 3; jj++ ){

			if( jj == Yi && ii == Xi){
				       //	#we leave
			}
			else {
				if( ii >= NG) ii = NG-1; if(jj >= NG) jj = NG-1;
				if( ii <  0)  ii = 0   ; if(jj < 0  ) jj = 0;

				ActualPoint = *(PointGrid +ii*NG + jj) ;

				if( *(NumPoints +ActualPoint*6 + 5) != 0 && *(NumPoints +ActualPoint*6 + 5) < (NumCluster+1) ){
					NumClass                    = *(NumPoints +ActualPoint*6 + 5);
					*(ScoreNeighbours+NumClass) = *(ScoreNeighbours+NumClass) + 1 ;
				} //#endif
			} //#endif
			ComptRecjj++;

		} //#endforjj
		ComptRecii++ ;
	} //#endforii

	for( IndNumClass = 1; IndNumClass <= NumCluster; IndNumClass++ ){
		if( *(ScoreNeighbours+ IndNumClass) >= *Knn )
			*(ClassPoints+i) = IndNumClass;
	} //#endforIndNumClass

	free(ScoreNeighbours);

} //#finfori

free(PointGrid); PointGrid=0;

/*fprintf(f, " matchgrid \n" );
for( i = 0; i < NPoint; i++ ){
fprintf(f, " %d \t %lf \t %lf \t %lf\n", i, (double) *(ClassPoints+i) , *(DataMatrix + i*Col  + *Cx), *(DataMatrix + i*Col  + *Cy));
} //#finfori
*/

//fclose(f);
//f = 0;

return ;  
}

//########################################################################
//# Compute Evaluation
//#
//########################################################################

//# Usage:
//#   Evaluation(Mat=matrice, NBClass=nb, Cx=1, Cy=2, ClassPoints=cl);
//#

void Evaluation_C ( double * DataMatrix, int * nlin, int * ncol, int *NCluster,
                double * ClassPoints, double * ListMis ) {

int             i                       = 0;
int             j                       = 0;
int             NumColClass             = 0;
int     *       ListSortedItemsByClass  = 0;
short   *       LabelClass              = 0;
short   *       CardLabelClass          = 0;
int             NBClass                 = 0;
int             NBCluster               = 0;
short   *       BestClass	        = 0;
int     **       ListVecClass	        = 0;
int             IndClass                = 0;
int             IndListMis	        = 0;
int             N                       = 0;
int             Counter                 = 0;
int             IndVecSorted            = 0;
double          Precision               = 0;
int             L                       = 0;
int             Elem                    = 0;
int             ElemPrec                = 0;
int             SizeVec                 = 0;

NBCluster       = *NCluster;
N               = *nlin;
NumColClass     = *ncol;

//f = fopen("D:\\R\\library\\svcR\\data\\sortieyy.txt", "a+");

//#output points class
//fprintf(f, " %d \n", NBCluster );
//fprintf(f, " %d \n", N );
//fprintf(f, " %d \n", NumColClass );

/*fprintf(f, " Evaluation_C ClassPoints \n" );
for( i = 0; i < N; i++ ){
fprintf(f, " %d \t %lf \n", i, (double) *(ClassPoints+i) );
} //#finfori
*/

if( N < 2 && NBCluster < 1 ){
	//fclose(f);
	//f = 0;
	return;
}

ListSortedItemsByClass  = (int*) calloc( N , sizeof(int) ); 

for( i = 0; i < N; i++  ){
        *(ListSortedItemsByClass+i) = (int) *(DataMatrix+i*NumColClass+ NumColClass-1);
} //endfori 

qsort( (int*)ListSortedItemsByClass , N , sizeof(int) , cmpint); 

//fclose(f); f = 0; return;
NBClass = ListSortedItemsByClass[N-1]; 

//#output points class

/*fprintf(f, " Evaluation_C ListSorted \n" );
for( i = 0; i < N; i++ ){
 fprintf(f, " %d \t %d \n", i, *(ListSortedItemsByClass+i) );
} //#finfori
*/

LabelClass              = (short*) calloc( NBClass*10 , sizeof(short) );
CardLabelClass          = (short*) calloc( NBClass*10 , sizeof(short) );

NBClass                 = 1;
LabelClass[NBClass]     = (short) ListSortedItemsByClass[0];
CardLabelClass[NBClass] = (short) 1;

for( i = 1; i < N; i++  ){
	if( ListSortedItemsByClass[i] == ListSortedItemsByClass[i-1] ){
		CardLabelClass[NBClass] = (short) CardLabelClass[NBClass] + 1;
	}
	else {
		NBClass		        = NBClass + 1;
		LabelClass[NBClass]     = (short) ListSortedItemsByClass[i];
        CardLabelClass[NBClass] = (short) 1;
	} //#endif

} //#finfori


//#output points class
/*fprintf(f, "\n\n Evaluation_C LabelClass \n" );
for( i = 1; i <= NBClass; i++ ){
 fprintf(f, " %d \t %d \t %d \n", i, LabelClass[i], CardLabelClass[i] );
} //#finfori
*/

BestClass          = (short*) calloc( 10*NBClass , sizeof(short) );
ListVecClass       = (int*)   calloc( 10*NBClass , sizeof(int) );


for( IndClass = 1; IndClass <= NBClass; IndClass++  ){

        ListVecClass[ IndClass ]                = (int*) calloc( 1 , sizeof(int) );
        *( (int*)ListVecClass[IndClass]+0)      = (int) 0;

} //endforIndClass

for( i = 0; i < N; i++ ){
        Elem = (int) *(DataMatrix +i*NumColClass +NumColClass-1);
	L = (int) *( (int*)(ListVecClass[Elem]) + 0 ) ;
        ListVecClass[ Elem ]  =  realloc( (int*)ListVecClass[ Elem ] , (L+2)*sizeof(int) );
        *( (int*)(ListVecClass[Elem]) + 0 )     = *( (int*)(ListVecClass[Elem]) + 0 ) + 1;
	*( (int*)(ListVecClass[Elem])+ (L+1))   = (int)ClassPoints[i];
} //#finFor

for( IndClass = 1; IndClass <= NBClass; IndClass++  ) {

        SizeVec = (int) *( (int*)(ListVecClass[ (short) LabelClass[IndClass] ])+ 0);
        qsort( (int*) (ListVecClass[(short)LabelClass[IndClass]]+1) , SizeVec , sizeof(int) , cmpint);

        //#output points class
      //  fprintf(f, " \n\nListVecClass \n" );
      //  for( i = 0; i < (int) *( (int*)(ListVecClass[ LabelClass[IndClass] ])+ 0); i++ ){
      //          Elem = (int) *( (int*)(ListVecClass[ LabelClass[IndClass] ])+ (i+1)) ;
      //          fprintf(f, " %d \t %d \n", i, Elem );
      //  } //#finfori
     
	Counter    = 1;

        *(BestClass+IndClass*2+ 2) = (short) 1;
        *(BestClass+IndClass*2+ 1) = (short) *( (int*)(ListVecClass[ (short)LabelClass[IndClass] ])+ 1);

	if( SizeVec >= 2 )
	for(  IndVecSorted = 2; IndVecSorted <= SizeVec; IndVecSorted++ ) {

                Elem     = *( (int*)(ListVecClass[ (short)LabelClass[IndClass] ])+ IndVecSorted);
                ElemPrec = *( (int*)(ListVecClass[ (short)LabelClass[IndClass] ])+ IndVecSorted - 1);

		if( Elem == ElemPrec ) {

			Counter = Counter + 1;
			if( Counter > *(BestClass+IndClass*2+ 2) ){

				*(BestClass+IndClass*2+ 2) = (short) Counter  ;
				*(BestClass+IndClass*2+ 1) = (short) ElemPrec ;

                        } //#finif
		}
		else {
                        Counter         = 1;
		} //#finif

	} //#finIndVecSorted
}//#finIndClass

//#output points class
//fprintf(f, "\n\n BestClass \n" );
//for( IndClass = 1; IndClass <= NBClass; IndClass++ ){
// fprintf(f, "IndClass %d \t %d \t %d \n", IndClass, *(BestClass+IndClass*2+ 2), *(BestClass+IndClass*2+ 1) );
//} //#finfori

for( IndClass = 1; IndClass <= NBClass; IndClass++  )
        if( *(BestClass+IndClass*2+ 1) != 0 )
	        Precision = Precision + *(BestClass+IndClass*2+ 2);
Precision = ( Precision/ (double)N ) * 100;

//#extraction of misclassified items
for( i = 0; i < N; i++ ) {
	for( IndClass = 1; IndClass <= NBClass; IndClass++ ) {
	        if( (int)*(DataMatrix+i*NumColClass+NumColClass-1) == (short) LabelClass[IndClass] ){

		        if( *(BestClass+IndClass*2+ 1) != ClassPoints[i] || *(BestClass+IndClass*2+ 1) == 0 ) {
					ListMis[IndListMis] = (double)i ;
					IndListMis          = IndListMis + 1;
			} //#endIf
                        
		} //#endIf
	} //#endforIndVec
} //#finFor

//#output points class

//fprintf(f, "\n\n Precision \t %lf \t IndListMis %d \n", Precision , IndListMis);
//for( i = 0; i < N; i++ ){
// fprintf(f, " %d \t %lf \n", i, ListMis[i] );
//} //#finfori


free(ListSortedItemsByClass);   ListSortedItemsByClass  = 0;
free(BestClass);                BestClass               = 0;
free(ListVecClass);             ListVecClass            = 0;
free(LabelClass);				LabelClass				= 0;
free(CardLabelClass);			CardLabelClass			= 0;
//fclose(f);
//f = 0;

return ;

}
//########################################################################
/***********************************************************
 Basic statistical, input-output and matrix manipulation

 Authors. Peter Mueller, Stephen Morris, David Rossell
 Last modified. 12 2007
***********************************************************/

#if !defined(STAT_H_INCLUDED)
#define STAT_H_INCLUDED

#include <math.h>
#include <stdio.h>
#include <vector>
using namespace std;

#if !defined(M_PI)
#define M_PI (3.1415926535897932385)
#endif

#if !defined(M_PI_2)
#define M_PI_2 (1.570796326794896619231321691640)
#endif

#if !defined(LOG_M_2PI)
#define LOG_M_2PI (1.8378770664093453)
#endif

#if !defined(SQ_M_PI_2)
#define SQ_M_PI_2 (2.5066282746310002416123552)
#endif

#if !defined(LOG_M_PI)
#define LOG_M_PI (1.1447298858494)
#endif

#if !defined(SIGN)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif


/**************************************************************/
/* Functions to compute means & variances                     */
/**************************************************************/

double meani(int *x, int lim);
double wmeani(int *x, int lim, double *w);
double meanx(double *x, int lim); //mean of x[0..lim]
double wmeanx(double *x, int lim, double *w);
double vari(int *x, int lim, int unb);
double wvari(int *x, int lim, double *w);
double varx(double *x, int lim, int unb);
double wvarx(double *x, int lim, double *w);
double cv(double *x, int ini, int fi); //coefficient of variation i.e. SD/mean
double cvinv(double *x, int ini, int fi); //coefficient of variation of 1/x

void colMeans(double *m, double *x, int nrow, int ncol);
void colVar(double *m, double *x, int nrow, int ncol);
void colCV(double *cv, double *x, int nrow, int ncol);
void colCVinv(double *cv, double *x, int nrow, int ncol); //CV of 1/x


/************************************************************************
                         BASIC BAYESIAN MODELS
************************************************************************/

void nn_bayes(double *mpo, double **Spo, double **Spo_inv, int p, double r1, double *mpr, double **Spr_inv, double r2, double *y, double **Slik_inv);  //Posterior of multiv normal mean with normal prior
void nn_bayes_rand(double *theta, int p, double r1, double **Spr_inv, double *mpr, double r2, double **Slik_inv, double *y); //Single draw from posterior of multiv normal mean with normal prior
double nn_integral(double *x, double *rx, double **Vxinv, double *detVx, double *mpr, double *rpr, double **Vprinv, double *detVpr, int *p, int *logscale); //Normal-Normal integral (useful to compute Bayes factors etc.)

void lm (double *b, double **XtX, double **invXtX, double *Xty, double *s, double *ypred, double *y, double **X, int *n, int *p, int *useXtX); //classical multiple linear regression
void lmbayes (double *bpost, double *spost, double *b, double **Vb, double *a_s, double *b_s, double **XtX, double **invXtX, double *Xty, int *B, double *y, double **X, int *n, int *p, int *useXtX, double *mpr, double **Spr_inv, double *tauprior, double *nu0, double *s0); //Bayesian multiple linear regression
void lmbayes_knownvar (double *bpost, double *b, double **Vb, double **XtX, double **invXtX, double *Xty, double *sigma, int *B, double *y, double **X, int *n, int *p, int *useXtX, double *mpr, double **Spr_inv, double *tauprior); //same as lmbayes with known variance sigma^2


/**************************************************************/
/* Input/output functions (interface)                         */
/**************************************************************/

FILE *openIn(char *);
FILE *openOut(char *);
//void scanFloat(char *, float *);
//void scanDouble(char *, double *);
//void scanInt(char *, int *);
//void fscanDouble(FILE *,char *, double *);
//void fscanInt(FILE *,char *, int *);
//void scanLong(char *, long *);
// 
//void scanFloatArray(char *,float *, int);
//void scanArray(char *,float *, int);
//void scanDoubleArray(char *,double *, int);
//void scanString(char *txt, char *s, int n);
//void fscanString(FILE *,char *txt, char *s, int n);
//void fscanDoubleArray(FILE *,double *, int);
//void scanDoubleMatrix(char *, double **, int, int);
//void fscanDoubleMatrix(FILE *ifile, double **x,int r,int c);
//void scanIntArray(char *,int *, int);
//void fscanIntArray(FILE *ifile, int *x, int n);

void writeInt(int);
void writeLong(long i);
void writeFloat(float);
void writeDouble(double);

void writeIntArray(int *,int, int);
void fwriteIntArray(FILE *, int *,int, int);
void fwriteIntMatrix(FILE *f, int **x, int rows, int cols);
void writeIntMatrix(int **x, int rows, int cols);
void writeDoubleArray(double *,int, int);
void writeDoubleMatrix2(double **, int , int);
void fwriteDoubleArray(FILE *, double *,int, int);
void fwriteDoubleMatrix2(FILE *, double **, int , int);
void writeDoubleMatrix(double **, int, int);
void writeFloatArray(float *, int, int);
void writeArray(float *, int, int); 

void fserror(char *proc, char *act, char *what);

/**************************************************************/
/* Debug messages etc. (mess)                                 */
/**************************************************************/

void errorC(char *,char *, int);
void err_msg(char *fct, char *txt, int n1, int n2, int n3);


/****************************************************************/
/* C Memory allocation                                          */
/****************************************************************/

float   *vector2(int,int);
double  *dvector(int,int);
double  **dmatrix(int,int,int,int);
//double  ***darray_3(int, int);
//double ***darray3(int n, int p, int q);
double ***darray3(int n1,int n2,int n3);   //allocate 3-way double array [0..n1-1][0..n2-1][0..n3-1]
int     *ivector(int,int);
int     **imatrix(int,int,int,int);
//int ***iarray_3(int lo, int hi);
//int ***iarray3(int p1, int p2, int p3);
int ***iarray3(int n1,int n2,int n3);   //allocate 3-way int array [0..n1-1][0..n2-1][0..n3-1]

void free_vector(float  *,int,int);
void free_dvector(double  *,int,int);
void free_ivector(int  *,int,int);
void free_dmatrix(double  **,int,int,int,int);
void free_imatrix(int  **,int,int,int,int);
void free_darray3(double ***a, int n1, int n2, int n3);
void free_iarray3(int ***a, int n1, int n2, int n3);

void nrerror(char *proc, char *act, char *what);

/**************************************************************/
/* Mathematical functions                                     */
/**************************************************************/

double gamln(double*);  //log-Gamma function
double gamln1(double*);  //auxiliary function called by gamln
double ldoublefact(double x);  //log-double factorial(x)
double digamma(double x);                          //from S Poetry (by Patrick J. Burns)
double trigamma(double x);
double polygamma(double x, long n, double low, double high, long terms, double nfact); //from S Poetry
double lnbeta(double a, double b); //log of Beta function
double betacf(double a, double b, double x); //continued fraction for incomplete Beta function
double lnchoose(double n, double k);
double choose(double n, double k);

double logit(double x);
double ilogit(double x);

double dsign(double x);  //returns 1.0 if x>=0, -1.0 if x<0
double isign(int x); //returns 1.0 if x>0, 0 if x==0, -1.0 if x<0

/**************************************************************/
/* Vector algebra (vector)                                    */
/**************************************************************/

void grid (double x0, double xn, int n, double *x);
void rA(double r,double **A, double **B, int rowini, int rowfi, int colini, int colfi);  //matrix*scalar
void A_plus_B(double **A, double **B, double **C, int rowini, int rowfi, int colini, int colfi); //matrix + matrix
void rA_plus_sB(double r,double **A,double s,double **B,double **C,int rowini,int rowfi,int colini,int colfi); //matrix*scalar + matrix*scalar
void rAx_plus_sBy(double r, double **A, double *x, double s, double **B, double *y, double *z, int rowini, int rowfi, int colini, int colfi); //scalar*matrix*vector + scalar*matrix*vector
void Ax_plus_y(double **A, double *x, double *y, double *z, int ini, int fi); //matrix*vector+vector
void xA(double *x,double **A,double *z, int ini, int fi);  //Multiply vector * matrix
void Ax(double **A,double *x,double *z,int rowini,int rowfi,int colini,int colfi);  //matrix * vector
void Avecx(double *A, double *x, double *z, int rowini, int rowfi, int colini, int colfi); //same but A is in vector format
void Atvecx(double *A, double *x, double *z, int rowini, int rowfi, int colini, int colfi); //same for A' (row/col indexes refer to A')
double xtAy(double *x, double **A, double *y, int ini, int fi); //t(vector)*matrix*vector

double quadratic_xtAx(double *x, double **A, int ini, int fi); //t(vector)*matrix*vector for quadratic forms (A symmetric)
double quadratic_xseltAselxsel(double *x, double *A, int *ncol, int *nsel, int *sel); // same but A is formatted as vector & only a subset of x,A is to be used
double quadratic_xtAselx(double *x, double *A, int *ncolA, int *nsel, int *sel); //same but subset is only for A
double quadratic_xseltAxsel(double *x, double **A, int ini, int *nsel, int *sel); //same but subset is only for x

void Atx(double **A,double *x,double *z,int rowini,int rowfi,int colini,int colfi); //t(matrix)*vector
void AB(double **A, int rowiniA, int rowfiA, int coliniA, int colfiA, double **B, int rowiniB, int rowfiB, int coliniB, int colfiB, double **C); //matrix*matrix, stored in C
void AtB(double **A, int rowiniA, int rowfiA, int coliniA, int colfiA, double **B, int rowiniB, int rowfiB, int coliniB, int colfiB, double **C); //t(matrix)*matrix, stored in C
void a_plus_b(double *a, double *b, double *c, int ini, int fi); //Vector sum i.e. c[i]=a[i]+b[i]
void a_prod_b(double *a, double *b, double *c, int ini, int fi); //Vector prod i.e. c[i]=a[i]*b[i]
void a_prod_b_sel(double *a, double *b, double *c, int *lengtha, int *nsel, int *sel); //same but only using indexes in sel
void a_zero(double *, int); //Set vector to zero
void R_zero(double **, int, int); //Set matrix to zero
void ddiag(double **A, int ini, int fi); //Diagonal matrix

int iabs(int x);   //absolute value of an integer
int imax_xy(int x, int y);
int imin_xy(int x, int y);
double max_xy(double x, double y);
double min_xy(double x, double y);
void minvec(double *x, int ini, int fi, double *xmin, int *minpos); //min of a vector and position at which min occurs
void maxvec(double *x, int ini, int fi, double *xmax, int *maxpos); //max of a vector and position at which max occurs

void choldc(double **a, int n, double **aout, bool *posdef);   //Cholesky decomposition
void choldc_inv(double **a, int n, double **aout, bool *posdef); //Inverse of Cholesky decomposition (input is original matrix)
void cholS_inv(double **cholS, int n, double **cholSinv); //Inverse of Cholesky decomposition (input is Cholesky matrix, faster)
void choldc_inv_internal(double **cholS, int n);
double choldc_det(double **chols, int n); //Determinant of a symmetric def+ using its Cholesky decomp
void inv_posdef(double **a, int n, double **aout, bool *posdef); //Inverse of a symmetric, positive definite matrix a[1..n][1..n]
void inv_posdef_upper(double **a, int n, double **aout, bool *posdef); //Same but only returns upper triangular elements
void invdet_posdef(double **a, int n, double **aout, double *det_a, bool *posdef); //Inverse and determinant of positive def matrix
void inv_posdef_chol(double **invchol, int n, double **aout); //Inverse given cholesky decomposition

void ludc(double **a, int n, int *indx, double *d); //LU decomposition 
void lu_solve(double **a, int n, int *indx, double b[]); //Solve A*x=b 
void lu_inverse(double **a, int n, double **aout); //Inverse of A[1..n][1..n]
double lu_det(double **a, int n); //Determinant of A[1..n][1..n]

void eigenvals(double **a, int n, double *vals);  //eigenvalues of a symmetric matrix
void tred2(double **a, int n, double d[], double e[], bool getVecs);
void tqli(double d[], double e[], int n, double **z, bool getVecs);
double pythag(double a, double b);

int dcompare(const void *a, const void *b);               
void dvecsort(double *v, int size);                           //sort a vector using qsort from stdlib
void dindexsort(double *x, int *index, int ilo, int ihi, int incr); //sort a vector of indexes using self-written quicksort routine
void iindexsort(int *x, int *index, int ilo, int ihi, int incr); //like dindexsort but for integers


/**************************************************************/
/* Random sampling                                            */
/**************************************************************/

void samplei_wr(int *x, int popsize, int n); //sample wo replacement from a vector of integers
void sampled_wr(double *x, int popsize, int n); //same for vector of doubles

/**************************************************************/
/* Probability distributions                                  */
/**************************************************************/

// Several
void setseed(long, long);
int rdisc(double *probs, int nvals);
double gamdev(double );
double dpoisson(int x, double mu, int logscale);

// Binomial, Multinomial
int rbinomial(int , double );
double dbinomial(int x, int n, double p, int logscale);
double dnegbinomial(int x, double r, double p, int logscale); //r: number of failures; p: success prob
void rmultinomial(int, int, double *, int *);

// Binomial related
double bbPrior(int k, int p, double alpha, double beta, int logscale); //Beta-Binomial, does not include the binomial coefficient
double dpoissonbin(int x, vector<double> *successProbs, vector<int> *nvars, int logscale, vector<double> *Tvector, vector<double> *poibinProbs); //Poisson-Binomial pmf. Result in poibinProbs[x]

// Uniform
double runif();
double dunifC(double x, double a, double b);
int runifdisc(int min, int max);

// Beta-Dirichlet
double rbetaC(double , double );
double pbetaC(double x, double pin, double qin); //quantile from a Beta(pin,qin)
void rdirichlet(double *w, double *alpha, int *p);
double ddirichlet(double *w, double *alpha, int *p, int *logscale); //Dirichlet density

// Normal
double dnormC(double y, double m, double s, int logscale); //density of Normal(m,s^2)
double dnormC_jvec(double *y, int n, double m, double s, int logscale); //joint density of y[0]...y[n-1] under Normal(m,s^2), i.e. returns scalar
double dmvnormC(double *y, int n, double *mu, double **cholsinv, double det, int logscale); //density of multivariate Normal
double	qnormC(double cdf, double m, double s);  //quantile from Normal(m,s^2)
double	pnormC(double y, double m, double s);  //cdf of Normal(m,s^2)
double rnormC(double mu, double s); //draw from univariate Normal(mu,s^2)
double rnorm_trunc(double ltrunc, double rtrunc, double m, double s); //draw trunc Normal given trunc points
double rnorm_trunc_prob(double lprob, double rprob, double m, double s); //draw trunc Normal given trunc probs
void rmvnormC(double *y, int n, double *mu, double **chols); //draw from multivariate Normal

// Moments
double mnorm(double order, double m, double sd); //raw moment of N(m,sd) of order "order"

// T Student
double dtC(double y, double mu, double s, int nu); //density of t with nu df
double dtmixC(double y, double *mu, double *s, double *probs, int nu, int ncomp, int logscale); //density of t_nu(mu[i],s[i]^2) mixtures with ncomp components
double dmvtC(double *y, int n, double *mu, double **cholsinv, double det, int nu, int logscale); //density of multivariate t
double rtC(int nu); //draw from univariate t with nu degrees of freedom
double rtmixC(double *mu, double *s, double *probs, int nu, int ncomp); //draw from mixture of t_nu(mu[i],s[i]^2)
double rt_trunc(int nu, double ltrunc, double rtrunc); //draw from truncated t given trunc points
double rt_trunc_prob(int nu, double lprob, double rprob);  //draw from truncated t given trunc probs
double qtC(double p, int nu);  //quantile from t-Student with nu degrees of freedom
double ptC(double x, int nu);  //CDF of t-Student with nu degrees of freedom
void rmvtC(double *y, int n, double *mu, double **chols, int nu); //draw from multivar T with nu degrees of freedom

// Gamma & Inverse gamma
double rgammaC(double a, double b); //a: shape; b: location; mean=a/b
double dgammaC(double x, double a, double b); //a: shape; b: location; mean=a/b
double dinvgammaC(double x, double a, double b); //a: shape; b: location; mean of x= b/(a-1)

// Non-local priors
double dmomNorm(double y, double m, double tau, double phi, int r, int logscale); //Normal MOM prior (power is 2*r)

/* More random variate stuff (dcdflib, from CMU statlib "www.stat.cmu.edu") */
double fifdint(double);
void cdfnor(int*,double*,double*,double*,double*,double*,int*,double*);
double spmpar(int*);
void cumnor(double*,double*,double*);
double dinvnr(double *p,double *q);
double stvaln(double*);
double devlpl(double [],int*,double*);
extern int ipmpar(int*);                      /* code in ipmpar.c */

/*even more stuff (ranlib) */
extern double genunf(double low,double high);
extern double gengam(double a,double r);
extern double sgamma(double a);
extern double snorm(void);
double fsign( double num, double sign );
extern double sexpo(void);
extern long mltmod(long a,long s,long m);
extern double ranf(void);
extern void gscgn(long getset,long *g);
extern void setall(long iseed1,long iseed2);   /* code in com.c */
extern void initgn(long isdtyp);               /* code in com.c */
extern long ignlgi(void);                      /* code in com.c */
extern void inrgcm(void);                      /* code in com.c */


/**************************************************************/
/* Integration                                                */
/**************************************************************/

double midpnt(double (*func)(double), double a, double b, int n); //nth stage refinement of integral of func from a to b (evenly spaced in x)
double midinf(double (*funk)(double), double aa, double bb, int n); //nth stage refinement of integral of func from aa to bb (evenly spaced in 1/x)
double qromo(double (*func)(double), double a, double b, double (*choose)(double(*)(double), double, double, int)); //Romberg integr on open interval (a,b)

/**************************************************************/
/* Interpolation, extrapolation and splines                   */
/**************************************************************/

void polint (double xa[], double ya[], int n, double x, double *y, double *dy); //interpolates via polynomials
double bspline_singlex(double x, int j, int degree, double *knots); //jth B-spline basis eval at single value x
void bspline(double **W, double *x, int *nx, int *degree, double *knots, int *nknots); //B-spline basis eval at vector of values x
void bspline_vec(double *W, double *x, int *nx, int *degree, double *knots, int *nknots); //same as bspline but returns a vector, so that it can be called from R
void mspline(double **W, double *x, int *nx, int *degree, double *knots, int *nknots); //M-spline basis eval at vector of values x
void mspline_vec(double *W, double *x, int *nx, int *degree, double *knots, int *nknots); //same as mspline but returns a vector, so that it can be called from R


/**************************************************************/
/* Function optimization                                      */
/**************************************************************/

double univmin(double ax,double bx,double cx,double (*f)(double),double tol,double *xmin,int itmax); //univariate minim
double dunivmin(double ax,double bx,double cx,double (*f)(double),double (*df)(double),double tol,double *xmin,int itmax);
void minimize(double th[],double **dirini,int n,double ftol,int *iter,double *fret,double (*f)(double []),int itmax);//multivar minim
void dirmin(double p[], double xi[], int n, double *fret, double (*func)(double []), int itmax, double dirminEPS); //minim in 1 direction
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,double (*func)(double)); //find bracketing triplets

#endif

#define MAXLOCI 1000000
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* data definition- holds sample size of each species, segsites in each species, D, length of locus */
struct hkaData{
  int na, nb, length, lArray[30], sArray[30];
  double d, sa, sb;
  char name[81];
};

struct hkaParams{
  int *lengths;
};

double harmSum(int n);
double harmSumSquare(int n);
void getData(int argc, char *argv[]);
void usage();
double calculateHKA(gsl_vector *params,gsl_vector *chiSquares, gsl_vector *expSA, gsl_vector *expD);
int chopByWhite(char *in, char *outArray[], int outSize);
int setHKAFunction(const gsl_vector *x, void *p, gsl_vector *f);
int solveSystem(gsl_vector *params);

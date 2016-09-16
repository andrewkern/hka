/* hka chi-square stat estimation 
based on original Hudson, Kreitman, Aguade 1987 implementation.

The parameters of the model are calculated by solving the system of equations numerically.

A. D. Kern  6/2005
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "hka.h"
#include "ctype.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>



/* currently the maximum number of loci is at 10^7, and the maximum locus name is 80 characters long */
struct hkaData data[MAXLOCI];
int locusNumber = 0;
int lineageNumber = 1;

int main(int argc, char *argv[]){
  int i = 0;
  gsl_vector  *chiSquares, *expectedPolyA, *expectedD, *devs;
  gsl_vector  *params; //0-locusNumber-1 = 4Nu's, locusNumber = tHat, locusNumber + 1 = fHat
  double x;

  getData(argc, argv);
  if(lineageNumber > 1){
    params = gsl_vector_alloc(locusNumber + 2);
  }
  else{
    params = gsl_vector_alloc(locusNumber + 1);
  }

  //solve the system of equations
  solveSystem(params);
 
  //calculate chiSquare Stats
  chiSquares =gsl_vector_alloc(locusNumber);
  expectedPolyA = gsl_vector_alloc(locusNumber);
  expectedD = gsl_vector_alloc(locusNumber);
  devs =  gsl_vector_alloc(locusNumber);

  x =  calculateHKA(params, chiSquares, expectedPolyA, expectedD);

  //print some stuff
  if(lineageNumber > 1){
    printf("chiSquared = %f\nt_hat = %f\nf_hat = %f\n",x, gsl_vector_get(params,locusNumber), gsl_vector_get(params,locusNumber + 1));
  }
  else{
    printf("chiSquared = %f\nt_hat = %f\nf_hat = 1.0\n",x, gsl_vector_get(params,locusNumber));
  }
    for(i = 0; i < locusNumber; i++){
      printf("locus %s \ttheta:\t %f\texpSA:\t %f\tobsSA:\t %f\texpD:\t %f\tobsD:\t %f\t chiSquared:\t%f\n",data[i].name,   gsl_vector_get(params,i),gsl_vector_get(expectedPolyA,i), data[i].sa, gsl_vector_get(expectedD,i), data[i].d, gsl_vector_get(chiSquares,i));
  } 
      
  gsl_vector_free(params);
  return(0);
}

void usage(){
  printf("usage:\nhka infile\n");
  exit(1);
}


/* parses the data file and reads in options (currently no options boss) */
void getData(int argc, char *argv[]){
  FILE *infile;
  int  l, na, nb, sa, sb, args, flag, i,lineCount;
  double d;
  char string[81];  //here's the character length limit
  char *sArray[30], line[500];

  if (argc < 2){
    usage();
  }
  else{
    infile = fopen(argv[1],"r");
    if (infile == NULL){
      fprintf(stderr,"Error opening infile!!!\n");
      exit(1);
    }
    flag = 0;
    args = 2;
    while(args < argc){
      switch(argv[args][1]){
      case 'm' : //missing data mode
	flag = 1;
      }
      args++;
    }
    if (flag){
      //missing data mode expect- na, d,lArray, sArray 
      while(fgets(line, 100, infile)){
	lineCount = chopByWhite(line,sArray,30);
	//	data[locusNumber].name = (char) sArray[0];
	data[locusNumber].na = atoi(sArray[1]);
	data[locusNumber].d = (double) atof(sArray[2]);
	for(i = 0; i < data[locusNumber].na - 2; i++){
	  data[locusNumber].lArray[i] = atoi(sArray[i+3]);
	  data[locusNumber].sArray[i] = atoi(sArray[i+data[locusNumber].na - 1 + 3]);
	}
	
	printf("%d\n",data[locusNumber].sArray[2]);
      }
    }
    else{
      
      /* collect hka data from infile- expect length, na, nb, sa, sb, d */
      while(fscanf(infile, "%80s %d %d %d %d %d %lf", string,  &l, &na, &nb, &sa, &sb, &d) != EOF){
	strcpy(data[locusNumber].name, string);
	data[locusNumber].length = l;
	data[locusNumber].na = na;
	data[locusNumber].nb = nb;
	if (nb > 1){
	  lineageNumber = 2;
	}
	data[locusNumber].sa = (double) sa;
	data[locusNumber].sb = (double) sb;
	data[locusNumber].d = (double) d;
	data[locusNumber].lArray[0] = 666;
	locusNumber++;
      }
    }
    fclose(infile);
  }
}

/* harmSum returns the denominator portion of Watterson's estimator */                         
double harmSum(int n){
  int i;
  double sum = 0.0;

  for(i = 1; i < n; i++){
    sum += 1.0 /  i;
  }
  return(sum);
}

/* harmSumSquare is for calculating the variance of Watterson's estimator */
double harmSumSquare(int n){
  int i;
  double sum = 0.0;

  for(i = 1; i < n; i++){
    sum += 1.0 / ( i * i);
  }
  return(sum);
}

/*setFunction- this sets up the multiroot function for solving */
int setHKAFunction(const gsl_vector *x, void *p, gsl_vector *f){
  double sumS, trueS, sumD, trueD, yi, trueYi;
  double sumSb, trueSb;
  int i;

  /*set up all the equations. the first are the sum of segSites and sum of divergence respectively */
  if (lineageNumber == 1){
    sumS = yi = sumD = trueS = trueD =  0.0;
    for(i = 0; i < locusNumber; i++){
      sumS += data[i].length * harmSum(data[i].na) * gsl_vector_get(x,i);
      sumD += data[i].length * (gsl_vector_get(x,locusNumber) + 1) * gsl_vector_get(x,i); 
      trueS += data[i].sa;
      trueD += data[i].d;
    }
    gsl_vector_set(f,0,fabs(trueS - sumS));
    gsl_vector_set(f,1,fabs(trueD - sumD));
    //Sa+D eqns
    for(i = 0; i < locusNumber - 1; i++){
      trueYi = data[i].sa + data[i].d;
      yi = gsl_vector_get(x,i) *						\
	((data[i].length * gsl_vector_get(x,locusNumber) +		\
	  data[i].length + (data[i].length * harmSum(data[i].na))));
      gsl_vector_set(f,i+2,fabs(trueYi - yi));
    }
  } 
  else{
    sumS = sumSb = yi = sumD = trueS = trueSb = trueD =  0.0;
    for(i = 0; i < locusNumber; i++){
      sumS += data[i].length * harmSum(data[i].na) * gsl_vector_get(x,i);
      sumSb += data[i].length * harmSum(data[i].nb) * gsl_vector_get(x,i) * gsl_vector_get(x, locusNumber + 1);
      sumD += data[i].length * (gsl_vector_get(x,locusNumber) + ((1.0 + gsl_vector_get(x,locusNumber+1)) / 2.0)) * gsl_vector_get(x,i);
      trueS += data[i].sa;
      trueSb += data[i].sb;
      trueD += data[i].d;
    }
    gsl_vector_set(f,0,fabs(trueS - sumS));
    gsl_vector_set(f,1,fabs(trueSb - sumSb));
    gsl_vector_set(f,2,fabs(trueD - sumD));
    //Sa+Sb+D eqns
    for(i = 0; i < locusNumber - 1; i++){
      trueYi = data[i].sa + data[i].sb + data[i].d;
      yi = (data[i].length * (gsl_vector_get(x,locusNumber) + ((1.0 + gsl_vector_get(x,locusNumber+1)) / 2.0)) * gsl_vector_get(x,i)) \
	+ ( data[i].length * harmSum(data[i].na) * gsl_vector_get(x,i)) + \
	(data[i].length * harmSum(data[i].nb) * gsl_vector_get(x,i) * gsl_vector_get(x, locusNumber + 1));
      gsl_vector_set(f,i+3,fabs(trueYi - yi));
    }
  }
  return GSL_SUCCESS;
}

int solveSystem(gsl_vector *params){
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  gsl_vector *x;
  int status;
  size_t  iter = 0;
  void *point=NULL;
  size_t n;

  if (lineageNumber == 1){
    n = locusNumber + 1;
    x = gsl_vector_alloc(locusNumber + 1);
  }
  else{
    n = locusNumber + 2;
    x = gsl_vector_alloc(locusNumber + 2);
  }
 
  gsl_multiroot_function f = {&setHKAFunction, n, point};
  gsl_vector_set_all(x, 1);
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, n);
  gsl_multiroot_fsolver_set(s, &f, x);

  do{
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);
    if(status){
      break;
    }
    status = gsl_multiroot_test_residual(s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 100);
 
  gsl_vector_memcpy(params, s->x);
  gsl_vector_free(x);
  gsl_multiroot_fsolver_free (s);

  return 0;
}


/* calculateHKA- this is the main beast that calculates estimates of t_hat, f(the ratio of popn sizes),
 and theta's for each locus. it returns the chi-square statistic, and takes pointers to the other 
 relavent quantities (t_hat, etc.). */ 
double calculateHKA(gsl_vector *params,gsl_vector *chiSquares, gsl_vector *expectedPolyA, gsl_vector *expectedD){
  double  chiSquared, chiSquared_i;
  double expSA, expSB, varSA, varSB, expD, varD;
  double hsA, hsB;
  double devA, devB, devD,fHat;
  int i;

  if(lineageNumber > 1){
    fHat = gsl_vector_get(params,locusNumber+1);
  }
  else{
    fHat = 1.0;
  }
  chiSquared = 0.0;
  for(i=0; i < locusNumber; i++){
    hsA = harmSum(data[i].na);
    hsB = harmSum(data[i].nb);
    //calculate expectations and variances at the ith locus
    expSA = gsl_vector_get(params,i) * data[i].length * hsA;
    varSA = expSA + ((gsl_vector_get(params,i) * data[i].length) * (gsl_vector_get(params,i) * data[i].length) * harmSumSquare(data[i].na));
    expSB =  gsl_vector_get(params,i) * fHat * data[i].length * hsB;
    varSB = expSB + ((gsl_vector_get(params,i) * fHat * data[i].length) * (gsl_vector_get(params,i) * fHat * data[i].length) * harmSumSquare(data[i].nb)); // equals zero if no poly data from speciesB
    expD = (gsl_vector_get(params,i) * data[i].length) * (gsl_vector_get(params,locusNumber) + (0.5 * (1.0 + fHat)));
    varD = expD + (((gsl_vector_get(params,i) * data[i].length) * 0.5 * (1.0 + fHat)) * ((gsl_vector_get(params,i) *  data[i].length) * 0.5 * (1.0 + fHat)));
    
    //calculate deviations for species_a_theta, species_b_theta, and divergence. 
    // then add to chi-square stat
    devA = ((data[i].sa - expSA) *  (data[i].sa - expSA)) / varSA;
    devB = ((data[i].sb - expSB) *  (data[i].sb - expSB)) / varSB;  //careful here will not be number if no poly data for speciesB
    devD = ((data[i].d - expD) * (data[i].d - expD)) / varD;
    
    //is deviation for species B finite?
    if (isnan(devB)){
      chiSquared_i = devA + devD;
    }
    else{
      chiSquared_i = devA + devB + devD;
    }
    chiSquared += chiSquared_i;
    gsl_vector_set(chiSquares,i,chiSquared_i);
    gsl_vector_set(expectedPolyA, i, expSA);
    gsl_vector_set(expectedD, i, expD);
  }

  return(chiSquared);
}                    
                         
                         
/* took this out of Jim Kent's tree */                         
int chopByWhite(char *in, char *outArray[], int outSize)
/* Like chopString, but specialized for white space separators. */
{
int recordCount = 0;
char c;
for (;;)
    {
    if (outArray != NULL && recordCount >= outSize)
	break;

    /* Skip initial separators. */
    while (isspace(*in)) ++in;
    if (*in == 0)
	break;
    
    /* Store start of word and look for end of word. */    
    if (outArray != NULL)
	outArray[recordCount] = in;
    recordCount += 1;
    for (;;)
        {
        if ((c = *in) == 0)
            break;
        if (isspace(c))
            break;
        ++in;
        }
    if (*in == 0)
	break;
 
    /* Tag end of word with zero. */
    if (outArray != NULL)
	*in = 0;
    /* And skip over the zero. */
    in += 1;
    }
return recordCount;
}                         
                                                  


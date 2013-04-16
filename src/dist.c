/* SCCS @(#)dist.c  1.4  02/08/98 */
/*
** The four routines for dist splitting
*/

#include <stdio.h>
#include "rpart.h"
#include "rpartS.h"
#include "node.h"
#include "rpartproto.h"

static double *dsts;
/* static double *wts, *grandmean, *diffs, *tdiffs; */
static int *tsplit, *countn, *count;

int distinit(int n, double *y[], int maxcat, char **error, 
          double *parm, int *size, int who, double *wt)
    {
    if (who==1 && maxcat >0) {
    graycode_init0(maxcat);        
    tsplit= (int *) ALLOC(maxcat*2, sizeof(int));   
    count = (int *)ALLOC(maxcat, sizeof(int)); 
    countn = (int *)ALLOC(maxcat*maxcat, sizeof(int)); 
    dsts  = (double *)ALLOC(maxcat*maxcat, sizeof(double));  
        if (countn==0 || count==0 || dsts== 0) {
        *error = "Could not allocate memory in distinit";
        return(1);
        }
    }
    *size = 1;
    return(0);
    }

/*
** The dist evaluation function.  Return the mean and the ss.
*/

void distss(int n, double *y[],  double *value, double *risk, double *wt)
    {
    int i,k;
    double temp;
    temp = 0;
/*    twt = 0; */
    for (i=1; i<n; i++)  
    for (k=0; k<i; k++) 
        temp += *y[rp.n*k-k*(k+1)/2+i-k-1];
    temp = temp/n;    
    *value = temp;
    *risk = temp;
    }

/*
** The dist splitting function.  Find that split point in x such that
**  the sum of squares of y within the two groups is decreased as much
**  as possible.  It is not necessary to actually calculate the SS, the
**  improvement involves only means in the two groups.
*/

void dist(int n,    double *y[],  FLOAT *x,     int nclass, 
       int edge, double *improve, FLOAT *split, int *csplit, double myrisk, double *wt)
    {
    int i, j, k, kj;
    double temp, sumdiffs_sq;
    double left_sum, right_sum;
/*    double left_wt, right_wt;  */
    int left_n, right_n;
    double best, total;
    int direction = LEFT;
    int where = 0;

    right_n = n;
        
    if (nclass==0) {
    
    left_n=0;
    best=0;
    
    total=0;
    for (k=1; k<n; k++)
    for (j=0; j<k; j++) 
    total += *y[rp.n*j-j*(j+1)/2+k-j-1];
    total = total/n;    

    for (i=0; right_n>edge; i++) {
	    temp=0; sumdiffs_sq=0; left_n++;  right_n--;
		right_sum=0; left_sum=0;

		if (i==0) left_sum=0;
		else {
			for (k=1; k<=i; k++)
			for (j=0; j<k; j++)  
			left_sum += *y[rp.n*j-j*(j+1)/2+k-j-1];
			left_sum = left_sum/(i+1);
		}
		
		if (i==(n-1)) right_sum=0;
		else {
			for (k=i+2; k<n; k++) 
			for (j=i+1; j<k; j++) 
			right_sum += *y[rp.n*j-j*(j+1)/2+k-j-1];
			right_sum = right_sum/(n-i-1);
		}

        if (x[i+1] !=x[i] &&  left_n>=edge) {
	        temp = total-left_sum-right_sum;

        if (temp > best) {
            best = temp;
            where = i;
            if (left_sum > right_sum) direction = LEFT;
                      else    direction = RIGHT;
            }
        }
    }

    *improve =  best/ myrisk;
    if (best>0) {   /* found something */
        csplit[0] = direction;
        *split = (x[where] + x[where+1]) /2;
        }
    }

    else {
    
/*
**  Do the easy coding for now - gd
**  Take countn and dsts as square matrices and fold them over
**  Fix it up later !!! 
*/
    for (i=0; i<nclass; i++) {
        count[i] =0;
    for (j=0; j<nclass; j++) {
        countn[i+nclass*j] =0;
        dsts[i+nclass*j] =0;
    }
    }

    k = x[0]-1;
    count[k]++;

    for (i=1; i<n; i++) {
    k = x[i]-1;
    count[k]++;
    for (j=0; j<i; j++) {
        kj = x[j]-1;    
        countn[k+nclass*kj]++;
            dsts[k+nclass*kj] += *y[rp.n*j-j*(j+1)/2+i-j-1];       
    }
    }

    for (i=0; i<nclass; i++) 
    for (j=0; j<=i; j++) {
    if (i!=j) {
        countn[i+nclass*j]=countn[i+nclass*j]+countn[j+nclass*i];
            dsts[i+nclass*j]=dsts[i+nclass*j]+dsts[j+nclass*i];    
    }
    }

        for (i=0; i<nclass; i++) {
        if (count[i]==0) tsplit[i] = 0;
        else tsplit[i] = RIGHT;
    }

    total = 0;
    for (k=0; k<nclass; k++) 
        if (tsplit[k]!=0) {
        for (j=0; j<=k; j++) 
            if (tsplit[j]!=0) 
        total += dsts[k+nclass*j];
    }

    /*
    ** Now find the split that we want
    */

    best = 0;
    /*
    ** Insert gray code bit here
    */

/*  if (numclass==2) graycode_init2(nclass, count, rate);
**              else graycode_init1(nclass, count);
**
**     Just use graycode_init1 here -- gd
*/

    graycode_init1(nclass, count);

    while((i=graycode()) < nclass) {

/* item i changes groups */

    left_n =0;  right_n = 0;
    left_sum = 0; right_sum = 0; 

    if (tsplit[i]==LEFT)  tsplit[i]=RIGHT;
    else tsplit[i]=LEFT;
        
    for (k=0; k<nclass; k++) 
        if (tsplit[k]==LEFT) {
        for (j=0; j<=k; j++) 
            if (tsplit[j]==LEFT)   {        
            left_n += countn[k+nclass*j];
            left_sum += dsts[k+nclass*j]; 
            }
        }
        else if (tsplit[k]==RIGHT) {
        for (j=0; j<=k; j++) 
            if (tsplit[j]==RIGHT)   {       
            right_n += countn[k+nclass*j];
            right_sum += dsts[k+nclass*j];   
                }
        }

    left_n = (int) (sqrt(2*left_n+0.25)+0.5);   
    right_n = (int) (sqrt(2*right_n+0.25)+0.5); 
    
    if (left_n>=edge  &&  right_n>=edge) {
    temp = total/n - left_sum/left_n - right_sum/right_n;

        if (temp > best) {
                best=temp;
                if (left_sum > right_sum)
                for (j=0; j<nclass; j++) csplit[j] = tsplit[j];
            else
                for (j=0; j<nclass; j++) csplit[j] = -tsplit[j];
                }
        }
        }
    }
    *improve = best / myrisk;      /* % improvement */

  }


/* SCCS @(#)distpred.c  1.3  02/08/98 */
/*
** The error function for dist splitting
*/

double distpred(double *y, double *yhat)
    {
    int i;
    double temp;
    temp = 0;
    for (i=0; i<rp.n; i++) temp += (y[i] - yhat[i])*(y[i] - yhat[i]);
    return(temp);
    }



/* SCCS @(#)mrt.c   1.4  02/08/98 */
/*
** The four routines for mrt splitting
*/

#include <stdio.h>
#include "rpart.h"
#include "rpartS.h"
#include "node.h"
#include "rpartproto.h"

static double *mean, *grandmean, *wts, *diffs, *tdiffs;

static int *tsplit, *countn, *countwt;

int mrtinit(int n, double *y[], int maxcat, char **error,
          double *parm, int *size,    int who, double *wt)
    {
    grandmean  = (double *)ALLOC(rp.num_y, sizeof(double));
    if (grandmean == 0) {
        *error = "Could not allocate memory in mrtinit";
        return(1);
        }
    diffs  = (double *)ALLOC(rp.num_y, sizeof(double));
    if (diffs == 0) {
        *error = "Could not allocate memory in mrtinit";
        return(1);
        }
    tdiffs  = (double *)ALLOC(rp.num_y, sizeof(double));
    if (tdiffs == 0) {
        *error = "Could not allocate memory in mrtinit";
        return(1);
        }
    if (who==1 && maxcat >0) {
    graycode_init0(maxcat);
    tsplit= (int *) ALLOC(maxcat*2, sizeof(int));
    countn = (int *)ALLOC(3*maxcat, sizeof(int));
    countwt = (int *)ALLOC(3*maxcat, sizeof(int));
    mean  = (double *)ALLOC(maxcat*rp.num_y, sizeof(double));
        if (countn==0 || mean==0 ) {
        *error = "Could not allocate memory in mrtinit";
        return(1);
        }
    wts    = mean + maxcat;
/*  sums   = wts + maxcat;  */
    }
    *size = rp.num_y;
    return(0);
    }

/*
** The mrt evaluation function.  Return the mean and the ss.
*/

void mrtss(int n, double *y[],  double *value, double *risk, double *wt)
    {
    int i,j;
    double temp, twt;
    double ss, ssj;
    ss = 0;
    twt = 0;
/*
** move ss = 0 to total over all y vars -- gd.
*/
    for (i=0; i<n; i++) twt += wt[i];

    for (j=0; j<rp.num_y; j++) {
    ssj = 0;
    temp =0;

    for (i=0; i<n; i++) temp += y[i][j] * wt[i];
    grandmean[j] = temp/twt;

    for (i=0; i<n; i++) {
    temp = y[i][j] - grandmean[j];
    if (rp.dissim==1) ssj += temp * temp * wt[i];
    else if (rp.dissim==2) ssj += fabs(temp) * wt[i];
    }
    value[j] = grandmean[j];
/*
** drop out individual SS for now GD
**
**    value[j+rp.num_y] = ssj;
**
*/
    ss += ssj;
    }
    *risk = ss;
    }

/*
** The mrt splitting function.  Find that split point in x such that
**  the sum of squares of y within the two groups is decreased as much
**  as possible.  It is not necessary to actually calculate the SS, the
**  improvement involves only means in the two groups.
*/

void mrt(int n,    double *y[],  FLOAT *x,     int nclass,
       int edge, double *improve, FLOAT *split, int *csplit, double myrisk, double *wt)
    {
    int i,j,k;
    double sumdiffs_sq;
    double left_sum, right_sum;
    double left_wt, right_wt;
    int left_n, right_n;
    double best;
    double temp = 0;
    int direction = LEFT;
    int where = 0;

    /*
    ** Compute the grand mean, which will be subtracted from all of the
    **  data elements "on the fly".  This makes the hand calculator formula
    **  numerically stable.
    ** Also get the total n
    */

    right_n = n;
    right_wt = 0;
    for (i=0; i<n; i++) right_wt += wt[i];
    for (j=0; j<rp.num_y; j++) {
	    grandmean[j] = 0; tdiffs[j] = 0;
	    for (i=0; i<n; i++) grandmean[j] += y[i][j] * wt[i];
 	    grandmean[j] /= right_wt;
    }

    if (nclass==0) {
    right_sum = 0; left_sum = 0;
    left_n = 0; left_wt = 0;
    best = 0;

    for (i=0; right_n>edge; i++) {
        temp = 0; sumdiffs_sq = 0; left_n++;  right_n--;
        left_wt += wt[i]; right_wt -= wt[i];

        for (j=0; j<rp.num_y; j++) {
            diffs[j] = (y[i][j] - grandmean[j]) * wt[i];
            temp += diffs[j];
            tdiffs[j] += diffs[j];
            if (rp.dissim==1) sumdiffs_sq += tdiffs[j]*tdiffs[j];
            else if (rp.dissim==2) sumdiffs_sq += fabs(tdiffs[j]);
       }

       left_sum  +=temp;
       right_sum -=temp;

       if (x[i+1] != x[i] &&  left_n>=edge) {
       if (rp.dissim==1) temp = sumdiffs_sq/left_wt + sumdiffs_sq/right_wt;
       else if (rp.dissim==2) temp = 2.0*sumdiffs_sq;

	   if (temp > best) {
        best = temp;
        where = i;
        if (left_sum < right_sum) direction = LEFT;
                  else    direction = RIGHT;
       }

      }

    }

    *improve = best / myrisk;
    if (best>0) {   /* found something */
        csplit[0] = direction;
        *split = (x[where] + x[where+1]) /2;
        }
    }
    else {

    for (i=0; i<nclass; i++) {
        countn[i] = 0; countwt[i] = 0;
                for (j=0; j<rp.num_y; j++) mean[i+nclass*j] = 0;
        }

    for (i=0; i<n; i++) {
            k = x[i] -1;
        countn[k] ++;
                countwt[k] += wt[i];
            for (j=0; j<rp.num_y; j++) mean[k+nclass*j] += (y[i][j] - grandmean[j]) * wt[i];
        }

        for (i=0; i<nclass; i++) {
    for (j=0; j<rp.num_y; j++) {
       if (countwt[i]>0) mean[i+nclass*j] /= countwt[i];
     }
    }

    for (i=0; i<nclass; i++) {
        if (countn[i]==0) tsplit[i] = 0;
        else tsplit[i] = RIGHT;
    }

    /*
    ** Now find the split that we want
    */
    left_n = 0;  right_n = n;
    left_wt = 0;
        right_wt = 0;
        for (i=0; i<n; i++) right_wt += wt[i];
    best = 0;
    /*
    ** Insert gray code bit here
    */

/*  if (numclass==2) graycode_init2(nclass, countn, rate);
**              else graycode_init1(nclass, countn);
**
**     Just use graycode_init1 here -- gd
*/

    graycode_init1(nclass, countn);

    while((i=graycode()) < nclass) {

/* item i changes groups */

    if (tsplit[i]==LEFT) {
        tsplit[i]=RIGHT;
        right_n += countn[i];
        left_n -= countn[i];
        right_wt += countwt[i];
        left_wt -= countwt[i];
        for (j=0; j<rp.num_y; j++) tdiffs[j] += mean[i+nclass*j] * countwt[i];
        }
    else {
        tsplit[i]=LEFT;
        right_n -= countn[i];
        left_n += countn[i];
            right_wt -= countwt[i];
        left_wt += countwt[i];
        for (j=0; j<rp.num_y; j++) tdiffs[j] -= mean[i+nclass*j] * countwt[i];
    }

    if (left_n>=edge  &&  right_n>=edge) {
        left_sum = 0; right_sum = 0; sumdiffs_sq = 0;
        for (j=0; j<rp.num_y; j++) {
            left_sum += tdiffs[j];
            right_sum -= tdiffs[j];
        if (rp.dissim==1) sumdiffs_sq += tdiffs[j]*tdiffs[j];
        else if (rp.dissim==2) sumdiffs_sq += fabs(tdiffs[j]);
        }
        if (rp.dissim==1) temp = sumdiffs_sq/left_wt + sumdiffs_sq/right_wt;
        else if (rp.dissim==2) temp = 2.0*sumdiffs_sq;
        if (temp > best) {
                best = temp;
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


/* SCCS @(#)mrtpred.c   1.3  02/08/98 */
/*
** The error function for mrt splitting
*/

double mrtpred(double *y, double *yhat)
    {
    int j;
    double temp;
    temp = 0;
    if (rp.dissim==1) {
    for (j=0; j<rp.num_y; j++) temp += (y[j] - yhat[j])*(y[j] - yhat[j]);
    }
    else if (rp.dissim==2) {
    for (j=0; j<rp.num_y; j++) temp += fabs((y[j] - yhat[j]));
    }
    return(temp);
    }

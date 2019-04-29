/*  ffmygrid.c

   Usage in Matlab:  [k lx ly bx by] = ffmygrid(uu,vv,z,xi,yi,Lx,Ly,'Option')

   (uu,vv) is k space coordinates.
   z is the value at those points.
   xi,yi are Cartesian coordinates you wish to have data.
   Lx, Ly is convolution kernel size in grid.
   'Option' specifies the kind of kernel
   'KB' 2d separable Kaiser Bessel function
    this is based on 'a Fast Sinc Function Gridding Alogrithm
    for  Fourier Inversion in Computer Tomography'
      by J.K.O'sullivan, IEEE t-mi, Dec 1985

   lx ly bx by
      Sangwoo Lee, 03/02/2000
 */

#include "math.h"
#include "mex.h"
#include "omp.h"
#include "string.h"

/* function for the zero order modified Bessel funciton of the first kind
 */
double bessi0(double x) {
  double fans;
  double y, xx, ax, ans;
  x = fabs(x);
  xx = x;
  if ((ax = fabs(xx)) < 3.75) {
    y = x / 3.75;
    y *= y;
    ans =
        1.0 +
        y * (3.5156229 +
             y * (3.0899424 +
                  y * (1.2067492 +
                       y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
  } else {
    y = 3.75 / ax;
    ans = (exp(ax) / sqrt(ax)) *
          (0.39894228 +
           y * (0.1328592e-1 +
                y * (0.225319e-2 +
                     y * (-0.157565e-2 +
                          y * (0.916281e-2 +
                               y * (-0.2057706e-1 +
                                    y * (0.2635537e-1 +
                                         y * (-0.1647633e-1 +
                                              y * 0.392377e-2))))))));
  }
  fans = (double)ans;
  return fans;
}

/* This function returns the minimum value in an array of n elements*/
double min(double x[], int n) {
  int k;
  double min_x;

  min_x = x[0];
  for (k = 1; k <= n - 1; k++) {
    if (x[k] < min_x)
      min_x = x[k];
  }
  return min_x;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *xx, *yy, *xi, *yi, *zzr, *zzi, *lx, *ly;
  double *resultr, *resulti, *rlx, *rly, *rbx, *rby, Lx, Ly, Bx, By, Lxy;
  double *resultr_private, *resulti_private;
  double Kxstep, Kystep;
  double kernelVal;
  char *opt;
  int optlen;
  int i, j, k, counter = 0, starti, startj, endi, endj;
  float pi = 3.141592;
  int M, N, L, N1, N2;
  int flag = 0;
  double distx, disty, *weight, minxi, minyi;

  /* get input arguments */
  xx = mxGetPr(prhs[0]);
  yy = mxGetPr(prhs[1]);
  zzr = mxGetPr(prhs[2]);
  zzi = mxGetPi(prhs[2]);
  xi = mxGetPr(prhs[3]);
  yi = mxGetPr(prhs[4]);
  lx = mxGetPr(prhs[5]);
  ly = mxGetPr(prhs[6]);

  if (nrhs > 8 | nrhs < 6) {
    mexErrMsgTxt("\ninvalid number of input");
  }

  optlen = mxGetNumberOfElements(prhs[7]) + 1;
  opt = mxCalloc(optlen, sizeof(char));
  mxGetString(prhs[7], opt, optlen);

  /* Check for real or complex datatypes */
  if (mxIsComplex(prhs[0])) {
    mexErrMsgTxt("xx is not real");
  }
  if (mxIsComplex(prhs[1])) {
    mexErrMsgTxt("yy is not real");
  }
  if (mxIsComplex(prhs[2])) {
    flag = 1;
  }
  if (mxIsComplex(prhs[3])) {
    mexErrMsgTxt("xi is not real");
  }
  if (mxIsComplex(prhs[4])) {
    mexErrMsgTxt("yi is not real");
  }

  /* Check for double datatype. Single is not supported */
  if (!(mxIsDouble(prhs[0]))) {
    mexErrMsgTxt("xx is not double precision ");
  }
  if (!(mxIsDouble(prhs[1]))) {
    mexErrMsgTxt("yy is not double precision");
  }
  if (!(mxIsDouble(prhs[2]))) {
    mexErrMsgTxt("data is not double precision");
  }
  if (!(mxIsDouble(prhs[3]))) {
    mexErrMsgTxt("xi is not double precision");
  }
  if (!(mxIsDouble(prhs[4]))) {
    mexErrMsgTxt("yi is not double precision");
  }

  if ((L = mxGetM(prhs[3])) == 1) {
    L = mxGetN(prhs[3]);
  }
  if ((M = mxGetM(prhs[4])) == 1) {
    M = mxGetN(prhs[4]);
  }
  if ((N = mxGetM(prhs[0])) == 1) {
    N = mxGetN(prhs[0]);
  }
  if ((N1 = mxGetM(prhs[1])) == 1) {
    N1 = mxGetN(prhs[1]);
  }
  if ((N2 = mxGetM(prhs[2])) == 1) {
    N2 = mxGetN(prhs[2]);
  }

  if ((N != N1) | (N1 != N2)) {
    mexErrMsgTxt("input grid dimension does not match");
  }

  /*mexPrintf("\n%s",opt);
     mexPrintf("\nL=%d, M=%d, N=%d, N1=%d, N2=%d\n",L,M,N,N1,N2);
   */

  if (!(plhs[0] = mxCreateDoubleMatrix(M, L, mxCOMPLEX))) {
    mexErrMsgTxt("\nfailed to create a matrix");
  }
  if (!(plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL))) {
    mexErrMsgTxt("\nfailed to create a matrix");
  }
  if (!(plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL))) {
    mexErrMsgTxt("\nfailed to create a matrix");
  }
  if (!(plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL))) {
    mexErrMsgTxt("\nfailed to create a matrix");
  }
  if (!(plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL))) {
    mexErrMsgTxt("\nfailed to create a matrix");
  }
  if (!(plhs[5] = mxCreateDoubleMatrix(M, L, mxREAL))) {
    mexErrMsgTxt("\nfailed to create a matrix");
  }
  if (!(N == N1 & N == N2)) {
    mexErrMsgTxt("\ninput (x,y,z) is not matching");
  }

  resultr = mxGetPr(plhs[0]);
  resulti = mxGetPi(plhs[0]);
  rlx = mxGetPr(plhs[1]);
  rly = mxGetPr(plhs[2]);
  rbx = mxGetPr(plhs[3]);
  rby = mxGetPr(plhs[4]);
  weight = mxGetPr(plhs[5]);

  Kxstep = fabs(xi[1] - xi[0]);
  Kystep = fabs(yi[1] - yi[0]);

  Lx = Kxstep * (*lx);
  *rlx = Lx;

  Ly = Kystep * (*ly);
  *rly = Ly;

  /**rlx=Lx=Kxstep*(*lx);
     *rly=Ly=Kystep*(*ly);
   */

  /* *rbx=Bx=pi*Lx/2/fabs(xi[1]-xi[0]);
     *rby=By=pi*Ly/2/fabs(yi[1]-yi[0]);
     NOTE:  Lx/fabs() = *lx
   */

  *rbx = Bx = pi * 1.45 / 2 * (*lx);
  *rby = By = pi * 1.45 / 2 * (*ly);

  /* minxi=xi[0];
     minyi=yi[M-1];
     mexPrintf("min x = %g",minxi);
     mexPrintf("min y = %g",minyi);
   */

  /* minxi = min(xi,M);
     mexPrintf("min x = %g \n",minxi);
     minyi = min(yi,L);
     mexPrintf("min y = %g \n",minyi);*/

  minxi = -L / 4.0;
  minyi = -M / 4.0;
  /*  mexPrintf("min y = %g \n",minyi);
      mexPrintf("min x = %g \n",minxi);  */

  for (i = 0; i < L * M; i++) {
    resultr[i] = 0;
    resulti[i] = 0;
    weight[i] = 0;
  }

  if (!strcmp(opt, "KB") | (nrhs == 7)) {

#pragma omp parallel private(i, j, k, distx, disty, starti, startj, endi,      \
                             endj, kernelVal, resultr_private,                 \
                             resulti_private)
    {
      resultr_private = (double *)malloc(sizeof(double) * L * M);
      resulti_private = (double *)malloc(sizeof(double) * L * M);

      for (i = 0; i < L * M; i++) {
        resultr_private[i] = 0;
        resulti_private[i] = 0;
      }
#pragma omp for schedule(dynamic, 1)
      for (k = 0; k < N; k++) {
        starti = (int)ceil((xx[k] - Lx / 2 - minxi) / Kxstep);
        if (starti < 0) {
          starti = 0;
        }
        endi = (int)floor((xx[k] + Lx / 2 - minxi) / Kxstep);
        if (endi >= L) {
          endi = (L - 1);
        }
        /*mexPrintf("endi = %d \n",endi);*/
        for (i = starti; i < endi + 1; i++) {
          startj = (int)ceil((yy[k] - Ly / 2 - minyi) / Kystep);
          if (startj < 0) {
            startj = 0;
          }
          endj = (int)floor((yy[k] + Ly / 2 - minyi) / Kystep);
          if (endj >= M) {
            endj = (M - 1);
          }
          /* mexPrintf("yy = %g, start = %d, end = %d \n",yy[k],startj,endj);*/
          for (j = startj; j < endj + 1; j++) {
            /*    mexPrintf("%d %d \n",i,j);  */
            distx = fabs(xi[i] - xx[k]);
            /*  mexPrintf(" %g  ",distx); */
            disty = fabs(yi[j] - yy[k]);
            /*   mexPrintf(" %g  ",disty); */
            /*if (distx > (Lx / 2)) {
              mexPrintf("dx = %g, Lx = %g", distx, Lx);
            }
            if (disty > (Ly / 2)) {
              mexPrintf("dy = %g, Ly = %g", disty, Ly);
            }
            */
            /* These operations were repeated, so using a variable to avoid the
             * extra cycles */
            kernelVal = 1.0 / Lx / Ly *
                        bessi0(Bx * sqrt(1 - 4 * distx * distx / Lx / Lx)) *
                        bessi0(By * sqrt(1 - 4 * disty * disty / Ly / Ly));

            resultr_private[M * i + j] += zzr[k] * kernelVal;
            if (flag) {
              resulti_private[M * i + j] += zzi[k] * kernelVal;
            }
            /*       weight[M*i+j]+=1/Lx/Ly*bessi0(Bx*sqrt(1-4*distx*distx/Lx/Lx))*bessi0(By*sqrt(1-4*disty*disty/Ly/Ly));
               if ((i=44) && (j=0)){
               mexPrintf("dx = %g, dy = %g",distx, disty);}*/
          }
        }
      }

#pragma omp critical
      {
        for (i = 0; i < L * M; i++) {
          resultr[i] += resultr_private[i];
          resulti[i] += resulti_private[i];
        }
      }

      free(resultr_private);
      free(resulti_private);
    }

    /*    for(i=0;i<L;i++){
       for(j=0;j<M;j++){
           if(weight[M*i+j] >= 0.000001){
        resultr[M*i+j]=resultr[M*i+j]/weight[M*i+j];
        resulti[M*i+j]=resulti[M*i+j]/weight[M*i+j];
       }
       }
       }
     */

  } else {
    mexErrMsgTxt("\ninvalid kernel name");
  }

  mxFree(opt);
}

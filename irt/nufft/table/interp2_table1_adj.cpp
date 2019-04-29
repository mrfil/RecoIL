// interp2_table1_adj.c
// *adjoint* of 2D periodic interpolation using table lookup.
//
// forward direction: (for m = 1,...,M)
// f(t_m) = \double_sum_{k1,2=0}^{K1,2-1} c[k1,k2]
//	h_1( (t_m1 - k1) mod K1 ) h_2( (t_m2 - k2) mod K2 )
//
// adjoint direction: (for k1,2=0,...,K1,2-1) (note complex conjugate!)
// c[k1,k2] = \sum_{m=1}^M f(t_m)
//	h_1^*( (t_m1 - k1) mod K1 ) h_2^*( (t_m2 - k2) mod K2 )
//
// Interpolators h_1,2 are nonzero (and tabulated) for -J/2 <= t <= J/2.
//
// Copyright 2004-4-2 Jeff Fessler and Yingying Zhang, University of Michigan

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "def,table2.h"


/*
* interp2_table0_complex_per_adj()
*/
void interp2_table0_complex_per_adj(
double *r_ck,		/* [K1,K2] out */
double *i_ck,
const int K1,
const int K2,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,
const double *r_h2,	/* [J2*L2+1,1] in */
const double *i_h2,
const int J1,
const int J2,
const int L1,
const int L2,
const double *p_tm,	/* [M,2] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm;

	//Let's declare variables properly here. --AMC
	double t1, t2, fmr, fmi, p2, coef2r, coef2i, v1r, v1i, v2r, v2i, p1, coef1r, coef1i;
	int koff1, k2, jj1, jj2, n2, k2mod, k12mod, k1, n1, k1mod, kk;
	double *r_ck_private, *i_ck_private;
	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	i_h1 += ncenter1;
	}
	{
	const int ncenter2 = floor(J2 * L2/2.);
	r_h2 += ncenter2;
	i_h2 += ncenter2;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*sizeof(*i_ck));

#pragma omp parallel default(none) shared(p_tm, r_fm, i_fm, r_h1, i_h1,\
	r_ck, i_ck, r_h2, i_h2) private(r_ck_private, i_ck_private, \
		t1, t2, mm, jj1, jj2, fmr, fmi, koff1, k1, k2, kk, p1, p2, n1, n2, coef2r, coef2i, \
		 coef1r, coef1i, k1mod, k2mod, k12mod, v2r, v1r, v1i, v2i)
		 {
			 r_ck_private = (double *)calloc(K1*K2,sizeof(*r_ck_private));
			 i_ck_private = (double *)calloc(K1*K2,sizeof(*i_ck_private));

			 /* interp */
			 #pragma omp for schedule(dynamic)
  	 	 for (mm=0; mm < M; mm++) {

				 /* Trying to replace this mess with array notation below - Works!
				 t2 = p_tm[M];
				 t1 = *p_tm++;
				 fmr = *r_fm++;
				 fmi = *i_fm++;
				 */
				 t2 = p_tm[M+mm];
				 t1 = p_tm[mm];
				 fmr = r_fm[mm];
				 fmi = i_fm[mm];

				 koff1 = 1 + floor(t1 - J1 / 2.);
				 k2 = 1 + floor(t2 - J2 / 2.);

				 for (jj2=0; jj2 < J2; jj2++, k2++) {
					 p2 = (t2 - k2) * L2;
					 n2 = /* ncenter2 + */ iround(p2);
					 coef2r = r_h2[n2];
					 coef2i = i_h2[n2];
					 k2mod = mymod(k2, K2);
					 k12mod = k2mod * K1;

					 v2r = coef2r * fmr + coef2i * fmi;
					 v2i = coef2r * fmi - coef2i * fmr;
					 k1 = koff1;

					 for (jj1=0; jj1 < J1; jj1++, k1++) {
						 p1 = (t1 - k1) * L1;
						 n1 = /* ncenter1 + */ iround(p1);
						 coef1r = r_h1[n1];
						 coef1i = i_h1[n1];
						 k1mod = mymod(k1, K1);
						 kk = k12mod + k1mod; /* 2D array index */

						 r_ck_private[kk] += coef1r * v2r + coef1i * v2i;
						 i_ck_private[kk] += coef1r * v2i - coef1i * v2r;
					 } /* j1 */
				 } /* j2 */
  	 }
		 #pragma omp critical
		 {
		 	for (mm = 0; mm < K1*K2; mm++) {
				r_ck[mm] += r_ck_private[mm];
		 		i_ck[mm] += i_ck_private[mm];
		  }
	 	 }
		 free(r_ck_private);
		 free(i_ck_private);
	 } // End OMP parallel region
}


/*
* interp2_table0_real_per_adj()
* 2D, 0th-order, real, periodic
*/
void interp2_table0_real_per_adj(
double *r_ck,		/* [K1,K2] out */
double *i_ck,
const int K1,
const int K2,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *r_h2,	/* [J2*L2+1,1] in */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
const int flip2,
#endif
const int J1,
const int J2,
const int L1,
const int L2,
const double *p_tm,	/* [M,2] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	}
	{
	const int ncenter2 = floor(J2 * L2/2.);
	r_h2 += ncenter2;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1, jj2;
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	const int koff1 = 1 + floor(t1 - J1 / 2.);
	int k2 = 1 + floor(t2 - J2 / 2.);

	for (jj2=0; jj2 < J2; jj2++, k2++) {
		const double p2 = (t2 - k2) * L2;
		const int n2 = /* ncenter2 + */ iround(p2);
		double coef2r = r_h2[n2];
		const int wrap2 = floor(k2 / (double) K2);
		const int k2mod = k2 - K2 * wrap2;
		const int k12mod = k2mod * K1;

#ifdef Provide_flip
		if (flip2 && (wrap2 % 2))
			coef2r = -coef2r; /* trick: sign flip */
#endif

		const double v2r = coef2r * fmr;
		const double v2i = coef2r * fmi;
		int k1 = koff1;

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = /* ncenter1 + */ iround(p1);
		register double coef1r = r_h1[n1];
		const int wrap1 = floor(k1 / (double) K1);
		const int k1mod = k1 - K1 * wrap1;
		const int kk = k12mod + k1mod; /* 2D array index */

#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; /* trick: sign flip */
#endif

		r_ck[kk] += coef1r * v2r;
		i_ck[kk] += coef1r * v2i;
	} /* j1 */
	} /* j2 */
    }
}


/*
* interp2_table1_real_per_adj()
* 2D, 1st-order, real, periodic
*/
void interp2_table1_real_per_adj(
double *r_ck,		/* [K1,K2] out */
double *i_ck,
const int K1,
const int K2,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *r_h2,	/* [J2*L2+1,1] in */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
const int flip2,
#endif
const int J1,
const int J2,
const int L1,
const int L2,
const double *p_tm,	/* [M,2] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	}
	{
	const int ncenter2 = floor(J2 * L2/2.);
	r_h2 += ncenter2;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1, jj2;
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	const int koff1 = 1 + floor(t1 - J1 / 2.);
	int k2 = 1 + floor(t2 - J2 / 2.);

	for (jj2=0; jj2 < J2; jj2++, k2++) {
		const double p2 = (t2 - k2) * L2;
		const int n2 = floor(p2);
		const double alf2 = p2 - n2;
		register const double *ph2 = r_h2 + n2;
		double coef2r = (1 - alf2) * *ph2 + alf2 * *(ph2+1);
		const int wrap2 = floor(k2 / (double) K2);
		const int k2mod = k2 - K2 * wrap2;
		const int k12mod = k2mod * K1;

#ifdef Provide_flip
		if (flip2 && (wrap2 % 2))
			coef2r = -coef2r; /* trick: sign flip */
#endif

		const double v2r = coef2r * fmr;
		const double v2i = coef2r * fmi;
		int k1 = koff1;

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = floor(p1);
		const double alf1 = p1 - n1;
		register const double *ph1 = r_h1 + n1;
		register double coef1r = (1 - alf1) * *ph1 + alf1 * *(ph1+1);
		const int wrap1 = floor(k1 / (double) K1);
		const int k1mod = k1 - K1 * wrap1;
		const int kk = k12mod + k1mod; /* 2D array index */

#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; /* trick: sign flip */
#endif

		r_ck[kk] += coef1r * v2r;
		i_ck[kk] += coef1r * v2i;
	} /* j1 */
	} /* j2 */
    }
}

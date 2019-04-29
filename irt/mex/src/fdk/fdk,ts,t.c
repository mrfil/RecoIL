/*
* fdk,ts,t.c
* Threaded versions of FDK back-projection
* For detector index (t,s).
* Copyright 2008-10-09, Jeff Fessler, University of Michigan
*/
#include "defs-env.h"
#include "def,fdk.h"
#include "jf,thread1.h"

typedef struct {
	float *image;	// [nz nx ny] <- trick!
	const cbct_ig *ig; // image geometry
	const cbct_cg *cg; // cone-beam CT system geometry
	int na; // # of views
	cfloat	*proj;	// [nt ns na] <- trick! projection views
	cdouble	*beta;	// [na] source angles [radians]
} fdk_ts_s;


/*
* fdk_ts_back_init()
* interface routine for threaded versions
*/
static sof fdk_ts_back_init(void *in, cint id, cint nthread)
{
	fdk_ts_s *pa = (fdk_ts_s *) in;
	const cbct_ig *ig = pa->ig;
	const cbct_cg *cg = pa->cg;
	cint na = pa->na;
	cfloat *proj = pa->proj;
	cdouble *beta = pa->beta;
	(void) nthread;

	for (int ia=0; ia < na; ++ia, proj += cg->ns * cg->nt) // each view
		Call(fdk_ts_back1, (pa->image,
			ig->nx, ig->ny, ig->nz,
			ig->dx, ig->dy, ig->dz,
			ig->offset_x, ig->offset_y, ig->offset_z,
			ig->mask2, id + 1, // each thread does some voxels only
			cg->dso, cg->dsd, cg->dfs,
			cg->ns, cg->nt,
			cg->ds, cg->dt, cg->offset_s, cg->offset_t,
			proj, beta[ia]))
	Ok
}


/*
* fdk_ts_back_t()
* entry point for threaded FDK back-projector
*/
sof fdk_ts_back_t(
float *image,	// [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cint na, // # of views
cfloat	*proj,	// [nt ns na] <- trick! projection views
cdouble	*beta,	// [na] source angles [radians]
cint nthread, // # of threads
cint chat)
{
	fdk_ts_s st;
#define put(arg) st.arg = arg;
	put(image)
	put(ig)
	put(cg)
	put(na)
	put(proj)
	put(beta)
#undef put

	Bzero(image, ig->nx * ig->ny * ig->nz) // initialize image volume to 0

	Call(jf_thread1_top, (fdk_ts_back_init,
                NULL /* wrap up */, &st, nthread, Chat))
        Ok
}

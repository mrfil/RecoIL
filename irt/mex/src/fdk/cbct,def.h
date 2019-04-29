// cbct,def.h
// cone-beam CT definitions
// Copyright 2008-10-09, Jeff Fessler, University of Michigan

#ifndef jf_cbct_def_h
#define jf_cbct_def_h

#include "defs-env.h"

// initial version used "double" for scalars
#ifndef rscalar
#define rscalar float
#endif
#define crscalar Const rscalar

// image geometry
typedef struct {
	int nx; // image dimensions
	int ny;
	int nz;
	rscalar dx; // pixel size (can be negative)
	rscalar dy; // can be negative to cause flip
	rscalar dz; // cannot be negative
	rscalar offset_x; // center offset in pixels (usually 0)
	rscalar offset_y;
	rscalar offset_z;
	byte *mask2; // [nx ny] 2D support mask: 0 or 1 ... nthread
} cbct_ig;

// cone-beam geometry
typedef struct {
	rscalar dso; // distance from source to isocenter (infinity for ||)
	rscalar dsd; // distance from source to detector (infinity for ||)
	rscalar dfs; // distance from detector focal point to source
			// 0 for 3rd-gen CT, infinity for flat detector
	int ns; // # detector channels per row
	int nt; // # detector rows (along axial direction)
	rscalar ds; // horizontal (transaxial) ray spacing
	rscalar dt; // vertical (axial) ray spacing
	rscalar offset_s; // detector channel offset [pixels]
	rscalar offset_t; // detector vertical offset [pixels]
} cbct_cg;

// work space (one for each thread)
typedef struct {
	float *view; // [ns nt] or [nt ns]

	// for sf1:
	float *gamma; // [ns] azimuthal angles between each ray and center one
	float *vec; // [Max(ns,nt)]
	float *weight_s; // [ns]
	float *uppers; // [nx+1]
} cbct_work;


// cbct,mask2.c

extern sof cbct_mask_init(
byte *mask_int, // [nx ny] out: integers 0 or 1 ... nthread
cbyte *mask_bin, // [nx ny] in: binary, possibly null
cint nx, cint ny,
cint nthread,
cint chat);

typedef enum {
	cbct_back_error,
	cbct_back_zero, // zero before back-projection
	cbct_back_inc, // incrementing back-projection (no zeroing first)
} cbct_back_init;


// cbct,*,back.c

typedef sof cbct_any_back1_type(
float *image, // [nz nx ny] <- trick!
cint nx,
cint ny,
cint nz,
crscalar dx,
crscalar dy, // can be negative to cause flip
crscalar dz,
crscalar offset_x, // center offset in pixels (usually 0)
crscalar offset_y,
crscalar offset_z_shift, // offset_z - zshifts[ia]
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
cbyte mask_id, // 1 ... nthread
crscalar dso, // distance from source to isocenter
crscalar dsd, // distance from source to detector
crscalar dfs, // distance from focal point to source (0 or inf)
cint ns,
cint nt,
crscalar ds, // horizontal ray spacing
crscalar dt, // vertical ray spacing
crscalar beta, // source angle [radians]
crscalar offset_s, // channel offset [pixels]
crscalar offset_t, // vertical offset on detector [pixels]
cfloat *proj_in, // [ns*nt] <- trick! projection view at angle beta
ctruf i_is_ns_nt, // 1 if input is [ns nt] or 0 if [nt ns]
cbct_work *cw, // work space
cfloat scale, // scale input projection view by this factor before backprojecting
cint iz_start, // usually 0
cint iz_end); // do [iz_start, iz_end)

extern cbct_any_back1_type cbct_nn1_back1; // cbct,sf1,back.c
extern cbct_any_back1_type cbct_pd1_back1; // cbct,nn1,back.c
extern cbct_any_back1_type cbct_sf1_back1; // cbct,sf1,back.c


// cbct,any,back,t.c

extern sof cbct_any_back_t(
float *image, // [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cbct_work *cw,
cint na, // # of views
cfloat *proj, // [ns*nt na] <- trick! projection views
ctruf i_is_ns_nt, // 1 if [ns nt] or 0 if [nt ns]
cfloat *beta, // [na] source angles [radians]
cfloat *offset_s, // [na] possibly null
cfloat *offset_t, // [na] possibly null
cfloat *zshifts, // [na]
cint nthread, // # of threads
cint iblock, // for OS, which block
cint nblock, // for OS, # blocks
ctruf compact, // 0 (for now)
cchar *systype, // nn1 pd1 sf1 ...
cbct_back_init,
cfloat scale,
cint iz_start, // usually 0
cint iz_end, // do [iz_start, iz_end)
cint chat);


// cbct,*,proj.c

typedef sof cbct_any_proj1_type(
cfloat *image, // [nz nx ny] <- trick!
cint nx,
cint ny,
cint nz,
crscalar dx,
crscalar dy, // can be negative to cause flip
crscalar dz,
crscalar offset_x, // center offset in pixels (usually 0)
crscalar offset_y,
crscalar offset_z_shift, // offset_z - zshifts[ia]
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
crscalar dso, // distance from source to isocenter
crscalar dsd, // distance from source to detector
crscalar dfs, // distance from focal point to source (0 or inf)
cint ns,
cint nt,
crscalar ds, // horizontal ray spacing
crscalar dt, // vertical ray spacing
crscalar beta, // source angle [radians]
crscalar offset_s, // channel offset [pixels]
crscalar offset_t, // vertical offset on detector [pixels]
float *proj_out, // [ns*nt] <- trick! projection view at angle beta
ctruf o_is_ns_nt, // 1 if output should be [ns nt], 0 if [nt ns]
cbct_work *cw, // work space
cfloat scale); // scale output projections by this factor

extern cbct_any_proj1_type cbct_nn1_proj1; // cbct,nn1,proj.c
extern cbct_any_proj1_type cbct_pd1_proj1; // cbct,pd1,proj.c
extern cbct_any_proj1_type cbct_sf1_proj1; // cbct,sf1,proj.c


// cbct,any,proj,t.c

extern sof cbct_any_proj_t(
cfloat *image, // [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cbct_work *cw,
cint na, // # of views
float *proj, // [ns*nt na] <- trick! projection views
ctruf o_is_ns_nt, // 1 for [ns nt] or 0 for [nt ns]
cfloat *beta, // [na] source angles [radians]
cfloat *offset_s, // [na] possibly null
cfloat *offset_t, // [na] possibly null
cfloat *zshifts, // [na]
cint nthread, // # of threads
cint iblock, // for OS, which block
cint nblock, // for OS, # blocks
ctruf compact, // 0 (for now)
cchar *systype, // nn1 pd1 sf1 ...
cfloat scale,
cint chat);


// cbct,work.c

extern cbct_work *cbct_work_alloc(
const cbct_ig *ig,
const cbct_cg *cg,
cint nthread);

extern sof cbct_work_free(cbct_work *cw, cint nthread);

extern void cbct_view_transpose(
float *po, // [n2 n1]
cfloat *pi, // [n1 n2]
cint n1, cint n2);

extern cfloat *cbct_back_view_prep(
cfloat *proj, // [ns nt] or [nt ns]
ctruf i_ns_nt, // [ns nt] input?
ctruf o_ns_nt, // [ns nt] output?
float *work, // [ns*nt] work space
cint ns,
cint nt,
cfloat scale);

extern float cbct_pd1_scale(
const cbct_cg *cg,
const cbct_ig *ig,
cint chat);


// cbct,sf1,misc.c

extern void cbct_sf1_sort_taus(float *taus);
extern float cbct_sf1_weight_value_s(cfloat *taus, cfloat bc);

#endif // jf_cbct_def_h

// fdk,ts.c
// Feldkamp aka FDK backprojection for arc/flat detector.
// For detector index (t,s).
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

#include "defs-env.h"
#include "def,fdk.h"

void fdk_ts_help(void)
{
	printf("\n\
\n\
	image = function('fdk,ts,back', nx,ny,nz, dx,dy,dz, \n\
		offset_x, offset_y, offset_z, mask2, \n\
		dso, dsd, ds, dt, offset_s, offset_t, proj, beta, nthread)\n\
\n\
		image output is single [nz nx ny] <- trick!\n\
		nx,ny,nz: (int32) image size\n\
		dx,dy,dz: (double) voxel size\n\
		offset_x,_y,_z: (double) center offset in pixels (usually 0)\n\
		mask2: (uint8) [nx ny] 2D support mask\n\
		dso: (double) distance from source to isocenter\n\
		dsd: (double) distance from source to detector\n\
		dfs: (double) distance from focal point to source (0 or inf)\n\
		ds: (double) horizontal ray spacing\n\
		dt: (double) vertical ray spacing\n\
		offset_s: (double) channel offset [pixels]\n\
		offset_t: (double) vertical offset on detector [pixels]\n\
		nthread: (int32) # of processors\n\
		proj: (single) [nt ns na] (trick!) projection view for each beta\n\
		beta: (double) [na] source angle(s) [radians]\n\
\n");
}


//
// fdk_ts_back1()
// The FDK backprojection is *added* to the image, so the user must zero it!
//
sof fdk_ts_back1(
float *image, // [nz nx ny] <- trick!
cint nx,
cint ny,
cint nz,
cfloat dx, // voxel size
cfloat dy, // can be negative to cause flip
cfloat dz,
cfloat offset_x, // image volume center offset in pixels (usually 0)
cfloat offset_y,
cfloat offset_z,
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
cbyte mask_id, // 1 ... nthread
cfloat dso, // distance from source to isocenter
cfloat dsd, // distance from source to detector
cfloat dfs, // distance from focal point to source (0 or inf)
cint ns, // projection view dimensions
cint nt,
cfloat ds, // horizontal ray spacing (view sample spacing)
cfloat dt, // vertical ray spacing (view sample spacing)
cfloat offset_s, // channel offset [pixels]
cfloat offset_t, // vertical offset on detector [pixels]
cfloat *proj, // [nt ns] <- trick! projection view at angle beta
cfloat beta) // source angle [radians]
{
	cfloat wx = (nx-1)/2. + offset_x;
	cfloat wy = (ny-1)/2. + offset_y;
	cfloat wz = (nz-1)/2. + offset_z;
	cfloat ws = (ns-1)/2. + offset_s;
	cfloat wt = (nt-1)/2. + offset_t;
	cfloat sinb = sin(beta);
	cfloat cosb = cos(beta);

	truf is_arc = False;
	if (dfs == 0)
		is_arc = True;
	else if (!Isinf(dfs))
		Warn("dfs not done - junk!")

//	if (dz / dt < 0) Fail("need dz/dt > 0")

	// loop over pixels in top slice
	for (int iy = 0; iy < ny; ++iy) {
		cfloat yy = dy * (iy - wy);
	 for (int ix = 0; ix < nx; ++ix, image += nz) {
		if (mask2[ix + iy*nx] != mask_id) // each thread does its part
			continue;

		cfloat xx = dx * (ix - wx);
		cfloat xbeta = xx * cosb + yy * sinb;
		cfloat ybetas = dso - (-xx * sinb + yy * cosb);
		cfloat mag = dsd / ybetas;
		cfloat ss = is_arc ? (dsd * atan2(xbeta, ybetas))
				: (mag * xbeta);
		cfloat ss_bin = ss / ds + ws;
		cint is = (int) Floorf(ss_bin); // index of nearest neighbor in "s"

		if (is < 0 || is >= ns-1) // out of FOV
			continue;

		else {
			cfloat w2 = is_arc ? // fan-beam image domain weighting
			(Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta))) : Sqr(mag);
			cfloat wr = ss_bin - is; // horizontal bilinear
			cfloat wl = 1. - wr; // interpolation factors
			register float *pi = image;
			register cfloat *pp1 = proj + is * nt;
			register cfloat *pp2 = proj + (is+1) * nt;

			int iz_min, iz_max;
#if 1
			if (dz / dt < 0) {
				Fail("todo: not done")
				iz_min = Ceilf(wz - wt * dt / (mag * dz));
				iz_max = Ceilf(wz + ((nt-2 - wt) * dt) / (mag * dz));
			}
			else {
				iz_min = Ceilf(wz - wt * dt / (mag * dz));
				iz_max = Ceilf(wz + ((nt-1 - wt) * dt) / (mag * dz));
			}

			iz_min = Max(iz_min, 0);
			iz_max = Min(iz_max, nz);

#else
			iz_min = 0;
			iz_max = nz;
#endif

			pi += iz_min;
			for (int iz = iz_min; iz < iz_max; ++iz, ++pi) { // slice loop
				cfloat zz = dz * (iz - wz);
				cfloat tt = mag * zz;
				cfloat tt_bin = tt / dt + wt;
				cint it = (int) Floorf(tt_bin); // nearest nbr in "t"

#if 0
				if (it < 0 || it >= nt-1) // out of FOV
					Fail4("bug min,max=%d,%d iz=%d it=%d",
						iz_min, iz_max, iz, it)
#endif

				cfloat wu = tt_bin - it;
				cfloat wd = 1. - wu;
				cfloat p1 = wl * pp1[it]
					+ wr * pp2[it]; // interpolate
				cfloat p2 = wl * pp1[it+1]
					+ wr * pp2[it+1]; // horizontal

				// final vertical interpolation:
				*pi += w2 * (wu * p1 + wd * p2);
			}
		}
	 }
	}

	Ok
}

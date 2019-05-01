![](https://github.com/acerjanic/RecoIL/blob/master/Logo.png)

RecoIL 
========

Tools for prototyping non-Cartesian MR image reconstruction
-----------------------------------------------------------
    

Features
--------

- Contains methods for regularized SENSE iterative image reconstruction using the CG solver framework derived from the MIRT.
- Supports several methods for field correction including a fast Time Segmentation model (Sutton and Fessler, 2003) as well as a DFT model incorporating field correction.
- Contains methods for phase corrected SENSE reconstructions for reconstructing multishot data with motion induced phase error such as diffusion or magnetic resonance elastography (MRE).
- Has an API (recoInfo) for abstracting away the details of your raw data files with a wrapper for Siemens TWIX (.dat) and experimental support for [ISMRMRD](https://github.com/ismrmrd/ismrmrd).


Installation
------------

Install RecoIL by adding ```initializePaths.m``` to your MATLAB startup. You can find more information about how to setup your MATLAB path on startup [at the Mathworks website](https://www.mathworks.com/help/matlab/ref/startup.html).


Contribute
----------

- Issue Tracker: github.com/mrfil/recoil/issues
- Source Code: github.com/mrfil/recoil


Support
-------

RecoIL has been tested with MATLAB R2017b on MacOS X and Linux. It should be compatible with modern MATLAB versions, but problems have been known to crop up from time to time. If you find an incompatibility, please let us know.

There are a few MEX files that may or may not be compiled correctly for Windows. YMMV.

If you are having issues, please let us know!
Please open an issue at [github.com/mrfil/recoil](https://wwww.github.com/mrfil/recoil).


Acknowledgements
----------------

We use the object format of the [MIRT](https://web.eecs.umich.edu/~fessler/code) as well as the NUFFT and support functions developed by Jeff Fessler and his students at the University of Michigan. Our ```solve_pwls_pcg.m``` is derived from ```pwls_pcg_1.m```. 

The contrib directory contains a number of tools developed by others. Of note, we have incorporated Brian Hargreave's vds code for spiral generation for our starter example.


License
-------

The project is licensed under the NCSA/University of Illinois Open Source license.

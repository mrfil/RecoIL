.. _basic-page:

**************
Basic concepts
**************

RecoIL is a toolkit for prototyping iterative, non-Cartesian Magnetic Resonance Image (MRI) reconstruction algorithms. 
RecoIL provides an object-oriented framework for prototyping these algorithms in MATLAB.

What do we mean by reconstruction?
----------------------------------

In MR imaging, data is most commonly, especially in clinical imaging, acquired in the Fourier domain and reconstructed in via an inverse Fourier transform (IFT).
This kind of reconstruction is termed 'direct reconstruction' as data is transformed from the k-space domain into the image domain via one IFT.

Iterative image reconstruction
------------------------------

Iterative image reconstruction is an alternative framework for image reconstruction. 
Rather than attempting to reconstruct images in one step, iterative reconstructions create a problem formulation,
in terms of known data, problem constraints, such as prior knowledge, and a system matrix.

The system matrix is a linear matrix operator that captures the imaging matrix. 
The most basic system matrix in MRI is the Fourier matrix. A simple iterative reconstruction could be formulated as:

.. math: 
    \langle \alpha, \beta  \rangle 
    \in 
    \Biggl \lbrace 
    { 
    M,\text{ if } 
    {
        l(\underline{x}) = 
        \frac { p(\underline{x}|M ) } { p(\underline{x}|U) } 
        \geq
        \frac { p(U) }{ p(M) } }
    \atop 
    U, \text{ otherwise } 




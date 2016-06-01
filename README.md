fourDvel
========

Routine for inferring time-dependent, 3d velocity fields from geolocated displacement data

requirements
------------
- Fortran 90 compiler (only tested on gfortran)
- BLAS and LAPACK

installation
------------
   $ make 

optionally:

   $ make install [prefix=/install/directory]

If BLAS and LAPACK are not installed in your LD_LIBRARY_PATH: 

   $ make F90INCLUDE=-L/path/to/blas_lapack 

usage
-----
Copy fourDvel.cmd to a directory and fill out the necessary blanks. Then do:

   $ fourDvel fourDvel.cmd 




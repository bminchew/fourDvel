# fourDvel command file

Number of Scenes (-)                                        ::  3
Master Scene (-)                                            ::  1

Output folder                       :: .
Output file prefix                  :: fourDvel
Output null value                   :: -100

# Optional outputs
Output diag(G^T*W*G)^-1)            :: n
Output offdiag(G^T*W*G)^-1)         :: n
Output diag(G^T*G)^-1               :: n
Output offdiag(G^T*G)^-1            :: n
Output GDOP sqrt(tr((G^TG)^-1))     :: n
Output sqrt(tr((G^T*W*G)^-1))       :: n

Output velocity magnitude           :: n
Output average correlation          :: n

Estimate DEM correction             :: n

# Optional output bounding box (CSV list in order SNWE)
# if not give, domain is calculated to include all input data
# grid spacing is always taken to be master scene spacing
Output bounding box (CSV order SNWE) :: none

# Optional comma-separated list of periods (max list length = 25)
Sinusoidal period (days)            :: 0.5, 14.76

Regularize by angular frequency     :: n
Regularization parameter            :: 10
Regularization reference period for horizontal components :: none
Regularization reference period for vertical component    :: none

############################################################
### Parameters for generating synthetic data

Looks for synthetic data [x,y]               :: 10, 10

Include DEM error :: y
Magnitude of DEM error (meters) :: 50

Sinusoidal amplitudes (corresponding to periods above)   ::  0.2, 1.0
Sinusoidal phase shifts (multiples of pi/4) :: 1, 2

############################################################

Scene :: 1
   LOS displacement file            ::
   Samples in LOS displacement (-)  ::
   Displacement null value :: 0
   Displacement conversion factor :: 1

   Correlation file                          ::
   Correlation null value                    :: 0

   # LOS files must be bit (by-pixel) interleaved in order ENU
   LOS file        :: Icelnd_32005_12039_12040_HH.los.grd
   LOS null value  :: -100

   Upper left corner Latitude (deg)      ::  65.14143312
   Upper left corner Longitude (deg)     :: -19.39213940
   Latitude Spacing (deg/pixel)          :: -0.000055560
   Longitude Spacing (deg/pixel)         ::  0.000111100

   Number of looks                       ::  36

   Perpendicular baseline (m)            ::  0.0
   Platform altitude (m)                 :: 651625.5
   Flat-Earth incidence angle file (deg from vertical) ::
   Flat-Earth incidence null value :: 0

   # Times must be given in decimal days since some reference day
   Master scene acquisition time (days)  ::
   Subordinate scene acq. time (days)    ::

   # Optional azimuth information (next 3 lines)
   # Azimuth offsets must have the same units as LOS displacement
   Azimuth offset file   :: none
   Azimuth offset null value            :: 0
   Azimuth conversion factor :: 1
   Platform heading file (in degrees east of north) :: none
   Platform heading null value :: 0


End of command file

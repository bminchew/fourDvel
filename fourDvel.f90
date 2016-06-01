   program fourDvel

!***************************************************************************
!***************************************************************************
!**
!**   FILE NAME: fourDvel.f90
!**
!**   DATE WRITTEN: February 2015
!**
!**   PROGRAMMER:
!**                          Brent Minchew
!**                California Institute of Technology
!**                       bminchew@caltech.edu
!**
!**   FUNCTIONAL DESCRIPTION: Time-dependent, 3D velocity inversion
!**
!**   INPUTS:  Phase or LOS displacement, correlation, LOS file, + others
!**            (LOS must be pixel-by-pixel interleave in order ENU)
!**
!**   OUTPUTS: 1) Velocity field in component order and reference
!**                  frame as given in LOS file
!**            2) [optional] error estimation for each component
!**
!**   ROUTINES CALLED:  dgelsy, dgesvd, dgetri, dgetrf
!**
!**   NOTES: must be called with appropriate command file
!**
!**   UPDATE LOG:
!**
!**   Date Changed        Reason Changed
!**   ------------       ----------------
!**   Feb. 2015    Modified from original code: vector_disp.f
!**   Mar. 2015    Added bbox option to allow user to control output region
!**   May  2015    Added flat-Earth incidence angle support
!**   Sep  2015    Added heading angle file support
!**   Sep  2015    Added frequency regularization
!**   Sep  2015    Changed LAPACK routines for least squares solver
!**   Feb  2016    Memory map input files
!**
!**   TO DO:
!**
!**   - Replace heading input file by calculating heading from LOS
!**   - Replace correlation input with variance input
!**
!***************************************************************************
!***************************************************************************
   use iso_c_binding
   implicit none

   interface
   type(c_ptr) function mmap(addr,len,prot,flags,fildes,off) bind(c,name='mmap')
   use iso_c_binding
   integer(c_int), value :: addr
   integer(c_size_t), value :: len
   integer(c_int), value :: prot
   integer(c_int), value :: flags
   integer(c_int), value :: fildes
   integer(c_size_t), value :: off
   end function mmap

   integer(c_int) function munmap(addr,len) bind(c,name='munmap')
   use iso_c_binding
   type(c_ptr), value :: addr
   integer(c_size_t), value :: len
   end function munmap
   end interface

   character*900           cmdfile,otpref,otfold,eastnm,nornm,upnm,msenm,gtgnm,sinstr
   character*900           str,str2,label,dumname,vmagnm,numnm,acornm,gdopnm,bboxstr
   character*30            curdate
   character*8             ichar
   character*1             errstr
   character*2             periodchar
   character(len=900),dimension(2000)::  unwnm,cornm,losnm,aznm,incnm,hdgnm

   integer                 ii,jj,kk,mm,pp,zz,ierr,ios,ar,r,k,ns,iid,dumi,tp,tf,alaz
   integer                 numscenes,master,cnt,doline,errest,nonul,maxcol,vm,rnkcheck
   integer                 stv(13),today(3),now(3),indx(3),mseout,scncnt,eoc,getwdp
   integer                 info,lwork,nrhs,corout,gdopest,gtgout,gtgoff,ocoff,rnk
   integer                 numvars,numperiods,dodemcor,periodiid,colsdom,demon,f_len_los
   integer                 startcol,endcol,periodgtgiid,perioderriid,regularize
   integer                 use_mmap,statstat,usrbbox,rnkdgelsy,linecol,mmiid,f_len
   real*4                  dum1,dum2,lne,nulo,seps
   real*8                  bbox(4),latdom,londom,r2d,d2r,alpha,reffreqhor,reffreqvrt
   real*8                  pi,twopi,ddum1,radearth,rho0,amp_error,phase_error,rcond
   real*8                  mse,dcor,dscor,sdum,rnktol,sumerr,cosinc,sininc,deps
   real*8                  cose,cosn,cosu,sine,sinn,sinu,re_eq,re_pl,coslat,sinlat
   real*8                  deltime,dlos(3),delcos,delsin,periodicfreqs(25),periods(25)
   real*4,allocatable::    unw(:,:),cor(:,:),los(:,:),east(:),north(:),up(:)
   real*4,allocatable::    dumv1(:),dumv2(:),azoff(:,:),hdgval(:,:),ehdg(:),nhdg(:),peghdg(:)
   real*4,allocatable::    erre(:),errn(:),erru(:),nulaz(:),nulhdg(:),gdopvec(:),incangle(:,:)
   real*4,allocatable::    nulu(:),nulc(:),nullos(:),nulinc(:),unwconv(:),aziconv(:)
   real*4,allocatable::    nlooks(:),msev(:),vmag(:),numout(:),avecor(:)
   real*4,allocatable::    gtgeast(:),gtgnorth(:),gtgup(:),gtgoen(:),gtgonu(:),gtgoeu(:)
   real*4,allocatable::    gtgdem(:),gtgsin(:,:),errdem(:),errsin(:,:)
   real*4,allocatable::    ocoen(:),oconu(:),ocoeu(:),wdop(:),gtg(:,:),gtgc(:,:),dem(:)
   real*4,allocatable::    phze(:,:),phzn(:,:),phzu(:,:),ampe(:,:),ampn(:,:),ampu(:,:)
   real*8,allocatable::    latin(:),lonin(:),latspc(:),lonspc(:),latsth(:),lonest(:)
   real*8,allocatable::    bperp(:),plthgt(:),timemast(:),timeserv(:),dtempmat(:,:)
   real*8,allocatable::    dgmat(:,:),ddvec(:),dvmat(:,:),dxvec(:)
   real*8,allocatable::    ccoef(:),work(:),dgtg(:,:),dgtgcopy(:,:),gdopmat(:,:),drtr(:,:)
   integer,allocatable::   lines(:),cols(:),colsrt(:),azlogic(:),ipiv(:),hdglogic(:)
   integer,allocatable::   scol(:),bcol(:),latoff(:),lonoff(:),eline(:),sline(:)
   integer,allocatable::   elinesrt(:),logic(:),jpvt(:)
   integer,allocatable::   losid(:),unwid(:),corid(:),azid(:),incid(:),hdgid(:)     ! file ids

   !!! memory map variables
   type(c_ptr)        ::   cptr
   integer(c_int)     ::   c_munmapret
   integer(c_size_t)  ::   c_len,c_len_los,c_off
   !real*4,target,allocatable:: dumtar1(:)
   real*4,pointer::       dumpoint(:)

   integer,parameter :: PROT_READ=1
   integer,parameter :: MAP_PRIVATE=2

   character(len=11),parameter:: fu='unformatted'
   character(len=6), parameter:: di='direct'
   character(len=7), parameter:: re='replace'

   !!! setup arrays of pointers
   type ptrtype
       real*4,pointer :: p(:)
   end type ptrtype

   type(ptrtype),allocatable :: unwptr(:),corptr(:),losptr(:)
   type(ptrtype),allocatable :: aziptr(:),hdgptr(:),incptr(:)

   print*,' ';print*,' '
   if (iargc().ne.1) then
      print*,'  ~~~  4D Velocity Field Inversion  ~~~  '
      print*,' '
      print*,' Usage: fourDvel input_command_file '
      print*,' '
      print*,' '
      stop
   endif

   call getarg(1,cmdfile)
   open(3,file=cmdfile,status='old')

   rnktol = 1.d-7
   rnkcheck = 1
   ios = 0
   lne = 1.d0
   errest = 0
   pi = 4.d0*atan(1.d0)
   twopi = 2.d0*pi
   r2d = 180.d0/pi
   d2r = pi/180.d0
   numscenes = 5000
   ns = 0
   mseout = 0
   eoc = 0
   vm = 0
   alaz = 0
   corout = 0
   gdopest = 0
   gtgout = 0
   gtgoff = 0
   ocoff = 0
   getwdp = 0
   numvars = 3
   dodemcor = 0
   demon = 0
   usrbbox = 0
   periodicfreqs = 0.d0
   scncnt = 0
   bbox = -999.
   doline = -1
   re_eq = 6378137.0d0  !!! earth semi-major axis (m)
   re_pl = 6356752.3d0  !!! earth semi-minor axis (m)
   deps = epsilon(cosinc)
   seps = epsilon(nulo)
   regularize = 0
   alpha = 0.d0
   reffreqhor = -1.d0
   reffreqvrt = -1.d0
   c_off = 0
   use_mmap = 1

   do while (ios.eq.0.and.ns.le.numscenes.and.scncnt.le.numscenes.and.eoc.eq.0)
      read(3,'(A)',iostat=ios) str
      k = scan(str,'::')
      label = trim(adjustl(str(1:k-1)))
      str=str(k+2:)

      select case (label)

      case('Number of Scenes (-)')
         read(str,*,iostat=ios) numscenes
         if (numscenes.lt.3) then
            print*,'Minimum of 3 scenes are required to do the inversion. Given ',numscenes
            stop
         endif
         allocate(lines(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in lines'
         allocate(cols(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in cols'
         allocate(latin(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in latin'
         allocate(lonin(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in lonin'
         allocate(latspc(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in latspc'
         allocate(lonspc(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in lonspc'
         allocate(nlooks(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nlooks'
         allocate(nulc(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nulc'
         allocate(nulu(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nulu'
         allocate(nullos(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nullos'
         allocate(nulaz(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nulaz'
         allocate(nulhdg(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nulhdg'
         allocate(nulinc(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nulinc'
         allocate(azlogic(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in azlogic'
         allocate(hdglogic(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in hdglogic'
         allocate(ehdg(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in ehdg'
         allocate(nhdg(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nhdg'
         allocate(peghdg(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in peghdg'
         allocate(bperp(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in bperp'
         allocate(timemast(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in timemast'
         allocate(timeserv(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in timeserv'
         allocate(plthgt(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in plthgt'
         allocate(unwconv(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in unwconv'
         allocate(aziconv(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in aziconv'
         allocate(unwptr(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in unwptr'
         allocate(corptr(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in corptr'
         allocate(losptr(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in losptr'
         allocate(aziptr(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in aziptr'
         allocate(hdgptr(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in hdgptr'
         allocate(incptr(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in incptr'
         azlogic = 0
         hdglogic = 0
         bperp = 0.d0
         plthgt = -1.d0
         unwconv = 1.d0
         aziconv = 1.d0
         nulhdg = -9999.

      case('Master Scene (-)')
         read(str,*,iostat=ios) master
      case('Output file prefix')
         otpref = str
      case('Output folder')
         otfold = str
      case('Output error estimates')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            errest = 1
         endif
      case('Output diag(G^T*W*G)^-1)')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            errest = 1
         endif
      case('Output model covariance (G^TWG)^-1)')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            errest = 1
         endif
      case('Output offdiag(G^T*W*G)^-1)')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            ocoff = 1
         endif
      case('Output diag(G^TG)^-1')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gtgout = 1
         endif
      case('Output diag(G^T*G)^-1')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gtgout = 1
         endif
      case('Output offdiag(G^T*G)^-1')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gtgoff = 1
         endif
      case('Output offdiag(G^TG)^-1')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gtgoff = 1
         endif
      case('Output sqrt(tr((G^T*W*G)^-1))')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            getwdp = 1
         endif
      case('Output GDOP sqrt(tr((G^TG)^-1))')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gdopest = 1
         endif
      case('Output velocity magnitude')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            vm = 1
         endif
      case('Output average correlation')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            corout = 1
         endif

      case('Output null value')
         if (trim(adjustl(str)).eq.'NaN'.or.trim(adjustl(str)).eq.'nan') then
            lne  = 0.d0
            nulo = 0.d0
         elseif (trim(adjustl(str)).eq.'Inf'.or.trim(adjustl(str)).eq.'inf') then
            lne  = 0.d0
            nulo = 1.d0
         elseif (trim(adjustl(str)).eq.'-Inf'.or.trim(adjustl(str)).eq.'-inf') then
            lne  = 0.d0
            nulo = -1.d0
         else
            read(str,*,iostat=ios) nulo
            lne = 1.d0
         endif
         nulo = nulo/lne

      case('Check G rank')
        errstr = trim(adjustl(str))
        if (errstr.eq.'n'.or.errstr.eq.'N') then
           rnkcheck = 0
        endif

      case('Rank (normalized) tolerance')
         read(str,*,iostat=ios) rnktol

      case('Estimate DEM correction')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
             demon = 1
         endif

      case('Output bounding box (CSV order SNWE)')
         if (trim(adjustl(str)).ne.'none') then
             call parse_bbox(str,bbox,size(bbox),dumi)
             if (dumi.ne.size(bbox)) then
                print*,'Incomplete bounding box list. Given ',str
                usrbbox = 0
             else
                usrbbox = 1
             endif
         endif

      case('Memory map input files')
        errstr = trim(adjustl(str))
        if (errstr.eq.'n'.or.errstr.eq.'N') then
            use_mmap = 0
        endif

      case('Sinusoidal period (days)')
         sinstr = str
         if (trim(adjustl(str)).ne.'none') then
            call get_periodics(str,periods,size(periods),numperiods)
            periodicfreqs(:numperiods) = twopi/(periods(:numperiods))
         else
            numperiods = 0
         endif

      case('Regularize by angular frequency')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
             regularize = 1
         endif
      case('Regularization parameter')
         if (trim(adjustl(str)).ne.'none'.and.trim(adjustl(str)).ne.'None' &
                & .and.trim(adjustl(str)).ne.'NONE') then
             read(str,*,iostat=ios) alpha
         else
             alpha = 0.d0
         endif
      case('Regularization reference period for horizontal components')
         if (trim(adjustl(str)).ne.'none'.and.trim(adjustl(str)).ne.'None' &
                 & .and.trim(adjustl(str)).ne.'NONE') then
              read(str,*,iostat=ios) ddum1
              reffreqhor = twopi/ddum1
         else

              reffreqhor = -1.d0
         endif
      case('Regularization reference period for vertical component')
         if (trim(adjustl(str)).ne.'none'.and.trim(adjustl(str)).ne.'None' &
                 & .and.trim(adjustl(str)).ne.'NONE') then
              read(str,*,iostat=ios) ddum1
              reffreqvrt = twopi/ddum1
         else
              reffreqvrt = -1.d0
         endif


      case('Scene')
         read(str,*,iostat=ios) ns
      case('Samples in LOS displacement (-)')
         if (ns.gt.0) then
            read(str,*,iostat=ios) cols(ns)
            scncnt = scncnt + 1
         endif
      case('LOS displacement file')
         if (ns.gt.0) unwnm(ns) = str
      case('Correlation file')
         if (ns.gt.0) cornm(ns) = str
      case('LOS file')
         if (ns.gt.0) losnm(ns) = str
      case('Azimuth offset file')
         if (trim(adjustl(str)).ne.'none'.and.trim(adjustl(str)).ne.'None' &
               & .and.trim(adjustl(str)).ne.'NONE'.and.ns.gt.0) then
            aznm(ns) = str
            azlogic(ns) = 1
            alaz = 1
         endif
      case('Flat-Earth incidence angle file (deg from vertical)')
         if (trim(adjustl(str)).ne.'none'.and.trim(adjustl(str)).ne.'None' &
                & .and.trim(adjustl(str)).ne.'NONE'.and.ns.gt.0) then
            incnm(ns) = str
         endif
      !case('Peg heading (degrees east of north)')
      !     if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) peghdg(ns)
      case('Platform heading file (in degrees east of north)')
          if (trim(adjustl(str)).ne.'none'.and.trim(adjustl(str)).ne.'None' &
                & .and.trim(adjustl(str)).ne.'NONE'.and.ns.gt.0) then
              hdgnm(ns) = str
              hdglogic(ns) = 1
          endif
      case('Platform heading null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nulhdg(ns)
      case('Azimuth offset null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nulaz(ns)
      case('Flat-Earth incidence null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nulinc(ns)
      case('Upper left corner Latitude (deg)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) latin(ns)
      case('Upper left corner Longitude (deg)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) lonin(ns)
      case('Latitude Spacing (deg/pixel)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) latspc(ns)
      case('Longitude Spacing (deg/pixel)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) lonspc(ns)
      case('Correlation null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nulc(ns)
      case('Displacement null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nulu(ns)
      case('Displacement conversion factor')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) unwconv(ns)
      case('Azimuth conversion factor')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) aziconv(ns)
      case('LOS null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nullos(ns)
      case('Number of looks')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nlooks(ns)
      case('Perpendicular baseline (m)')
         if(ns.le.numscenes.and.ns.gt.0.and.trim(adjustl(str)).ne.'None'.and. &
            & trim(adjustl(str)).ne.'none') read(str,*,iostat=ios) bperp(ns)
      case('Platform altitude (m)')
         if(ns.le.numscenes.and.ns.gt.0.and.trim(adjustl(str)).ne.'None'.and. &
            & trim(adjustl(str)).ne.'none') read(str,*,iostat=ios) plthgt(ns)
      case('Master scene acquisition time (days)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) timemast(ns)
      case('Subordinate scene acq. time (days)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) timeserv(ns)
      case('End of command file')
         eoc = 1
      case default
         continue
      end select
   enddo

!  set dem logical
   if (maxval(abs(bperp)).gt.1.d0) then
        dodemcor = demon
   endif

!  setup domain
   if (any(abs(bbox(1:2)).gt.90.d0).and.usrbbox.eq.1) then
        usrbbox = 0
        print *,'bounding box latitudes out of +/- 90 deg. bounds...given ',bbox(1:2)
   endif
   if (any(abs(bbox(3:4)).gt.360.d0).and.usrbbox.eq.1) then
        usrbbox = 0
        print *,'bounding box longitudes out of +/- 360 deg. bounds...given ',bbox(3:4)
   endif

   if (usrbbox.eq.0) then
        allocate(lonest(numscenes))
        lonest = lonin + lonspc*(cols-1)

        latdom = maxval(latin)
        londom = minval(lonin)
        colsdom = 1 + int(abs((londom-maxval(lonest))/lonspc(master)))
        doline = -1
        deallocate(lonest)
   else
        colsdom = 1 + int(abs((bbox(4)-bbox(3))/lonspc(master)))
        doline = 1 + int(abs((bbox(1)-bbox(2))/latspc(master)))
        latdom = max(bbox(1),bbox(2))
        londom = min(bbox(3),bbox(4))
   endif

!  determine necessary allocation widths
   maxcol = maxval(cols)
!  now allocate the rest
   allocate(bcol(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in bcol'
   allocate(scol(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in scol'
   allocate(eline(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in eline'
   allocate(sline(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in sline'
   allocate(latoff(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in latoff'
   allocate(lonoff(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in lonoff'
   allocate(losid(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in losid'
   allocate(unwid(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in unwid'
   allocate(corid(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in corid'
   allocate(ccoef(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in ccoef'
   allocate(logic(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in logic'
   allocate(dumv1(10*maxcol),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in dumv1'
   !allocate(dumtar1(10*maxcol),stat=ierr)
    !       if(ierr.ne.0) print*,'allocation error in dumtar1'
   allocate(dumv2(maxcol),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in dumv2'

   allocate(unw(numscenes,colsdom),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in unw'
   allocate(cor(numscenes,colsdom),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in cor'
   allocate(los(numscenes,3*colsdom),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in los'
   allocate(east(colsdom),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in east'
   allocate(north(colsdom),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in north'
   allocate(up(colsdom),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in up'

   allocate(numout(colsdom),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in numout'
   if (mseout.eq.1) allocate(msev(colsdom))
   if (vm.eq.1) allocate(vmag(colsdom))
   if (gdopest.eq.1) allocate(gdopvec(colsdom))
   if (corout.eq.1) allocate(avecor(colsdom))
   if (getwdp.eq.1) allocate(wdop(colsdom))
   if (gtgout.eq.1) then
      allocate(gtgeast(colsdom))
      allocate(gtgnorth(colsdom))
      allocate(gtgup(colsdom))
      if (numperiods.gt.0) then
         allocate(gtgsin(6*numperiods,colsdom))
      endif
      if (dodemcor.eq.1) then
         allocate(gtgdem(colsdom))
      endif
   endif
   if (errest.eq.1) then
      allocate(erre(colsdom))
      allocate(errn(colsdom))
      allocate(erru(colsdom))
      if (numperiods.gt.0) then
          allocate(errsin(6*numperiods,colsdom))
      endif
      if (dodemcor.eq.1) then
          allocate(errdem(colsdom))
      endif
   endif
   if (ocoff.eq.1) then
      allocate(ocoen(colsdom))
      allocate(oconu(colsdom))
      allocate(ocoeu(colsdom))
   endif
   if (gtgoff.eq.1) then
      allocate(gtgoen(colsdom))
      allocate(gtgonu(colsdom))
      allocate(gtgoeu(colsdom))
   endif
   if (alaz.eq.1) then
      allocate(azoff(numscenes,colsdom))
      allocate(hdgval(numscenes,colsdom))
      allocate(azid(numscenes))
      allocate(hdgid(numscenes))
   endif
   if (numperiods.gt.0) then
      allocate(phze(numperiods,colsdom))
      allocate(phzn(numperiods,colsdom))
      allocate(phzu(numperiods,colsdom))
      allocate(ampe(numperiods,colsdom))
      allocate(ampn(numperiods,colsdom))
      allocate(ampu(numperiods,colsdom))
   endif
   if (dodemcor.eq.1) then
      allocate(dem(colsdom))
      allocate(incangle(numscenes,colsdom))
      allocate(incid(numscenes))
      incangle = 90.d0
   endif

!  open output files
   otpref = adjustl(otpref)
   otfold = adjustl(otfold)
   tp = len_trim(otpref)
   if (otpref(tp:tp).eq.'.') tp = tp-1
   tf = len_trim(otfold)
   if (otfold(tf:tf).eq.'/') tf = tf-1

   call system('mkdir -p '//otfold(1:tf))
   call system('mkdir -p '//otfold(1:tf)//'/errors/')

   eastnm = otfold(1:tf)//'/'//otpref(1:tp)//'.east'
   nornm  = otfold(1:tf)//'/'//otpref(1:tp)//'.north'
   upnm   = otfold(1:tf)//'/'//otpref(1:tp)//'.up'
   numnm  = otfold(1:tf)//'/'//otpref(1:tp)//'.num'
   open(11,file=eastnm,access=di,form=fu,status=re,recl=4*colsdom)
   open(12,file=nornm,access=di,form=fu,status=re,recl=4*colsdom)
   open(13,file=upnm,access=di,form=fu,status=re,recl=4*colsdom)
   open(19,file=numnm,access=di,form=fu,status=re,recl=4*colsdom)

   if (errest.eq.1) then
      eastnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.obscov.east'
      nornm  = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.obscov.north'
      upnm   = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.obscov.up'
      open(14,file=eastnm,access=di,form=fu,status=re,recl=4*colsdom)
      open(15,file=nornm,access=di,form=fu,status=re,recl=4*colsdom)
      open(16,file=upnm,access=di,form=fu,status=re,recl=4*colsdom)
   endif
   if (mseout.eq.1) then
      msenm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.mse'
      open(17,file=msenm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   if (vm.eq.1) then
      vmagnm = otfold(1:tf)//'/'//otpref(1:tp)//'.mag'
      open(18,file=vmagnm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   if (corout.eq.1) then
      acornm = otfold(1:tf)//'/'//otpref(1:tp)//'.avecor'
      open(20,file=acornm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   if (gdopest.eq.1) then
      gdopnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.gdop'
      open(21,file=gdopnm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   if (gtgout.eq.1) then
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.gtgi.east'
      open(22,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.gtgi.north'
      open(23,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.gtgi.up'
      open(24,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   if (gtgoff.eq.1) then
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.gtgi.en'
      open(25,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.gtgi.nu'
      open(26,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.gtgi.eu'
      open(27,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   if (ocoff.eq.1) then
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.obscov.en'
      open(28,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.obscov.nu'
      open(29,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.obscov.eu'
      open(30,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   if (getwdp.eq.1) then
      gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.obscov.dop'
      open(31,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   if (dodemcor.eq.1) then
       gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.dem'
       open(32,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
   endif

   periodiid = 151     ! starting input file id
   iid = periodiid
   do ii=1,numperiods
       write(periodchar,'(i2)') ii
       gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.east'
       open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
       iid = iid + 1
       gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.north'
       open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
       iid = iid + 1
       gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.up'
       open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
       iid = iid + 1
       gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.east'
       open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
       iid = iid + 1
       gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.north'
       open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
       iid = iid + 1
       gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.up'
       open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
       iid = iid + 1
   enddo

   periodgtgiid = iid
   if (gtgout.eq.1) then
       do ii=1,numperiods
           write(periodchar,'(i2)') ii
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.gtgi.east'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.gtgi.north'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.gtgi.up'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.gtgi.east'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.gtgi.north'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.gtgi.up'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
       enddo
       if (dodemcor.eq.1) then
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.gtgi.dem'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
       endif
   endif

   perioderriid = iid
   if (errest.eq.1) then
       do ii=1,numperiods
           write(periodchar,'(i2)') ii
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.obscov.east'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.obscov.north'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinamp'//trim(adjustl(periodchar))//'.obscov.up'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.obscov.east'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.obscov.north'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.sinphz'//trim(adjustl(periodchar))//'.obscov.up'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
       enddo
       if (dodemcor.eq.1) then
           gtgnm = otfold(1:tf)//'/errors/'//otpref(1:tp)//'.obscov.dem'
           open(iid,file=gtgnm,access=di,form=fu,status=re,recl=4*colsdom)
           iid = iid + 1
       endif
   endif

!  set up each scene
   iid = iid + 100     ! starting input file id
   do ii=1,numscenes
      unwid(ii) = iid
      iid = iid + 1
      corid(ii) = iid
      iid = iid + 1
      losid(ii) = iid
      iid = iid + 1

      dumname = trim(adjustl(unwnm(ii)))
      call stat(dumname,stv,statstat)
      if (statstat.ne.0) then
         print*,'error using stat with file ',dumname
         print*,'system returned status ',statstat
         stop
      endif
      f_len = stv(8)
      f_len_los = 3*f_len
      c_len = f_len
      c_len_los = f_len_los
      lines(ii) = f_len/(4*cols(ii))

      if (use_mmap.eq.1) then
          dumname = trim(adjustl(unwnm(ii)))
          mmiid = unwid(ii)
          open(mmiid,file=dumname,access=di,form=fu,status='old',recl=f_len)
          cptr = mmap(0,c_len,PROT_READ,MAP_PRIVATE,fnum(mmiid),c_off)
          call c_f_pointer(cptr,unwptr(ii)%p,[c_len])
          close(mmiid)

          dumname = trim(adjustl(cornm(ii)))
          mmiid = corid(ii)
          open(mmiid,file=dumname,access=di,form=fu,status='old',recl=f_len)
          cptr = mmap(0,c_len,PROT_READ,MAP_PRIVATE,fnum(mmiid),c_off)
          call c_f_pointer(cptr,corptr(ii)%p,[c_len])
          close(mmiid)

          dumname = trim(adjustl(losnm(ii)))
          mmiid = losid(ii)
          open(mmiid,file=dumname,access=di,form=fu,status='old',recl=f_len_los)
          cptr = mmap(0,c_len_los,PROT_READ,MAP_PRIVATE,fnum(mmiid),c_off)
          call c_f_pointer(cptr,losptr(ii)%p,[c_len_los])
          close(mmiid)
      endif

      if (azlogic(ii).eq.1) then
         if (hdglogic(ii).ne.1) then
             print *,'No heading file given for ',aznm(ii),'. Stopping'
             stop
         endif

         azid(ii) = iid
         iid = iid + 1
         hdgid(ii) = iid
         iid = iid + 1

         azoff(ii,:) = nulaz(ii)
         hdgval(ii,:) = nulhdg(ii)
         !ehdg(ii) = sin(peghdg(ii)*d2r)   depreciated
         !nhdg(ii) = cos(peghdg(ii)*d2r)   depreciated

         if (use_mmap.eq.1) then
             dumname = trim(adjustl(aznm(ii)))
             mmiid = azid(ii)
             open(mmiid,file=dumname,access=di,form=fu,status='old',recl=f_len)
             cptr = mmap(0,c_len,PROT_READ,MAP_PRIVATE,fnum(mmiid),c_off)
             call c_f_pointer(cptr,aziptr(ii)%p,[c_len])
             close(mmiid)

             dumname = trim(adjustl(hdgnm(ii)))
             mmiid = hdgid(ii)
             open(mmiid,file=dumname,access=di,form=fu,status='old',recl=f_len)
             cptr = mmap(0,c_len,PROT_READ,MAP_PRIVATE,fnum(mmiid),c_off)
             call c_f_pointer(cptr,hdgptr(ii)%p,[c_len])
             close(mmiid)
         endif
      endif
      if (dodemcor.eq.1) then
         if (plthgt(ii).gt.0.and.abs(bperp(ii)).ne.0.d0) then
            incid(ii) = iid
         else
            incid(ii) = -iid
         endif
         iid = iid + 1

         if (use_mmap.eq.1) then
             dumname = trim(adjustl(incnm(ii)))
             mmiid = incid(ii)
             open(mmiid,file=dumname,access=di,form=fu,status='old',recl=f_len)
             cptr = mmap(0,c_len,PROT_READ,MAP_PRIVATE,fnum(mmiid),c_off)
             call c_f_pointer(cptr,incptr(ii)%p,[c_len])
             close(mmiid)
         endif
      endif

!  calculate pixel offsets for each scene relative to the master
      latoff(ii) = idnint((latin(ii) - latdom)/latspc(ii)) ! positive southward
      lonoff(ii) = idnint((lonin(ii) - londom)/lonspc(ii)) ! positive eastward

      eline(ii) = lines(ii) + latoff(ii)
      sline(ii) = -latoff(ii)

!  set column offset conditions
      if (lonoff(ii).le.0) then      ! shifted left or not at all relative to master
         bcol(ii) = 0                ! no buffer on left
         scol(ii) = -lonoff(ii) + 1  ! starting column
      else                       !  shifted right relative to master
         bcol(ii) = lonoff(ii)   ! number of null-valued columns to buffer on left
         scol(ii) = 1
      endif
      ccoef(ii) = 2.d0*nlooks(ii)
   enddo

   if (doline.lt.0) then
        allocate(latsth(numscenes))
        latsth = latin + latspc*(eline-1)
        doline = 1 + int(abs((latdom-minval(latsth))/latspc(master)))
        deallocate(latsth)
   endif

   if (alpha.le.0.d0) then
        regularize = 0
   endif

   print*,' ~~~~~      fourDvel      ~~~~~'
   print*,' '
   print*,' Total lines     = ',doline
   print*,' Total samples   = ',colsdom
   print*,' Total scenes    = ',numscenes
   print*,' Total sinusoids = ',numperiods
   if (dodemcor.eq.1) then
        print*,' DEM correction:   on'
   else
        print*,' DEM correction:   off'
   endif
   if (regularize.eq.1) then
        print*,' Regularization:   on'
        !print*,'     with alpha = ',alpha
   else
        print*,' Regularization:   off'
   endif
   print*,' ';print*,' '

    !  write ENVI-style header file
    tp = len_trim(otpref)
    if (otpref(tp:tp).eq.'.') tp = tp-1
    tf = len_trim(otfold)
    if (otfold(tf:tf).eq.'/') tf = tf-1
    call fdate(curdate)
    open(10,file=otfold(1:tf)//'/'//otpref(1:tp)//'.hdr',status='replace')
    write(10,"(a)") 'ENVI'
    write(10,"(a,a,a)") 'description = {Generated with fourDvel on ',trim(adjustl(curdate)),'}'
    write(ichar,'(i8)') colsdom
    write(10,"(a,a)") 'samples = ',adjustl(ichar)
    write(ichar,'(i8)') doline
    write(10,"(a,a)") 'lines = ',adjustl(ichar)
    write(10,"(a)") 'bands = 1'
    write(10,"(a)") 'header offset = 0'
    write(10,"(a)") 'file type = ENVI Standard'
    write(10,"(a)") 'data type = 4'
    write(10,"(a)") 'interleave = bsq'
    write(10,"(a)") 'byte order = 0'
    write(10,123) 'map info = {Geographic Lat/Lon, 1, 1, ',londom,', ',latdom,', ',&
               & abs(lonspc(master)),', ',abs(latspc(master)),', WGS-84}'
    write(10,"(a,f10.2)") '.output null value = ',nulo
    write(10,"(a,a)") '.sinusoidal periods (days) = ',trim(adjustl(sinstr))
    if (dodemcor.eq.1) then
       write(10,"(a,a)") '.DEM correction: on'
    else
       write(10,"(a,a)") '.DEM correction: off'
    endif
    if (regularize.eq.1) then
       write(10,"(a,a)") '.Regularization: on'
       write(10,"(a,e8.2)") '.   Regularization parameter = ',alpha
       if (reffreqhor.gt.0.d0) then
          write(10,"(a,f8.5)") '.   Regularization reference period, horizontal = ',twopi/reffreqhor
       endif
       if (reffreqvrt.gt.0.d0) then
          write(10,"(a,f8.5)") '.   Regularization reference period, vertical = ',twopi/reffreqvrt
       endif
    else
       write(10,"(a,a)") '.Regularization: off'
    endif
    close(10)

!  Square things away
   dumv2 = 0.d0
   numvars = 3 + 6*numperiods + dodemcor
   lwork = 10*numvars

!  last-minute memory allocations and deallocations
   deallocate(nlooks)
   deallocate(hdglogic)
   deallocate(peghdg)
   allocate(gtg(numvars,numvars))
   allocate(gtgc(numvars,numvars))
   allocate(dgtg(numvars,numvars))
   allocate(dgtgcopy(numvars,numvars))
   allocate(gdopmat(numvars,numvars))
   allocate(drtr(numvars,numvars))
   allocate(dxvec(numvars))
   allocate(jpvt(numvars))
   allocate(work(lwork))

   if (regularize.eq.1.and.numperiods.gt.0) then
      call load_Cmprior(drtr,alpha,periodicfreqs,numvars,numperiods,reffreqhor,reffreqvrt)
   else
      drtr = 0.d0
   endif

!  Let's do it
   do ii=1,doline
      if(mod(ii,1).eq.0.or.ii.eq.doline.or.ii.eq.1) then
         write(*,"(A1,A,t8,I8)",advance='no') achar(13),'Line = ',ii
      endif

      cnt = 0
      do mm=1,numscenes
         if (ii.ge.latoff(mm).and.ii.le.eline(mm).and.(sline(mm)+ii).gt.0) then
            linecol = (ii+sline(mm)-1)*cols(mm)

            if (use_mmap.eq.1) then
                dumv1 = nulu(mm)
                dumv1(bcol(mm)+1:bcol(mm)+cols(mm)) = unwptr(mm)%p(linecol+1:linecol+cols(mm))
                unw(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)

                dumv1 = nulc(mm)
                dumv1(bcol(mm)+1:bcol(mm)+cols(mm)) = corptr(mm)%p(linecol+1:linecol+cols(mm))
                cor(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)

                dumv1 = nullos(mm)
                dumv1(3*bcol(mm)+1:3*(bcol(mm)+cols(mm))) = losptr(mm)%p(3*linecol+1:3*(linecol+cols(mm)))
                los(mm,:) = dumv1(3*(scol(mm)-1)+1:3*(colsdom+scol(mm)-1))
            else
                dumv1 = nulu(mm)
                dumname = trim(adjustl(unwnm(mm)))
                open(unwid(mm),file=dumname,access=di,form=fu,status='old',recl=4*cols(mm))
                read(unwid(mm),rec=ii+sline(mm)) dumv1(bcol(mm)+1:bcol(mm)+cols(mm))
                close(unwid(mm))
                unw(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)

                dumv1 = nulc(mm)
                dumname = trim(adjustl(cornm(mm)))
                open(corid(mm),file=dumname,access=di,form=fu,status='old',recl=4*cols(mm))
                read(corid(mm),rec=ii+sline(mm)) dumv1(bcol(mm)+1:bcol(mm)+cols(mm))
                close(corid(mm))
                cor(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)

                dumv1 = nullos(mm)
                dumname = trim(adjustl(losnm(mm)))
                open(losid(mm),file=dumname,access=di,form=fu,status='old',recl=3*4*cols(mm))
                read(losid(mm),rec=ii+sline(mm)) dumv1(3*bcol(mm)+1:3*(bcol(mm)+cols(mm)))
                close(losid(mm))
                los(mm,:) = dumv1(3*(scol(mm)-1)+1:3*(colsdom+scol(mm)-1))
            endif

            if(azlogic(mm).eq.1) then
                if (use_mmap.eq.1) then
                    dumv1 = nulaz(mm)
                    dumv1(bcol(mm)+1:bcol(mm)+cols(mm)) = aziptr(mm)%p(linecol+1:linecol+cols(mm))
                    azoff(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)

                    dumv1 = nulhdg(mm)
                    dumv1(bcol(mm)+1:bcol(mm)+cols(mm)) = hdgptr(mm)%p(linecol+1:linecol+cols(mm))
                    hdgval(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)
                else
                    dumv1 = nulaz(mm)
                    dumname = trim(adjustl(aznm(mm)))
                    open(azid(mm),file=dumname,access=di,form=fu,status='old',recl=4*cols(mm))
                    read(azid(mm),rec=ii+sline(mm)) dumv1(bcol(mm)+1:bcol(mm)+cols(mm))
                    close(azid(mm))
                    azoff(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)

                    dumv1 = nulhdg(mm)
                    dumname = trim(adjustl(hdgnm(mm)))
                    open(hdgid(mm),file=dumname,access=di,form=fu,status='old',recl=4*cols(mm))
                    read(hdgid(mm),rec=ii+sline(mm)) dumv1(bcol(mm)+1:bcol(mm)+cols(mm))
                    close(hdgid(mm))
                    hdgval(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)
                endif
            endif
            if (dodemcor.eq.1) then
                if (incid(mm).gt.0) then
                    if (use_mmap.eq.1) then
                        dumv1 = nulinc(mm)
                        dumv1(bcol(mm)+1:bcol(mm)+cols(mm)) = incptr(mm)%p(linecol+1:linecol+cols(mm))
                        incangle(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)
                    else
                        dumv1 = nulinc(mm)
                        dumname = trim(adjustl(incnm(mm)))
                        open(incid(mm),file=dumname,access=di,form=fu,status='old',recl=4*cols(mm))
                        read(incid(mm),rec=ii+sline(mm)) dumv1(bcol(mm)+1:bcol(mm)+cols(mm))
                        close(incid(mm))
                        incangle(mm,:) = dumv1(scol(mm):colsdom+scol(mm)-1)
                    endif
                endif
            endif
            cnt = cnt + 1
         endif
      enddo

      east = nulo
      north = nulo
      up = nulo
      numout = nulo

      if (mseout.eq.1) msev = nulo
      if (vm.eq.1) vmag = nulo
      if (corout.eq.1) avecor = nulo
      if (gdopest.eq.1) gdopvec = nulo
      if (getwdp.eq.1) wdop = nulo
      if (gtgout.eq.1) then
         gtgeast = nulo
         gtgnorth = nulo
         gtgup = nulo
         if (numperiods.gt.0) then
             gtgsin = nulo
         endif
         if (dodemcor.eq.1) then
             gtgdem = nulo
         endif
      endif
      if (gtgoff.eq.1) then
         gtgoen = nulo
         gtgonu = nulo
         gtgoeu = nulo
      endif
      if (ocoff.eq.1) then
         ocoen = nulo
         oconu = nulo
         ocoeu = nulo
      endif
      if (errest.eq.1) then
         erre = nulo
         errn = nulo
         erru = nulo
         if (numperiods.gt.0) then
             errsin = nulo
         endif
         if (dodemcor.eq.1) then
             errdem = nulo
         endif
      endif
      if (numperiods.gt.0) then
         phze = nulo
         phzn = nulo
         phzu = nulo
         ampe = nulo
         ampn = nulo
         ampu = nulo
      endif
      if (dodemcor.eq.1) then
         dem = nulo
      endif

      if (cnt.ge.numvars) then
         if (dodemcor.eq.1) then
            !!! precompute parameters for dem correction
            coslat = cos((latdom+(ii-1)*latspc(master))*d2r)
            sinlat = dsqrt(1.d0 - coslat**2)
            radearth = dsqrt(((coslat*re_eq**2)**2+(sinlat*re_pl**2)**2)/ &
                       & (((coslat*re_eq)**2+(sinlat*re_pl)**2)))
         endif

         !!! this will narrow omp loop bound and make it more efficient when
         !!! there are large null-spaces bounding the image
         startcol = -1
         endcol = colsdom
         do jj=1,colsdom
             nonul = 0
             do mm=1,numscenes
                dum1 = product(los(mm,(1+(jj-1)*3):(3+(jj-1)*3))-nullos(mm))*cor(mm,jj) ! test null vals
                if(unw(mm,jj).ne.nulu(mm).and.dum1.ne.0.d0.and.cor(mm,jj).ne.nulc(mm)) then
                   nonul = nonul + 1
                   if (azlogic(mm).eq.1.and.azoff(mm,jj).ne.nulaz(mm).and.hdgval(mm,jj).ne.nulhdg(mm)) then
                        nonul = nonul + 1
                   endif
                endif
             enddo
             numout(jj) = nonul

             if (nonul.ge.numvars) then
                if (startcol.lt.0) then
                    startcol = jj
                endif
                endcol = jj
             endif
         enddo
         if (startcol.lt.1) then
             startcol = 1
         endif

         !$omp parallel do &
         !$omp default(private) &
         !$omp shared(los,cor,unw,ccoef,nullos,nulu,nulc,unwconv,aziconv,regularize) &
         !$omp shared(azlogic,azoff,nulaz,hdgval,nulhdg,ehdg,nhdg,startcol,endcol) &
         !$omp shared(east,north,up,erre,errn,erru,numout,mseout,msev,radearth) &
         !$omp shared(vm,vmag,corout,avecor,gdopest,gdopvec,getwdp,wdop,rnkcheck) &
         !$omp shared(gtgout,gtgeast,gtgnorth,gtgup,gtgoff,gtgoen,gtgonu,gtgoeu) &
         !$omp shared(phze,phzn,phzu,ampe,ampn,ampu,numperiods,periodicfreqs,drtr) &
         !$omp shared(timemast,timeserv,dem,dodemcor,numvars,plthgt,bperp,alpha) &
         !$omp shared(errest,ocoff,ocoen,oconu,ocoeu,cols,master,numscenes,rnktol) &
         !$omp shared(errsin,gtgsin,errdem,gtgdem,incangle,nulinc,pi,twopi,d2r,r2d) &
         !$omp shared(reffreqhor,reffreqvrt,lwork)
         do jj=startcol,endcol
            nonul = 0
            logic = 0
            do mm=1,numscenes
               dum1 = product(los(mm,(1+(jj-1)*3):(3+(jj-1)*3))-nullos(mm))*cor(mm,jj) ! test null vals
               if(unw(mm,jj).ne.nulu(mm).and.dum1.ne.0.d0.and.cor(mm,jj).ne.nulc(mm)) then
                  nonul = nonul + 1
                  logic(mm) = 1
                  if (azlogic(mm).eq.1.and.azoff(mm,jj).ne.nulaz(mm).and.hdgval(mm,jj).ne.nulhdg(mm)) then
                     nonul = nonul + 1
                  endif
               endif
            enddo

            if (nonul.ge.numvars) then
               allocate(dgmat(nonul,numvars))
               allocate(dvmat(nonul,nonul))
               allocate(ddvec(nonul))
               allocate(dtempmat(numvars,nonul))

               dgmat = 0.d0
               dvmat = 0.d0
               dscor = 0.d0

               kk = 1
               do mm=1,numscenes
                  if (logic(mm).eq.1) then
                     dcor = (dble(cor(mm,jj)))**2.d0
                     dvmat(kk,kk) = ccoef(mm)*dcor/(1-dcor)
                     dlos(1) = dble(los(mm,(1+(jj-1)*3)))
                     dlos(2) = dble(los(mm,(2+(jj-1)*3)))
                     dlos(3) = dble(los(mm,(3+(jj-1)*3)))

                     deltime = timemast(mm) - timeserv(mm)
                     dgmat(kk,1) = dlos(1)*deltime
                     dgmat(kk,2) = dlos(2)*deltime
                     dgmat(kk,3) = dlos(3)*deltime
                     ddvec(kk)   = dble(unw(mm,jj)*unwconv(mm))
                     dscor = dscor + dble(cor(mm,jj))

                     if (numperiods.gt.0) then
                        call load_G_sinusoids(dgmat,dlos,periodicfreqs,timemast(mm), &
                            & timeserv(mm),kk,nonul,numvars,numperiods)
                     endif

                     if (dodemcor.eq.1) then
                         if (incangle(mm,jj).gt.10.d0.and.incangle(mm,jj).lt.80.d0.and.plthgt(mm).gt.0.and.&
                                &incangle(mm,jj).ne.nulinc(mm)) then
                            cosinc = cos(d2r*incangle(mm,jj))
                            sininc = sin(d2r*incangle(mm,jj))
                            rho0 = dsqrt((plthgt(mm) + radearth)**2 - (sininc*radearth)**2) &
                                & - radearth*cosinc
                            dgmat(kk,numvars) = -1.e3*bperp(mm)/(rho0*sininc) !!! scale for G condition
                         else
                            dgmat(kk,numvars) = 0.d0  !!! no sensitivity to dem or not enough info
                         endif
                     endif
                     kk = kk + 1

                     if(azlogic(mm).eq.1.and.azoff(mm,jj).ne.nulaz(mm).and.hdgval(mm,jj).ne.nulhdg(mm)) then
                        dvmat(kk,kk) = ccoef(mm)*dcor/(1-dcor)
                        dlos(1) = sin(dble(hdgval(mm,jj))*d2r)
                        dlos(2) = cos(dble(hdgval(mm,jj))*d2r)
                        dlos(3) = 0.d0

                        dgmat(kk,1) = dlos(1)*deltime
                        dgmat(kk,2) = dlos(2)*deltime
                        dgmat(kk,3) = dlos(3)*deltime
                        ddvec(kk) = dble(azoff(mm,jj)*aziconv(mm))
                        dscor = dscor + dble(cor(mm,jj))

                        if (numperiods.gt.0) then
                           call load_G_sinusoids(dgmat,dlos,periodicfreqs,timemast(mm), &
                               & timeserv(mm),kk,nonul,numvars,numperiods)
                        endif
                        if (dodemcor.eq.1) then
                            dgmat(kk,numvars) = 0.  !!! azimuth offsets are pure horizontal
                        endif
                        kk = kk + 1
                     endif
                  endif
               enddo

               !!! now that we have G, Cd, and d, it's on
               dtempmat = matmul(transpose(dgmat),dvmat)
               dgtg = matmul(dtempmat,dgmat) + drtr

               call getrank(dgtg,numvars,numvars,rnk,rnktol,rcond,info)
               !!! gives the option to not care about rank results...results will likely be noisy
               !!!  around the edges but there should be no difference in well observed areas
               if (rnkcheck.ne.1) then
                  rnk = numvars
               endif

               if (rnk.ge.numvars) then
                  dgtgcopy = dgtg
                  dxvec = matmul(dtempmat,ddvec)

                  !!! LAPACk QR pivoting (fairly robust and fast) least squares solver
                  call dgelsy(numvars,numvars,1,dgtgcopy,numvars,dxvec,numvars,jpvt,rcond, &
                                 & rnkdgelsy,work,lwork,info)

                  if (info.eq.0) then
                     east(jj)  = sngl(dxvec(1))
                     north(jj) = sngl(dxvec(2))
                     up(jj)    = sngl(dxvec(3))

                     do pp=1,numperiods
                        cose = dxvec(4+6*(pp-1))
                        cosn = dxvec(5+6*(pp-1))
                        cosu = dxvec(6+6*(pp-1))
                        sine = dxvec(7+6*(pp-1))
                        sinn = dxvec(8+6*(pp-1))
                        sinu = dxvec(9+6*(pp-1))

                        phze(pp,jj) = sngl(modulo(datan2(cose,sine),twopi))
                        phzn(pp,jj) = sngl(modulo(datan2(cosn,sinn),twopi))
                        phzu(pp,jj) = sngl(modulo(datan2(cosu,sinu),twopi))
                        ampe(pp,jj) = sngl(dsqrt(cose**2 + sine**2))
                        ampn(pp,jj) = sngl(dsqrt(cosn**2 + sinn**2))
                        ampu(pp,jj) = sngl(dsqrt(cosu**2 + sinu**2))
                     enddo

                     if (dodemcor.eq.1) then
                        dem(jj) = 1.e3*dxvec(numvars)
                     endif

                     if (corout.eq.1) then
                        avecor(jj) = sngl(dscor/dble(nonul))
                     endif
                     if (vm.eq.1) then
                        vmag(jj) = sngl(dsqrt(dxvec(1)**2+dxvec(2)**2+dxvec(3)**2))
                     endif
                     if (gtgout.eq.1.or.gdopest.eq.1.or.gtgoff.eq.1) then
                        gdopmat = matmul(transpose(dgmat),dgmat)  !!! might want to add drtr...think about it
                        call dinv(numvars,gdopmat)
                        if (gtgout.eq.1) then
                           gtgeast(jj)  = sngl(gdopmat(1,1))
                           gtgnorth(jj) = sngl(gdopmat(2,2))
                           gtgup(jj)    = sngl(gdopmat(3,3))
                           do pp=1,numperiods
                               gtgsin(1+6*(pp-1),jj) = sngl(amp_error(gdopmat(4+6*(pp-1)+3,4+6*(pp-1)+3),&
                                   & gdopmat(1+6*(pp-1)+3,1+6*(pp-1)+3),dble(phze(pp,jj))))
                               gtgsin(2+6*(pp-1),jj) = sngl(amp_error(gdopmat(5+6*(pp-1)+3,5+6*(pp-1)+3),&
                                   & gdopmat(2+6*(pp-1)+3,2+6*(pp-1)+3),dble(phzn(pp,jj))))
                               gtgsin(3+6*(pp-1),jj) = sngl(amp_error(gdopmat(6+6*(pp-1)+3,6+6*(pp-1)+3),&
                                   & gdopmat(3+6*(pp-1)+3,3+6*(pp-1)+3),dble(phzu(pp,jj))))
                               gtgsin(4+6*(pp-1),jj) = sngl(phase_error(gdopmat(4+6*(pp-1)+3,4+6*(pp-1)+3),&
                                   & gdopmat(1+6*(pp-1)+3,1+6*(pp-1)+3),dble(phze(pp,jj)),dble(ampe(pp,jj))))
                               gtgsin(5+6*(pp-1),jj) = sngl(phase_error(gdopmat(5+6*(pp-1)+3,5+6*(pp-1)+3),&
                                   & gdopmat(2+6*(pp-1)+3,2+6*(pp-1)+3),dble(phzn(pp,jj)),dble(ampn(pp,jj))))
                               gtgsin(6+6*(pp-1),jj) = sngl(phase_error(gdopmat(6+6*(pp-1)+3,6+6*(pp-1)+3),&
                                   & gdopmat(3+6*(pp-1)+3,3+6*(pp-1)+3),dble(phzu(pp,jj)),dble(ampu(pp,jj))))
                           enddo
                           if (dodemcor.eq.1) then
                              gtgdem(jj) = sngl(gdopmat(numvars,numvars))
                           endif
                        endif
                        if (gtgoff.eq.1) then
                           gtgoen(jj) = sngl(gdopmat(1,2))
                           gtgonu(jj) = sngl(gdopmat(2,3))
                           gtgoeu(jj) = sngl(gdopmat(1,3))
                        endif
                        if (gdopest.eq.1) then
                           sumerr = 0.d0
                           do pp=1,numvars
                              sumerr = sumerr + gdopmat(pp,pp)
                           enddo
                           gdopvec(jj) = sngl(dsqrt(sumerr))
                        endif
                     endif
                     if (errest.eq.1.or.ocoff.eq.1.or.getwdp.eq.1) then
                        call dinv(numvars,dgtg)
                        if (errest.eq.1) then
                           erre(jj) = sngl(dgtg(1,1))
                           errn(jj) = sngl(dgtg(2,2))
                           erru(jj) = sngl(dgtg(3,3))
                           do pp=1,numperiods
                              errsin(1+6*(pp-1),jj) = sngl(amp_error(dgtg(4+6*(pp-1)+3,4+6*(pp-1)+3),&
                                  & dgtg(1+6*(pp-1)+3,1+6*(pp-1)+3),dble(phze(pp,jj))))
                              errsin(2+6*(pp-1),jj) = sngl(amp_error(dgtg(5+6*(pp-1)+3,5+6*(pp-1)+3),&
                                  & dgtg(2+6*(pp-1)+3,2+6*(pp-1)+3),dble(phzn(pp,jj))))
                              errsin(3+6*(pp-1),jj) = sngl(amp_error(dgtg(6+6*(pp-1)+3,6+6*(pp-1)+3),&
                                  & dgtg(3+6*(pp-1)+3,3+6*(pp-1)+3),dble(phzu(pp,jj))))
                              errsin(4+6*(pp-1),jj) = sngl(phase_error(dgtg(4+6*(pp-1)+3,4+6*(pp-1)+3),&
                                  & dgtg(1+6*(pp-1)+3,1+6*(pp-1)+3),dble(phze(pp,jj)),dble(ampe(pp,jj))))
                              errsin(5+6*(pp-1),jj) = sngl(phase_error(dgtg(5+6*(pp-1)+3,5+6*(pp-1)+3),&
                                  & dgtg(2+6*(pp-1)+3,2+6*(pp-1)+3),dble(phzn(pp,jj)),dble(ampn(pp,jj))))
                              errsin(6+6*(pp-1),jj) = sngl(phase_error(dgtg(6+6*(pp-1)+3,6+6*(pp-1)+3),&
                                  & dgtg(3+6*(pp-1)+3,3+6*(pp-1)+3),dble(phzu(pp,jj)),dble(ampu(pp,jj))))
                           enddo
                           if (dodemcor.eq.1) then
                              errdem(jj) = sngl(dgtg(numvars,numvars))
                           endif
                        endif
                        if (ocoff.eq.1) then
                           ocoen(jj) = sngl(dgtg(1,2))
                           oconu(jj) = sngl(dgtg(2,3))
                           ocoeu(jj) = sngl(dgtg(1,3))
                        endif
                        if (getwdp.eq.1) then
                           sumerr = 0.d0
                           do pp=1,numvars
                              sumerr = sumerr + dgtg(pp,pp)
                           enddo
                           wdop(jj) = sngl(dsqrt(sumerr))
                        endif
                     endif
                     if (mseout.eq.1) then
                        if (nonul.gt.2) then
                           mse = 0.d0 ! depreciated: was: sum(dyvec**2)/dble(nonul) where dyvec came from dggglm
                        else
                           mse = 0.d0
                        endif
                        msev(jj) = sngl(mse)
                     endif
                  endif
               endif
               deallocate(dgmat,dvmat,ddvec,dtempmat)
            endif
         enddo
         !$omp end parallel do
      endif

      write(11,rec=ii) (east(mm),mm=1,colsdom)
      write(12,rec=ii) (north(mm),mm=1,colsdom)
      write(13,rec=ii) (up(mm),mm=1,colsdom)
      write(19,rec=ii) (numout(mm),mm=1,colsdom)

      if (errest.eq.1) then
         write(14,rec=ii) (erre(mm),mm=1,colsdom)
         write(15,rec=ii) (errn(mm),mm=1,colsdom)
         write(16,rec=ii) (erru(mm),mm=1,colsdom)
      endif
      if (mseout.eq.1) write(17,rec=ii) (msev(mm),mm=1,colsdom)
      if (vm.eq.1) write(18,rec=ii) (vmag(mm),mm=1,colsdom)
      if (corout.eq.1) write(20,rec=ii) (avecor(mm),mm=1,colsdom)
      if (gdopest.eq.1) write(21,rec=ii) (gdopvec(mm),mm=1,colsdom)
      if (gtgout.eq.1) then
         write(22,rec=ii) (gtgeast(mm),mm=1,colsdom)
         write(23,rec=ii) (gtgnorth(mm),mm=1,colsdom)
         write(24,rec=ii) (gtgup(mm),mm=1,colsdom)
         iid = periodgtgiid
         do pp=1,numperiods
            iid = iid + 1  !!! do propogation of errors above to output the proper errors
         enddo
      endif
      if (gtgoff.eq.1) then
         write(25,rec=ii) (gtgoen(mm),mm=1,colsdom)
         write(26,rec=ii) (gtgonu(mm),mm=1,colsdom)
         write(27,rec=ii) (gtgoeu(mm),mm=1,colsdom)
      endif
      if (ocoff.eq.1) then
         write(28,rec=ii) (ocoen(mm),mm=1,colsdom)
         write(29,rec=ii) (oconu(mm),mm=1,colsdom)
         write(30,rec=ii) (ocoeu(mm),mm=1,colsdom)
      endif
      if (getwdp.eq.1) then
         write(31,rec=ii) (wdop(mm),mm=1,colsdom)
      endif
      if (dodemcor.eq.1) then
         write(32,rec=ii) (dem(mm),mm=1,colsdom)
      endif
      iid = periodiid
      do pp=1,numperiods
         write(iid,rec=ii) (phze(pp,mm),mm=1,colsdom)
         iid = iid + 1
         write(iid,rec=ii) (phzn(pp,mm),mm=1,colsdom)
         iid = iid + 1
         write(iid,rec=ii) (phzu(pp,mm),mm=1,colsdom)
         iid = iid + 1
         write(iid,rec=ii) (ampe(pp,mm),mm=1,colsdom)
         iid = iid + 1
         write(iid,rec=ii) (ampn(pp,mm),mm=1,colsdom)
         iid = iid + 1
         write(iid,rec=ii) (ampu(pp,mm),mm=1,colsdom)
         iid = iid + 1
      enddo
      if (gtgout.eq.1.or.gdopest.eq.1.or.gtgoff.eq.1) then
         iid = periodgtgiid
         do pp=1,numperiods
            write(iid,rec=ii) (gtgsin(1+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (gtgsin(2+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (gtgsin(3+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (gtgsin(4+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (gtgsin(5+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (gtgsin(6+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
         enddo
         if (dodemcor.eq.1) then
            write(iid,rec=ii) (gtgdem(mm),mm=1,colsdom)
            iid = iid + 1
         endif
      endif
      if (errest.eq.1.or.ocoff.eq.1.or.getwdp.eq.1) then
         iid = perioderriid
         do pp=1,numperiods
            write(iid,rec=ii) (errsin(1+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (errsin(2+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (errsin(3+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (errsin(4+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (errsin(5+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
            write(iid,rec=ii) (errsin(6+6*(pp-1),mm),mm=1,colsdom)
            iid = iid + 1
         enddo
         if (dodemcor.eq.1) then
            write(iid,rec=ii) (errdem(mm),mm=1,colsdom)
            iid = iid + 1
         endif
      endif
   enddo

   if (use_mmap.eq.1) then !!! unmap mmap'ed files
      do ii=1,numscenes
         c_munmapret = munmap(c_loc(unwptr(ii)%p(1)),c_len)
         if (c_munmapret.ne.0) print *,'could not munmap unwptr'
         c_munmapret = munmap(c_loc(corptr(ii)%p(1)),c_len)
         if (c_munmapret.ne.0) print *,'could not munmap corptr'
         c_munmapret = munmap(c_loc(losptr(ii)%p(1)),c_len_los)
         if (c_munmapret.ne.0) print *,'could not munmap losptr'

         if (azlogic(ii).eq.1) then
            c_munmapret = munmap(c_loc(aziptr(ii)%p(1)),c_len)
            if (c_munmapret.ne.0) print *,'could not munmap aziptr'

            c_munmapret = munmap(c_loc(hdgptr(ii)%p(1)),c_len)
            if (c_munmapret.ne.0) print *,'could not munmap hdgptr'
         endif

         if (dodemcor.eq.1) then
            c_munmapret = munmap(c_loc(incptr(ii)%p(1)),c_len)
            if (c_munmapret.ne.0) print *,'could not munmap incptr'
         endif
      enddo
   endif

   print*,' ';print*,' '
   print*,'fourDvel is done'
   print*,' '

100   format(a8,i6)
110   format(a15,f9.2)
123   format(a,f16.10,a,f16.10,a,f14.10,a,f14.10,a)

   end program
!  *************************************************************************************
!  *************************************************************************************
   function amp_error(sinerr,coserr,phase)
   !!! expects sinerr and coserr to be straight from inv(C_m) (i.e. sigma**2)
   !!! returns sigma_amp**2
   real*8   amp_error
   real*8   sinerr,coserr,phase,ssphz,csphz

   ssphz = sin(phase)**2
   csphz = cos(phase)**2
   amp_error = (coserr*ssphz - sinerr*csphz)/(ssphz**2 - csphz**2)
   return
   end function

!  *************************************************************************************
!  *************************************************************************************
   function phase_error(sinerr,coserr,phase,amp)
   real*8   phase_error
   real*8   sinerr,coserr,phase,amp,ssphz,csphz,ampsq

   ssphz = sin(phase)**2
   csphz = cos(phase)**2
   ampsq = amp**2
   phase_error = (sinerr*ssphz - coserr*csphz)/(ampsq*(ssphz**2 - csphz**2))
   return
   end function

!  *************************************************************************************
!  *************************************************************************************
   subroutine load_G_sinusoids(dgmat,dlos,perfreq,tm,ts,kk,nonul,numv,nump)
   implicit none
   integer pp,kk,nump,nonul,numv
   real*8  delcos,delsin,dgmat(nonul,numv),dlos(3)
   real*8  perfreq(nump),tm,ts

   do pp=1,nump
       delcos = cos(perfreq(pp)*tm) - cos(perfreq(pp)*ts)
       delsin = sin(perfreq(pp)*tm) - sin(perfreq(pp)*ts)
       dgmat(kk,4+6*(pp-1)) = delcos*dlos(1)
       dgmat(kk,5+6*(pp-1)) = delcos*dlos(2)
       dgmat(kk,6+6*(pp-1)) = delcos*dlos(3)
       dgmat(kk,7+6*(pp-1)) = delsin*dlos(1)
       dgmat(kk,8+6*(pp-1)) = delsin*dlos(2)
       dgmat(kk,9+6*(pp-1)) = delsin*dlos(3)
   enddo
   end subroutine

!  *************************************************************************************
!  *************************************************************************************
   subroutine load_Cmprior(drtr,alpha,perfreq,numv,nump,reffreqhor,reffreqvrt)
   implicit none
   integer pp,numv,nump
   real*8  drtr(numv,numv),perfreq(nump),horfreq,vrtfreq,alpha,reffreqhor,reffreqvrt

   drtr = 0.d0
   if (reffreqhor.gt.0.d0) then
      horfreq = reffreqhor
   else
      horfreq = minval(perfreq)
   endif

   if (reffreqvrt.gt.0.d0) then
      vrtfreq = reffreqvrt
   else
      vrtfreq = maxval(perfreq)
   endif

   do pp=1,nump
      drtr(4+6*(pp-1),4+6*(pp-1)) = alpha*(perfreq(pp)/horfreq - 1.d0)**2
      drtr(5+6*(pp-1),5+6*(pp-1)) = alpha*(perfreq(pp)/horfreq - 1.d0)**2
      drtr(6+6*(pp-1),6+6*(pp-1)) = alpha*(vrtfreq/perfreq(pp) - 1.d0)**2
      drtr(7+6*(pp-1),7+6*(pp-1)) = alpha*(perfreq(pp)/horfreq - 1.d0)**2
      drtr(8+6*(pp-1),8+6*(pp-1)) = alpha*(perfreq(pp)/horfreq - 1.d0)**2
      drtr(9+6*(pp-1),9+6*(pp-1)) = alpha*(vrtfreq/perfreq(pp) - 1.d0)**2
   enddo
   !!! DEM term if applicable
   if ((3+6*nump).lt.numv) then
      !drtr(numv,1:numv-1) = alpha
      drtr(numv,numv) = alpha
   endif
   end subroutine

!  *************************************************************************************
!  *************************************************************************************
   subroutine getrank(a,n,m,rnk,st,rcond,info)
!  returns rank(A) for a double precision nXm matrix A
   implicit none
   integer n,m,rnk,info,i,lwork
   real*8  a(n,m),pa(n,m),s(m),u(n,m),vt(n,m),st,thresh,eps,rcond
   real*8,allocatable:: work(:)
   lwork = 10*(m+n)
   allocate(work(lwork))

   pa = a
   call dgesvd('N','N',n,m,pa,n,s,u,n,vt,n,work,lwork,info)
   if (info.gt.0) then
      !print *,'dgesvd did not converge by number of superdiagonals = ',info
      rnk = -1
   elseif (info.lt.0) then
      print *,'dgesvd: illegal entry in position ',abs(info)
      stop
   else
      !!! rcond = reciprocal condition number ~ min(eig)/max(eig) [inverse of numpy.linalg.cond]
      rcond = s(m)/s(1)
      s = s/s(1)
      if (st.gt.0) then
           thresh = st
      else
           eps = epsilon(st)
           thresh = n*eps
      endif
      rnk = 0
      do i=1,m
         if (s(i).gt.thresh) rnk = rnk + 1
      enddo
   endif
   end subroutine

!  *************************************************************************************
!  *************************************************************************************
   subroutine get_gtg_rcond(dgtg,n,rcond,info)
   implicit none
   integer n,m,info,lwork
   real*8  dgtg(n,n),dgtgcopy(n,n),rcond,s(n),u(n,n),vt(n,n)
   real*8,allocatable:: work(:)
   lwork = 20*n
   allocate(work(lwork))

   !!! rcond = reciprocal condition number ~ min(eig)/max(eig) [inverse of numpy.linalg.cond]
   dgtgcopy = dgtg
   call dgesvd('N','N',n,n,dgtgcopy,n,s,u,n,vt,n,work,lwork,info)
   rcond = s(n)/s(1)
   end subroutine

!  *************************************************************************************
!  *************************************************************************************
   subroutine dinv(n,mat)
!  returns inv(mat) for a double precision nXn matrix mat
   implicit none
   integer     n,info,ipiv(n)
   real*8      mat(n,n), work(n*n)
   call dgetrf(n,n,mat,n,ipiv,info)
   call dgetri(n,mat,n,ipiv,work,n*n,info)
   end subroutine

!  *************************************************************************************
!  *************************************************************************************
   subroutine getmse(v,b,g,gtg,m,n,mse)
!  calculates mean square error scalar mse
!
!  inputs:   v = inverse covariance matrix
!            b = data vector
!            g = linear operator
!            gtg = inv(transpose(g)*inv(v)*g)
!            m = number of observations
!            n = size of model (will be 3 for 3D vector field)

   implicit none

   integer  m,n
   real*4   gtg(n,n),v(m,m),g(m,n),b(m)
   real*4   va(m,n),vagtg(m,n),at(n,m),atv(n,m)
   real*4   big(m,m),par1(m),mse

   if (m.eq.n) then
      mse = 0
      return
   endif

   at = transpose(g)
   va = matmul(v,g)
   vagtg = matmul(va,gtg)
   atv = matmul(at,v)
   big = v - matmul(vagtg,atv)

   par1(1) = b(1)*big(1,1) + b(2)*big(2,1) + b(3)*big(3,1)
   par1(2) = b(1)*big(1,2) + b(2)*big(2,2) + b(3)*big(3,2)
   par1(3) = b(1)*big(1,3) + b(2)*big(2,3) + b(3)*big(3,3)

   mse  = (par1(1)*b(1) + par1(2)*b(2) + par1(3)*b(3))/(m-n)

   end subroutine

!  *************************************************************************************
!  *************************************************************************************
    subroutine get_periodics(str,pers,m,n)
        character*500    :: str
        integer          :: m,n,i,j,ios,pos2,pos1
        real*8           :: pers(m)

        pos1 = 1
        do n=1,m
            pos2 = index(str(pos1:),",")
            if (pos2.eq.0) then
                read(str(pos1:),*,iostat=ios) pers(n)
                exit
            endif
            read(str(pos1:pos1+pos2-2),*,iostat=ios) pers(n)
            pos1 = pos2 + pos1
        enddo
    end subroutine

!  *************************************************************************************
!  *************************************************************************************
    subroutine parse_bbox(str,bbox,m,n)
        character*500    :: str
        integer          :: n,m,ios,pos2,pos1
        real*8           :: bbox(4)

        pos1 = 1
        do n=1,m
            pos2 = index(str(pos1:),",")
            if (pos2.eq.0) then
                read(str(pos1:),*,iostat=ios) bbox(n)
                exit
            endif
            read(str(pos1:pos1+pos2-2),*,iostat=ios) bbox(n)
            pos1 = pos2 + pos1
        enddo
    end subroutine

!  *************************************************************************************
!  *************************************************************************************
    subroutine bsort_int(a,n)
        integer :: n,a(n),temp,i,j
        logical :: swapped

        do j = n-1, 1, -1
        swapped = .false.
        do i = 1, j
          if (a(i) > a(i+1)) then
            temp = a(i)
            a(i) = a(i+1)
            a(i+1) = temp
            swapped = .true.
          end if
        end do
        if (.not. swapped) exit
        end do
    end subroutine

!  *************************************************************************************
!  *************************************************************************************

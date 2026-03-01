! SBBX_to_nc.f
!
! This program uses 2 input files:
!  the first contains        SURFACE AIR TEMPERATURE anomalies,
!  the second          OCEAN MIXED LAYER TEMPERATURE anomalies.
!
! It reads monthly anomalies for 8000 equal area subboxes, finds the  anomaly
! with respect to 1951-1980 for a user-specified month/year.
!
! The results are replicated as a gridded netCDF data file. The anomaly variable
! is stored as a scaled int -- i.e., intval = nint(realval*100) -- in order to
! reduce the size of the output dataset.
!
! This program was originally written c. 2005 and run on Mac OS X and AIX5.
! More recently (2013), it was rewritten and compiled to run on Linux CentOS 6.
! On the Linux machine, C libraries and utilities were installed via yum, e.g.,
!
!   yum install netcdf
!   yum install netcdf-devel
!
! This program was then compiled and linked using the command
!
!   gfortran -fconvert=big-endian SBBX_to_nc.f -o SBBX_to_nc.exe \
!            -L/usr/lib64 -lnetcdf -lnetcdff -I/usr/include
!
! Although this program is written in FORTRAN-90, calls to netCDF library
! routines use the netCDF FORTRAN 77 API. See
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77/
!
! This version of the program expects six command-line arguments when the
! program is executed. More or less will cause the program to stop.
!
!    iland   = 0    -> do not include land surface data;
!            = 250  -> use TS250_DATA for input land surface data;
!            = 1200 -> use TS1200_DATA for input land surface data
!    iocean  = 0 -> do not include ocean data;
!            = 1 -> use HadR2_DATA for input ocean data;
!            = 3 -> use ERSST_DATA for input ocean data
!            = 4 -> use ERSSTv4_DATA for input ocean data
!            = 5 -> use ERSSTv5_DATA for input ocean data
!    ifyr    = first year to include in output file
!    ilyr    = last year to include in output file
!    ilmo    = last month in last year to include (1..12)
!    outfile = name of the output netCDF file (e.g., gistemp.nc)
!
! For example,
!
!    ./SBBX_to_nc.exe 1200 4 1880 2013 7 gistemp1200.nc
!
! Creates a file with data from Jan 1880 through Jul 2013 using TS1200_DATA
! and ERSSTv4_DATA as input.
!
! You can alter parameters IM and JM to achieve different gridding, although
! 2-degree gridding is optimal.
!
! Both land and ocean input files have the same structure:
! Record 1 starts with 8 integers I1-I8 and an 80-byte TITLE.
! All further records start with 7 integers N1-N7,
!             a real number R, followed by a data array (real*4).
! I1 or N1 is the length of the data array in the NEXT record.
! Unless its length is 0, each data-array contains a time series
! of monthly T-anomalies (C) starting with January of year I6 for
! one grid box. N2,N3,N4,N5 indicate the edges in .01 degrees of
! that grid box in the order: latitude of southern edge, latitude of
! northern edge, longitude of western edge, longit. of eastern edge.
! The world is covered with 8000 equal area grid boxes, so each
! file has 8001 records. I7 is the flag for miss data (9999).
!
! Please note that depending on your choices, the output file can be LARGE.
! An output file including data from 1880 to 2013 is about 52 MB.
!
      program SBBX_to_nc
!     implicit none

      include "netcdf.inc"

    ! Output resolution. Default is 2-deg.
      integer, parameter :: im = 180   ! lon. grid points, 180 = 2-deg resolution
      integer, parameter :: jm =  90   ! lat. grid points,  90 = 2-deg resolution

    ! Arrays for coordinate variable values
      real, allocatable :: lat(:), lon(:)
      integer, allocatable :: time(:), timeb(:)

    ! The temperature anomaly output array
      integer*2, allocatable :: tout(:,:,:)

    ! Input anomaly arrays.
      real, allocatable :: tin(:), tav(:)
      real, allocatable :: tino(:), tavo(:)
      real w_eag(8000)
      real, allocatable :: t_gcm(:,:)
      real, allocatable :: t_eag(:,:)

      integer istatus            ! Placeholder for status returned by all nc_* functions
      integer ncid               ! ID number for NC dataset, as assigned by nc_create
      integer latdid, londid     ! IDs for the dimensions, as assigned by nc_def_dim
      integer timedid, timebdid
      integer latvid, lonvid     ! IDs for the variables, as assigned by nc_def_var
      integer timevid, timebvid
      integer datavid
      integer bshape(2)          ! Shape of time boundary coordinate array
      integer dshape(3)          ! Shape of anomaly data array
      integer jday1800           ! All dates are given as days offset from Jan 1, 1800
      integer*2 att_missing(1)   ! Array to hold missing-data/fill value

    ! Some constants for the output file global metadata.
      character(*), parameter :: nctitle
     *              = "GISTEMP Surface Temperature Analysis"
      character(*), parameter :: ncinst
     *              = "NASA Goddard Institute for Space Studies"
      character(*), parameter :: ncsource
     *              = "http://data.giss.nasa.gov/gistemp/"
      character(*), parameter :: ncconv = "CF-1.6"

    ! And more info for the output file global metadata.
      character(8)  :: rundate
      character(10) :: runtime
      character(92) :: nchist

    ! User inputs and vars for handling command-line params
      character(80) :: buffer
      integer :: numargs
      integer :: iland, iocean
      integer :: ifyr, ilyr, ilmo
      integer :: ifyrb, ilyrb    ! First and last year of base period
      character(80) :: outfile

    ! ...
      character*80 title, titleo
      integer info(8), infoo(8)
      real :: bad, dl, dlo

      character(len=8), dimension(0:5) :: legendOcn = (/'no ocean',
     *                            'HadReyv2','        ','ERSSTv3 ',
     *                            'ERSSTv4 ','ERSSTv5 '/)

      integer :: ifyrm1         ! first year of output minus 1
      integer :: iyears         ! number of years (including partials) in output
      integer :: imonths        ! number of months in output
      integer :: navgb          ! number of years in base period
      real :: rland
      real :: flatres, flonres  ! lat and lon resolution of data
      real :: fsouth, fwest     ! first lat and lon array values
      integer :: imm, immb, imoff, imoffb, it, iendmo
      integer :: iyr, iyrb, iyrbgc, iyrbeg, iyrbgo, iyy
      integer :: i1tin, i1tino, monmc, iyrend
      integer :: i, j           ! do loop indexes

!
! The inputs
!
      iland    = 1200       ! 0 = no land, 250/1200 = radius
      iocean   =    5       ! 0 = no ocean, 1 = H/Rv2, 3 = NCDC/ERv3
                            ! 4 = NCDC/ERv4, 5 = NCDC/ERv5
      ifyr     = 1880       ! first year of output
      ilyr     = 2013       ! last year of output
      ilmo     =    7       ! last month of output
      ifyrb    = 1951       ! first year of base period
      ilyrb    = 1980       ! last year of base period

      numargs = iargc()

      if (numargs == 0) then
        print *, 'Need: iland,iocean, yr1,yr2,mo2, outfile, yrfb,yrlb'
        stop '...Stopped...'
      end if

      call getarg (1,buffer) ; read (buffer,*) iland
      if (numargs>1) then
         call getarg (2,buffer) ; read (buffer,*) iocean
      end if
      if (iocean > 0) outfile = trim(legendOcn(iocean))//'.nc'
      if (iland == 250)  outfile='gistemp250.nc'
      if (iland == 1200) outfile='gistemp1200.nc'
      if (numargs > 2) then
         call getarg (3,buffer) ; read (buffer,*) ifyr
      end if
      if (numargs > 3) then
         call getarg (4,buffer) ; read (buffer,*) ilyr
      end if
      if (numargs > 4) then
         call getarg (5,buffer) ; read (buffer,*) ilmo
      end if
      if (numargs > 5) call getarg (6,outfile)
      if (numargs > 6) then
         call getarg (7,buffer) ; read (buffer,*) ifyrb
      end if
      if (numargs > 7) then
         call getarg (8,buffer) ; read (buffer,*) ilyrb
      end if

      ifyrm1  = ifyr - 1
      iyears  = ilyr - ifyr + 1
      imonths = 12 * (iyears - 1) + ilmo

      if(iland==0) write (6,*) 'ILAND   = ', iland, 'no land data used'
      if(iland >0) write (6,*) 'Smoothing radius (km) = ', iland
      write (6,*) 'IOCEAN  = ', iocean,legendOcn(iocean)
      write (6,*) 'IFYR    = ', ifyr
      write (6,*) 'ILYR    = ', ilyr
      write (6,*) 'ILMO    = ', ilmo
      write (6,*) 'OUTFILE = ', outfile
      write (6,*) 'Base period  = ', ifyrb,'-',ilyrb

!
! Start looking at the input.
!
      if (iland == 0 .and. iocean == 0) then
        print *, 'Must select land, ocean or both'
        stop '...Stopped...'
      endif

                       rland =   100.  ! land&ocen is used, station reach 100 km
      if (iland  == 0) rland = -9999.  ! only ocean data are used
      if (iocean == 0) rland =  9999.  ! only land data are used

!
! Initialize the dataset, define global metadata attributes.
!
      print *, 'Creating dataset...'
      istatus = nf_CREATE (outfile, nf_CLOBBER, ncid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      print *, 'Writing global attributes...'
      istatus = nf_put_att_text (
     *     ncid, nf_global, 'title', len(nctitle), nctitle)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text (
     *     ncid, nf_global, 'institution', len(ncinst), ncinst)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text (
     *     ncid, nf_global, 'source', len(ncsource), ncsource)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text (
     *     ncid, nf_global, 'Conventions', len(ncconv), ncconv)
      if (istatus /= nf_noerr) call handle_err(istatus)

      call date_and_time(DATE=rundate, TIME=runtime)
      nchist = 'Created '
     *     // rundate(1:4) // '-' // rundate(5:6) // '-' // rundate(7:8)
     *     // ' '
     *     // runtime(1:2) // ':' // runtime(3:4) // ':' // runtime(5:6)
     *     // ' by SBBX_to_nc 2.0 - '
     *     // 'ILAND=****  '
     *     // 'IOCEAN=***       '
     *     // 'Base: YYYY-YYYY'

      if (iland == 0) then
         nchist(55:59) = 'none,'
      else if (iland == 250) then
         nchist(55:59) = '250, '
      else if (iland == 1200) then
         nchist(55:59) = '1200,'
      else
         print *, 'iland must be 0, 250 or 1200'
         stop '...Stopped...'
      end if
      if (iocean == 0) then
         nchist(68:76) = 'none,    '
      else if (iocean == 1) then
         nchist(68:76) = 'Had/R2,  '
      else if (iocean == 3) then
         nchist(68:76) = 'NCDC/ER3,'
      else if (iocean == 4) then
         nchist(68:76) = 'NCDC/ER4,'
      else if (iocean == 5) then
         nchist(68:76) = 'NCDC/ER5,'
      else
         print *, 'iocean must be 0, 1 or 3, 4, 5'
         stop '...Stopped...'
      end if
      write (nchist(84:92), '(i4,a1,i4)') ifyrb, '-', ilyrb

      istatus = nf_put_att_text (
     *     ncid, nf_global, 'history', len(nchist), nchist)
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Define the dimensions.
!
      print *, 'Defining dimensions...'
      istatus = nf_def_DIM(ncid, 'lat', jm, latdid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_def_DIM(ncid, 'lon', im, londid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_def_DIM(ncid, 'time', imonths, timedid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_def_DIM(ncid, 'nv', 2, timebdid)
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Define the coordinate variable values.
! Dates are measured as days since 1800-01-01.
!
      allocate ( lat(jm) )
      allocate ( lon(im) )
      allocate ( time(imonths) )
      allocate ( timeb(imonths*2) )

      flatres = 180. / jm
      fsouth =  -90. + 0.5 * flatres

      flonres = 360. / im
      fwest  = -180. + 0.5 * flonres
      offlon=0. ! grid starts at the date line

      do J = 1, jm
        lat(J) = fsouth + (J - 1.) * flatres
      end do

      do I = 1, im
        lon(I) = fwest  + (I - 1.) * flonres
      end do

      jday1800 = julday(1800,1,1)
      do iyy = 1, iyears
        iyr = ifyrm1 + iyy
        iendmo = 12
        if (iyr == ilyr) iendmo = ilmo
        do imm = 1, iendmo
          imoff  = (iyy-1)*12 + imm
          time(imoff) = julday(iyr,imm,15) - jday1800

          iyrb = iyr
          immb = imm + 1
          if (immb == 13) then
            iyrb = iyr+1
            immb = 1
          end if
          imoffb = imoff*2 - 1
          timeb(imoffb)   = julday(iyr, imm, 1) - jday1800
          timeb(imoffb+1) = julday(iyrb,immb,1) - jday1800
        end do
      end do

      print *, 'Defining coordinate variables...'

!
! Define latitude coordinate variable.
!
      istatus = nf_def_var(ncid,'lat',nf_float,1,latdid,latvid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,latvid,'standard_name',8,'latitude')
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,latvid,'long_name',8,'Latitude')
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,latvid,'units',13,'degrees_north')
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Define longitude coordinate variable.
!
      istatus = nf_def_var(ncid,'lon',nf_float,1,londid,lonvid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,lonvid,'standard_name',9,'longitude')
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,lonvid,'long_name',9,'Longitude')
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,lonvid,'units',12,'degrees_east')
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Define time coordinate variable.
!
      istatus = nf_def_var(ncid,'time',nf_int,1,timedid,timevid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,timevid,'long_name',4,'time')
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,timevid,'units',30,'days since 1800-01-01 00:00:00')
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,timevid,'bounds',9,'time_bnds')
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Define time coordinate bounds variable.
!
      bshape(1) = timebdid
      bshape(2) = timedid
      istatus = nf_def_var(ncid,'time_bnds',nf_int,2,bshape,timebvid)
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Define the anomaly variable.
!
      print *, 'Defining the temperature data variable...'
      dshape(1) = londid
      dshape(2) = latdid
      dshape(3) = timedid
      istatus = nf_def_var(ncid,'tempanomaly',nf_short,3,dshape,datavid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,datavid,'long_name',27,'Surface temperature anomaly')
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,datavid,'units',1,'K')
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_real(
     *     ncid,datavid,'scale_factor',NF_FLOAT,1,0.01)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_att_text(
     *     ncid,datavid,'cell_methods',10,'time: mean')
      if (istatus /= nf_noerr) call handle_err(istatus)

      att_missing(1) = 32767
      istatus = nf_put_ATT_int2(
     *     ncid,datavid,'_FillValue',
     *     nf_short,1,att_missing)
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Done defining dataset contents and metadata.
! Commit header information to disk.
!
      istatus = nf_enddef(ncid)
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Write the coordinate variables to the NC dataset
! Deallocate the memory for these vars when done.
!
      print *, 'Writing the coordinate variables to the dataset'

      istatus = nf_put_var_real(ncid, latvid, lat)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_var_real(ncid, lonvid, lon)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_var_int(ncid, timevid, time)
      if (istatus /= nf_noerr) call handle_err(istatus)

      istatus = nf_put_var_int(ncid, timebvid, timeb)
      if (istatus /= nf_noerr) call handle_err(istatus)

      deallocate ( lat )
      deallocate ( lon )
      deallocate ( time )
      deallocate ( timeb )

!
! Examine the input datasets.
! Read the relevant parts of the header records.
!
      if (iland > 0) then
         if (iland == 250) then
            open (8, file='TS250_DATA', form='unformatted',
     *            access='sequential')
         else if (iland == 1200) then
            open (8, file='TS1200_DATA', form='unformatted',
     *            access='sequential')
         else
            print *, 'iland must be 0, 250 or 1200'
            stop '...Stopped...'
         end if
         read (8) info, title
         write (6,*) title
         do i = 1, 8
           infoo(i)=info(i)   !  needed for land-only case
         end do
      end if

      if (iocean > 0) then
         if (iocean == 3) then
            open (9, file='ERSST_DATA', form='unformatted',
     *            access='sequential')
         else if (iocean == 1) then
            open (9, file='HadR2_DATA', form='unformatted',
     *            access='sequential')
         else if (iocean == 4) then
            open (9, file='ERSSTv4_DATA', form='unformatted',
     *            access='sequential')
         else if (iocean == 5) then
            open (9, file='ERSSTv5_DATA', form='unformatted',
     *            access='sequential')
         else
            print *, 'iocean must be 0, 1 or 3, 4'
            stop '...Stopped...'
         end if

         read (9) infoo, titleo
         write (6,*) titleo
         if (iland == 0) then
           do i = 1, 8
             info(i) = infoo(i)
           end do
         end if
      end if

      mnow   = info(1)       ! length of first data record (months)
      monm   = info(4)       ! max length of time series (months)
      iyrbeg = info(6)       ! beginning of time series (calendar year)
      bad    = info(7)       ! bad data value = 9999

      mnowo  = infoo(1)
      monmo  = infoo(4)
      iyrbgo = infoo(6)

!
! Align land and ocean time series
!
      iyrbgc = min(iyrbgo, iyrbeg)            ! use earlier of the 2 start years

      i1tin  = 1 + 12 * (iyrbeg - iyrbgc)     ! land  offset in combined period
      i1tino = 1 + 12 * (iyrbgo - iyrbgc)     ! ocean offset in combined period
      monmc  = max(monm+i1tin-1,monmo+i1tino-1)  ! use later of the 2 ends
      iyrend = iyrbgc - 1 + monmc / 12        ! last calendar year of timeseries

!
! Loop over subboxes - find output data
!
      print *, 'Reading and gridding temperature data...'

!**** Initializations for equal-area grid
      CALL get_ij(0,0,1,ij)
      call eantrp0 (im,jm,offlon,180./flatres,bad)

      allocate ( tin(monmc), tino(monmc), tav(monmc), tavo(monmc) )
      allocate ( t_eag(8000,imonths) ) ! dT on equal-area grid
      t_eag = bad                      ! set all to missing

      do 6100 n = 1, 8000
        tin  = bad     ! set all months to missing initially
        tino = bad

        dl = 9999.          ! in case only ocn data are read in

!**** read in time series TIN/tino of monthly means: land/ocean data

        if (rland.ge.0.) then
          call sread (8, tin(i1tin), mnow, lats, latn, lonw, lone,
     *                dl, next)
          mnow = next ! mnow/next: length of current/next time series
        end if

        wocn = 0.                 ! weight for ocean data
        if (iocean > 0) then   !  read in ocean data
          call sread (9, tino(i1tino), mnowo, lats, latn, lonw, lone,
     *                dlo, nexto)
          mnowo = nexto
          if (dl   > rland) wocn=1. ! dl:subbox_center->nearest station (km)
          if (lats > 7500)  wocn=0. ! disregard SST near sea ice
        end if
!
! At this point the 2 time series TIN,tino are all set and can be used to
! compute means, trends, etc. As an example, we find the requested anomaly:
!
        do 6090 iyy = 1, iyears
          iyr = ifyrm1 + iyy
          iendmo = 12
          if (iyr == ilyr) iendmo = ilmo
        do 6090 imm = 1, iendmo
          imoff = (iyy-1)*12 + imm

! Find the mean over the base period
          tavb    = 0.
          tavbo   = 0.
          tav(1)  = bad
          tavo(1) = bad

!**** m1,m1b: Location in combined series of 1st month needed
          m1  = 12 * (iyr   - iyrbgc) + imm
          m1b = 12 * (ifyrb - iyrbgc) + imm

          navgb = ilyrb+1-ifyrb
          if (rland.ge.0.) then
!*            collect selected month for each base period year
            call avg(tin(m1b), 12, navgb,    1, bad, 1, tav)
!*            find mean over the base period for the selected month
            call avg(tav,   navgb,    1, navgb, bad, 1, tav)
          end if

          if (iocean > 0) then  ! do same for ocean data
            call avg(tino(m1b), 12,navgb,    1, bad, 1, tavO)
            call avg(tavO,  navgb,    1, navgb, bad, 1, tavO)
          end if

          tavb  = tav(1)              ! tavb default: land value
          tavbo = tavo(1)

!
! Put the requested anomaly into tav(1), then t_eag(ij,imoff) (equal-area grid)
!
          tav(1)  = bad
          tavo(1) = bad

          if (tavb /= bad) tav(1) = tin(m1)         ! tav(1): land value

          if (rland < 9999. .and. tavbo /= bad) tavo(1) = tino(m1)

          if (tav(1) == bad .or. rland < 0.) then   ! disregard land data
             tavb   = tavbo                          ! tavb: ocean value
             tav(1) = tavo(1)                        ! tav(1): ocean value
          end if

          if (tavo(1) /= bad) then            ! switch tavb/tav(1) to
              tavb   = tavb   * (1.-wocn) + tavbo   * wocn  ! ocn value if appropriate
              tav(1) = tav(1) * (1.-wocn) + tavo(1) * wocn
          end if
!
! Replicate tav(1) at the appropriate places in the output array
!
          if (tav(1) /= bad) then
             tav(1) = tav(1) - tavb
             call get_ij ((lats+latn)/2,(lonw+lone)/2,0,ij)
             t_eag(ij,imoff) = tav(1)
          end if

 6090   continue

 6100 continue
! End of loop over subboxes

!
! Initialize the output array, filling it with the bad/missing value.
!
      allocate ( tout(im,jm,imonths) )
      tout = 32767

!
! Interpolate to the output regular grid
!
      allocate ( t_gcm(im,jm) )
      do it=1,imonths
        do i=1,8000
          w_eag(i) = 1.
          if(t_eag(i,it) == bad) w_eag(i)=0.
        end do
        call eantrp (w_eag,t_eag(1,it),t_gcm)
        do j=1,jm
          do i=1,im
            if (t_gcm(i,j) == 9999.) then
                tout(i,j,it) = 32767
            else
                tout(i,j,it) = nint (t_gcm(i,j)*100)
            end if
        end do ; end do
      end do

!
! Write output grid to the dataset.
!
      print *, 'Writing the temperature data to the dataset...'
      istatus = nf_put_var_int2(ncid, datavid, tout)
      if (istatus /= nf_noerr) call handle_err(istatus)

!
! Done. Close the NC dataset.
!
      print *, 'Closing the dataset...'
      istatus = nf_close(ncid)
      if (istatus /= nf_noerr) call handle_err(istatus)

      stop
      end

!*********
!*********
!
      function julday(yyyy, mm, dd)
      integer julday, yyyy, mm, dd, greg

      parameter (greg=15+31*(10+12*1582))

      integer jy, jm, ja
      jy = yyyy
      if (jy == 0) then
        print *, 'julday: there is no year zero'
        stop '...Stopped...'
      endif
      if (jy < 0) jy = jy + 1
      if (mm > 2) then
        jm = mm + 1
      else
        jy = jy - 1
        jm = mm + 13
      endif
      julday = int(365.25*jy) + int(30.6001*jm) + dd + 1720995
      if (dd+31*(mm+12*yyyy).ge.greg) then
        ja = int (0.01*jy)
        julday = julday + 2 - ja + int(0.25*ja)
      endif
      return
      end function julday

!*********
!*********
!
      subroutine handle_err(status)
      integer status
      if (status /= nf_noerr) then
        print *, 'ERROR #', status
        print *, NF_STRERROR(status)
        stop '...Stopped...'
      endif
      end subroutine handle_err

!*********
!*********
!
      subroutine sread (ndisk,array,len, N1,N2,N3,N4, dstn,lnext)
      real array(len)
      read (ndisk) lnext, N1,N2,N3,N4, NR1,NR2, dstn, array
      return
      end subroutine sread

!*********
!*********
!
      subroutine avg (array, km, navg, lav, bad, lmin, dav)

      real array(km,navg), dav(navg)

      do 100 n = 1, navg
        sum   = 0.
        kount = 0
        do 50 L = 1, lav
          if (array(L,n) == bad) GO TO 50
          sum   = sum + array(L,n)
          kount = kount + 1
   50   continue
        dav(n) = bad
        if (kount.ge.lmin) dav(n) = sum / kount
  100 continue
      return
      end subroutine avg

!*****
!*****
!
      SUBROUTINE EANTRP0 (INB,JNB,OFFIB,DIVJB, SKIB)
!     IMPLICIT real*4 (A-H, O-Z)
!**** modified by R.Ruedy from a program written by Dr. Gary L. Russell
!****                                                      at NASA/GISS
!**** EANTRP performs a horizontal interpolation of per unit area or per
!**** unit mass quantities defined on grid A, calculating the quantity
!**** on grid B.  B grid values that cannot be calculated because the
!**** covering A grid boxes have WTA = 0, are set to the value of SKIP.
!**** The area weighted integral of the quantity is conserved.
!****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (TWOPI=6.283185307179586477d0)
!     REAL*4  WTA(INA,JNA),A(INA,JNA),B(INB,JNB),
      REAL*4  WTA(8000),      A(8000),      B(*),
     *        OFFIB,DIVJB, SKIB,SKIP
      real*8  SINA(0:80),SINB(0:361),ibefor(81),
     *        FMIN(720,4),FMAX(720,4),GMIN(361),GMAX(361)
      integer IMIN(720,4),IMAX(720,4),JMIN(362),JMAX(362)
      save SINA,SINB,ibefor,FMIN,FMAX,GMIN,GMAX,IMIN,IMAX,JMIN,JMAX
      LOGICAL*4 QMPOLE
      DATA IMB,JMB/2*0/, SKIP/0/
!****
!     IMA = 40*izone(iband) izone=1,2,3,4,4,3,2,1
      JMA = 80
      IMB = INB
      JMB = JNB
      SKIP = SKIB

      IF(IMB.lt.1 .or. IMB.gt.720 .or. JMB.lt.1 .or. JMB.gt.361)
     *   GO TO 400
!****
!**** Partitions in the I direction
!**** RIA = longitude in degrees of right edge of grid box IA on grid A
!**** RIB = longitude in degrees of right edge of grid box IB of grid B
!**** IMIN(IB) = box on grid A containing left edge of box IB on B
!**** IMAX(IB) = box on grid A containing right edge of box IB on B
!**** FMIN(IB) = fraction of box IMIN(IB) on A that is left of box IB
!**** FMAX(IB) = fraction of box IMAX(IB) on A that is right of box IB
!****
      do izone=1,4
      DIA = 360d0/(izone*40)
      DIB = 360d0/IMB
      IA  = 1
      RIA = IA*DIA - 360
      IB  = IMB
      DO 150 IBP1=1,IMB
      RIB = (IBP1-1+OFFIB)*DIB
  110 IF(RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GO TO 110
!**** Right edge of A box IA and right edge of B box IB coincide
  130 IMAX(IB,izone) = IA
      FMAX(IB,izone) = 0.
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBP1,izone) = IA
      FMIN(IBP1,izone) = 0.
      GO TO 150
!**** A box IA contains right edge of B box IB
  140 IMAX(IB,izone) = IA
      FMAX(IB,izone) = (RIA-RIB)/DIA
      IMIN(IBP1,izone) = IA
      FMIN(IBP1,izone) = 1.-FMAX(IB,izone)
  150 IB = IBP1

      IMAX(IMB,izone) = IMAX(IMB,izone) + 40*izone
!       WRITE (0,*) 'Zone=',izone,' IMA=',40*izone
!       WRITE (0,*) 'IMIN=',(IMIN(I,izone),I=1,IMB)
!       WRITE (0,*) 'IMAX=',(IMAX(I,izone),I=1,IMB)
!       WRITE (0,*) 'FMIN=',(FMIN(I,izone),I=1,IMB)
!       WRITE (0,*) 'FMAX=',(FMAX(I,izone),I=1,IMB)
      end do
!****
!**** Partitions in the J direction
!****
!**** RJA = latitude in radians at top edge of box JA on grid A
!**** SINA(JA) = sine of latitude of top edge of box JA on grid A
      SINA(0)  = -1. ; sband = SINA(0) ! sine of northern edge of band
      do izone=1,4
        dsband = .1d0 * izone
        do jzs=1,10
          ja = jzs + (izone-1)*10
          SINA(JA) = sband + .1d0*jzs*dsband
        end do
        sband = sband + dsband
      end do
      SINA(40)  = 0.
      do ja=41,80
        sina(ja) = -sina(80-ja)
      end do
      ja=1
      ibefor(1) = 0
      do izone=1,8
      iz = izone ; if(iz>4) iz = 9-izone
      do jzs=1,10
      ibefor(ja+1) = ibefor(ja) + 40*iz
      ja=ja+1
      end do
      end do
!       WRITE (0,*) 'JA, sina(ja)           ',0,sina(0)
!       do ja=1,80
!         WRITE (0,*) 'JA, sina(ja), ibef(ja)=',ja,sina(ja),ibefor(ja)
!       end do
!       WRITE (0,*) 'JA,           ibef(ja)=',81,ibefor(81)
!**** RJB = latitude in radians at top edge of box JB on grid B
!**** SINB(JB) = sine of latitude of top edge of box JB on grid B
      OFFJB = (DIVJB-JMB)/2.
      DJB   = .5*TWOPI/DIVJB
      DO 220 JB=1,JMB-1
      RJB = (JB+OFFJB)*DJB - .25*TWOPI
  220 SINB(JB) = DSIN(RJB)
      SINB(0)  = -1.
      SINB(JMB)=  1.
!****
!**** JMIN(JB) = index of box of A that contains bottom edge of box JB
!**** JMAX(JB) = index of box of A that contains top edge of box JB
!**** GMIN(JB) = fraction of box JMIN(JB) on A grid that is below box JB
!**** GMAX(JB) = fraction of box JMAX(JB) on A grid that is above box JB
!****
      JMIN(1) = 1
      GMIN(1) = 0.
      JA = 1
      DO 350 JB=1,JMB-1
  310 IF(SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GO TO 310
!**** Top edge of A box JA and top edge of B box JB coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0.
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0.
      GO TO 350
!**** A box JA contains top edge of B box JB
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA)-SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB)-SINA(JA-1)
  350 CONTINUE
      JMAX(JMB) = JMA
      GMAX(JMB) = 0.
!       WRITE (0,*) 'JMIN=',(JMIN(J),J=1,JMB)
!       WRITE (0,*) 'JMAX=',(JMAX(J),J=1,JMB)
!       WRITE (0,*) 'GMIN=',(GMIN(J),J=1,JMB)
!       WRITE (0,*) 'GMAX=',(GMAX(J),J=1,JMB)
!       WRITE (0,*) 'Ibef=',(Ibefor(J),J=1,JMB)
      RETURN
!****
!**** Invalid parameters or dimensions out of range
!****
  400 WRITE (0,940) IMB,JMB,OFFIB,DIVJB, SKIP
      STOP 400
  940 FORMAT ('0Arguments received by EANTRP0 in order:'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid boxes from',
     *                    ' IDL to left edge of grid box I=1'/
     *  E24.8,' = DIVJB   = number of whole grid boxes from SP to NP'/
     *  E24.8,' = SKIP    = value to be put in B array when B',
     *  ' grid box is subset of A grid boxes with WTA = 0'/
     *  '0These arguments are invalid or out of range.')
!****

      ENTRY EANTRP (WTA,A,B)
!****
!**** EANTRP performs the horizontal interpolation
!**** Input: WTA = weighting array for values on the A grid
!****          A = per unit area or per unit mass quantity
!**** Output:  B = horizontally interpolated quantity on B grid
!****
      QMPOLE = .FALSE.
      GO TO 500

      ENTRY EANTRPP (WTA,A,B)
!****
!**** EANTRPP is similar to EANTRP but polar values are replaced by
!**** their longitudinal mean
!****
      QMPOLE = .TRUE.
!****
!**** Interpolate the A grid onto the B grid
!****
  500 DO 520 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      DO 520 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0.
      VALUE = 0.
      DO 510 JA=JAMIN,JAMAX
      izone = (ja + 9)/10
      if(izone > 4) izone=9-izone
      ima = izone*40
      IAMIN = IMIN(IB,izone)
      IAMAX = IMAX(IB,izone)
      G = SINA(JA)-SINA(JA-1)
      IF(JA.eq.JAMIN)  G = G - GMIN(JB)
      IF(JA.eq.JAMAX)  G = G - GMAX(JB)
      DO 510 IAREV=IAMIN,IAMAX
      IA  = 1+MOD(IAREV-1,IMA)
      IJA = IA + Ibefor(JA)
      F   = 1.
      IF(IAREV.eq.IAMIN)  F = F - FMIN(IB,izone)
      IF(IAREV.eq.IAMAX)  F = F - FMAX(IB,izone)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
  510 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = SKIP
      IF(WEIGHT.ne.0.)  B(IJB) = VALUE/WEIGHT
  520 continue
!****
!**** Replace individual values near the poles by longitudinal mean
!****
      IF(.NOT.QMPOLE)  RETURN
      DO 630 JB=1,JMB,JMB-1
      WEIGHT = 0.
      VALUE  = 0.
      DO 610 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      IF(B(IJB).eq.SKIP)  GO TO 610
      WEIGHT = WEIGHT + 1.
      VALUE  = VALUE  + B(IJB)
  610 continue
      BMEAN = SKIP
      IF(WEIGHT.ne.0.)  BMEAN = VALUE/WEIGHT
      DO 620 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
  620 B(IJB) = BMEAN
  630 continue
      RETURN
      END SUBROUTINE EANTRP0

      subroutine get_ij(lat,lon,ifirst,ij)
!     IMPLICIT real*4 (A-H, O-Z)
! Given the latitude and longitude of the center, get_ij finds the
! position of the grid box if they are arranged in the order:
!     latitudes 90S to 90N, and 180W to 180E within each latitude

      integer,save :: lats(81),nbefor(81)

      if(ifirst==1) then
! Set LATS(j),   the j-th Southern latitude edge in .01 degrees
! and IBEFOR(j), the number of grid boxes South of LATS(j)
        xbypi = 9000d0/asin(1d0)        ! 100 * 180/pi
        j=1 ; nbefor(1)=0 ; sband = -1. ! sine of southern edge of band
        do iband=1,8
          iz = iband ; if(iz>4) iz = 9-iband
          dsband = .1d0 * iz
          do jzs=1,10
            nbefor(j+1) = nbefor(j) + 40*iz
            LATS(J) = nint(xbypi*asin(sband + .1d0*(jzs-1)*dsband))
            j = j + 1
          end do
          sband = sband + dsband
        end do
        LATS(81)  = -lats(1)   ! 9000
!          do j=1,81
!            WRITE (0,*) 'J,lats,nbefor=',j,lats(j),nbefor(j)
!          end do
        return
      end if

!**** find "j" 1 -> 80 using a binary search
      jl = 1 ; jh = 80
   10 jm = (jl+jh)/2
      if(lat<lats(jm)) then
        jh = jm-1
      else
        jl = jm
      end if
      if(lat>lats(jl+1).and.lat<lats(jh)) go to 10
      j=jl
      if(lat>lats(jh)) j=jh

!**** find zone and "i"
      izon = (j+9)/10          ! 1 2 3 4 4 3 2 1
      if(izon>4) izon = 9-izon
      idlon = 900/izon         ! 100 * 360/(40*izon)
      i = 1 + ( lon + 18000 )/idlon

      ij = i + nbefor(j)

      return
      end subroutine get_ij

!*********
!*********
!
! Done.
!



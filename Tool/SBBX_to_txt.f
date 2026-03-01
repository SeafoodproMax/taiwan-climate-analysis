C**** SBBX_to_txt.f ; to compile: f90-compiler SBBX_to_txt.f
C****  options needed: -fconvert=big-endian -frecord-marker=4

C**** This program expects 1 (trimmed) file:
C**** Record 1 starts with 8 integers I1-I8 and an 80-byte TITLE.
C**** All further records start with 7 integers N1-N7,
C****             a real number R, followed by a data array (real*4).
C**** I1 or N1 is the length of the data array in the NEXT record.
C**** Unless its length is 0, each data-array contains a time series
C**** of monthly T-anomalies (C) starting with January of year I6 for
C**** the grid box. N2,N3,N4,N5 indicate the edges in .01 degrees of
C**** that grid box in the order: latitude of southern edge, latitude of
C**** northern edge, longitude of western edge, longit. of eastern edge.
C**** The world is covered with 8000 equal area grid boxes, so each
C**** file has 8001 records. I7 is the flag for missing data (9999).
C****
      CHARACTER*80 TITLE
      INTEGER INFO(8)
      REAL*4, allocatable :: TIN(:)

C**** Read in file name, open in- and output file
      if(iargc() .ne. 1) then
        write(*,*) 'Usage: SBBX_to_txt sbbx-file_name'
        stop
      end if
      call getarg(1,title)
      open(8,file=title,form='unformatted',status='old',
     *    access='sequential',err=900)
      open(12,file=trim(title)//'.txt',form='formatted')

C**** Read, display, and use the relevant parts of the header record
      READ(8) INFO,TITLE
      WRITE(6,*) TITLE
      Mnow=INFO(1)       ! length of first data record (months)
      MONM=INFO(4)       ! max length of time series (months)
      IYRBEG=INFO(6)     ! beginning of time series (calendar year)
      BAD=INFO(7)
      write(*,*) 'missing_data flag:',INFO(7)
      IYREND=IYRBEG-1+MONM/12        ! last calendar year of timeseries

C**** Output file headers
        write(12,'(a80i5,a1,i4,a)') title,iyrbeg,'-',iyrend,
     *   ' varying base period'
        write(12,'(a)') ' latitude-rnge longitude-range year     Jan'//
     *  '     Feb     Mar     Apr     May     Jun     Jul     Aug'//
     *  '     Sep     Oct     Nov     Dec'

C**** Loop over subboxes
      allocate (tin(monm))
      DO 100 N=1,8000
C**** Read in time series TIN of monthly means
      tin = bad

      CALL SREAD (8,TIN,Mnow,LATS,LATN,LONW,LONE,DL,next)
      Mnow=next

      i1=1
      do iy=iyrbeg,iyrend
        iok=0
        do i=i1,i1+11
          if(tin(i).ne.bad) iok=iok+1
        end do
        if(iok>0) write(12,'(2f7.2,2f8.2,i5,12f8.2)')
     *    .01*LATS,.01*latn,.01*LONW,.01*lone, iy,(tin(i),i=i1,i1+11)
        i1=i1+12
      end do

  100 CONTINUE
C**** End of loop over subboxes
      STOP
  900 write(*,*) 'readd error for file ',title
      stop
      END

      SUBROUTINE SREAD (NDISK,ARRAY,LEN, N1,N2,N3,N4, DSTN,LNEXT)
      REAL ARRAY(LEN)
      READ(NDISK) LNEXT, N1,N2,N3,N4, NR1,NR2, DSTN, ARRAY
      RETURN
      END

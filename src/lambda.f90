subroutine lambda_setup(ncores,param)
! Allocate, zero, and deallocate various arrays needed for MST & lambda.
! File I/O for obj star masses, lambda values, edge lengths for CDFs, etc.

  USE directories_module
  use lambdaparams_module
  implicit none
  
  integer, intent(in) :: ncores
  double precision, intent(in) :: param(1:ncores)
  integer :: i,j
  integer :: nedge       ! number of edges in MST (nmst-1)
  integer :: nlam        ! number of lambda types
  character(len=5) :: Nmedstr
  LOGICAL :: dirExists
  
!lambda: uses average random mst length & total obj length
!lambda bar: uses mean edge lengths (arithmetic mean)
!            Equivalent to lambda as 1/n cancels
!lambda rms: uses root mean square. p=2 in 'generalised mean' eq.
!lambda smr: uses square mean root (if that's a thing...). p=1/2.
!lambda har: generalised mean with p=-1
!lambda til: uses median edge length
!lambda star: uses mst made of median edge length, then add total MST
!gamma: uses geometric mean
!lambda ln: uses generalised f-mean, where f^(-1) is ln
!lambda N median: takes the 2 or 3 median points (2 if even number in MST)
!                 then finds the mean of these

!======================================================
! Set # of stars in MST, # of random MSTs, & # of CDFs
!======================================================
  
  nmst=20 !number of stars in the MST
  nmed=3 !when using x number of median points, findlamNmed
         !odd for odd nedge, even nmst & v/v)
  
!For large nmst it can take a while to build the MST,
! so reduce the number of times we do the loop -
! higher nmst MSTs are less stochastic anyway...
  nloop = 1000
  IF(nmst >= 100) nloop = 50
  
  
  ! create directory for lambda data:
  lampath = trim(paramdir)//'/lambda'
  INQUIRE(file = trim(lampath), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(lampath))
  END IF
  
  nCDF = 20 !Number of CDFs plotted for random MSTs; must be =< nloop
  
  ! create directory for CDF data (edge lengths in MST for different subsets)
  CDFPath = trim(paramdir)//'/CDFdata'
  INQUIRE(file = TRIM(CDFPath), exist = dirExists)
!(Works for gfortran. For ifort: ...directory=newDir,exist...)
  IF (.NOT. dirExists) THEN
     WRITE(6,'(a)') "Creating new directory: '"//TRIM(CDFPath)//"'"
     CALL system('mkdir -p '//TRIM(CDFPath))
  END IF
  
!======================================================
  
  
! Lengths of the edges of object stars:
  ALLOCATE(edgelengths(1:nmst-1))


!*******************************************************!
! Mass segregation
!
!########################################  
!Set up different types types of average
!######################################## 

! Find the degree of mass segregation, lambda, using
! average MST edge lengths, with different types of average.
! Average edge lengths for object stars and
! random stars are saved for each snapshot.

  nlam=0 !number of different types of lambda
  
  findlam = .FALSE.
  findlambar = .TRUE.
  findlamrms = .TRUE.
  findlamsmr = .FALSE.
  findlamhar = .FALSE.
  findlamtil = .FALSE.
  findlamNmed = .TRUE.
  findlamstar = .FALSE.
  findgam = .TRUE.
  findlamln = .FALSE.
  
  IF (findlam) THEN
     nlam = nlam + 1
     
! lambda: final lambda value calculated (output)
! l_low: +ve error bar
! l_up: -ve error bar
! l_avranmst: median average edge length of random MTSs
! l_objmst: average edge length of object MST
! l_allmsts: lambda values for all the random MSTs that are generated
     allocate(l_allmsts(1:nloop))
     
     lambda=0.
     l_up=0.
     l_low=0.
     l_avranmst=0.
     l_objmst=0.
     l_allmsts=0.
  END IF

  IF (findlambar) THEN
     nlam = nlam + 1
     allocate(lbar_allmsts(1:nloop))
     lambda_bar=0.
     l_up_bar=0.
     l_low_bar=0.
     lbar_avranmst=0.
     lbar_objmst=0.
     lbar_allmsts=0.
  END IF

  IF (findlamrms) THEN
     nlam = nlam + 1
     allocate(lrms_allmsts(1:nloop))
     lambda_rms=0.
     l_up_rms=0.
     l_low_rms=0.
     lrms_avranmst=0.
     lrms_objmst=0.
     lrms_allmsts=0.
  END IF

  IF (findlamsmr) THEN
     nlam = nlam + 1
     allocate(lsmr_allmsts(1:nloop))
     lambda_smr=0.
     l_up_smr=0.
     l_low_smr=0.
     lsmr_avranmst=0.
     lsmr_objmst=0.
     lsmr_allmsts=0.
  END IF

  IF (findlamhar) THEN
     nlam = nlam + 1
     allocate(lhar_allmsts(1:nloop))
     lambda_har=0.
     l_up_har=0.
     l_low_har=0.
     lhar_avranmst=0.
     lhar_objmst=0.
     lhar_allmsts=0.
  END IF
  
  IF (findlamtil) THEN
     nlam = nlam + 1
     allocate(ltil_allmsts(1:nloop))
     lambda_til=0.
     l_up_til=0.
     l_low_til=0.
     ltil_avranmst=0.
     ltil_objmst=0.
     ltil_allmsts=0.
  END IF
  
  IF (findlamNmed) THEN
     
     !Nmed checks:
     nedge=nmst-1
     !------------
     ! Even nmst: (odd edges)
     if (.not. MOD(nedge,2)==0) then
        if (MOD(Nmed,2)==0) stop 'Nmed must be odd'
        if (.not. Nmed .gt. 1) stop 'Nmed must be greater than 1'
        if (Nmed .gt. nedge) stop 'Nmed must be less than nedge'
     else
     !-----------
     ! Odd nmst: (even edges)
        if (.not. MOD(Nmed,2)==0) stop 'Nmed must be even'
        if (.not. Nmed .gt. 2) stop 'Nmed must be greater than 2'
        if (Nmed .gt. nedge) stop 'Nmed must be less than nedge'
     end if
     
     write(Nmedstr,'(I5)') Nmed
     
     nlam = nlam + 1
     allocate(lNmed_allmsts(1:nloop))
     lambda_Nmed=0.
     l_up_Nmed=0.
     l_low_Nmed=0.
     lNmed_avranmst=0.
     lNmed_objmst=0.
     lNmed_allmsts=0.
  END IF
  
  IF (findlamstar) THEN
     nlam = nlam + 1
     allocate(lstar_allmsts(1:nloop))
     lambda_star=0.
     l_up_star=0.
     l_low_star=0.
     lstar_avranmst=0.
     lstar_objmst=0.
     lstar_allmsts=0.
  END IF
  
  IF (findgam) THEN
     nlam = nlam + 1
     allocate(lgam_allmsts(1:nloop))
     lambda_gam=0.
     l_up_gam=0.
     l_low_gam=0.
     lgam_avranmst=0.
     lgam_objmst=0.
     lgam_allmsts=0.
  END IF
  
  IF (findlamln) THEN
     nlam = nlam + 1
     allocate(lln_allmsts(1:nloop))
     lambda_ln=0.
     l_up_ln=0.
     l_low_ln=0.
     lln_avranmst=0.
     lln_objmst=0.
     lln_allmsts=0.
  END IF
  
  IF(nlam==0) THEN
     write(6,*) 'All lambda methods set to FALSE; need at least one TRUE'
     stop
  END IF
  
!====================================
  
  write(6,*)"       Calculating lambda..."
  
  
  call find_lambda(ncores,param)

  write(6,*)"       ...done"
  
  
!**********************
! Lambda outputs
!**********************
  
  IF (findlam) THEN
! Median random MST edge lengths &
! object edge lengths used for each lambda, and lambda values with errors:
     OPEN(12,file=trim(lampath)//'/MST_lambda.dat',status='replace')
     WRITE(12,13) l_avranmst,l_objmst,lambda,l_low,l_up
     CLOSE(12)
     
! Lambda CDF: lambda bar values for every random MST at each snapshot
     open(12,file=TRIM(CDFpath)//'/allMSTs_lam.dat',status='replace')
     write(12,14) l_allmsts(1:nloop) 
     close(12)
     
  END IF
  
  
! Lambda bar: average MST, object MST, lambda bar value, and errors
  IF (findlambar) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_lambar.dat',status='replace')
     WRITE(12,13) lbar_avranmst,lbar_objmst, &
             & lambda_bar,l_low_bar,l_up_bar
     CLOSE(12)
     
! Lambda bar CDF: lambda bar values for every random MST at each snapshot
     open(12,file=TRIM(CDFpath)//'/allMSTs_lambar.dat',status='replace')
     write(12,14) lbar_allmsts(1:nloop)
     close(12)
     
  END IF
  
  
  IF (findlamrms) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_lamrms.dat',status='replace')
     WRITE(12,13) lrms_avranmst,lrms_objmst, &
             & lambda_rms,l_low_rms,l_up_rms
     CLOSE(12)
     
     open(12,file=TRIM(CDFpath)//'/allMSTs_lamrms.dat',status='replace')
     write(12,14) lrms_allmsts(1:nloop)
     close(12)
     
  END IF

  
  IF (findlamsmr) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_lamsmr.dat',status='replace')
     WRITE(12,13) lsmr_avranmst,lsmr_objmst, &
          & lambda_smr,l_low_smr,l_up_smr
     CLOSE(12)
     
     open(12,file=TRIM(CDFpath)//'/allMSTs_lamsmr.dat',status='replace')
     write(12,14) lsmr_allmsts(1:nloop)
     close(12)
     
  END IF

  
  IF (findlamhar) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_lamhar.dat',status='replace')
     WRITE(12,13) lhar_avranmst,lhar_objmst, &
          & lambda_har,l_low_har,l_up_har
     CLOSE(12)
     
     open(12,file=TRIM(CDFpath)//'/allMSTs_lamhar.dat',status='replace')
     write(12,14) lhar_allmsts(1:nloop)
     close(12)
     
  END IF

  
  IF (findlamtil) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_lamtil.dat',status='replace')
     WRITE(12,13) ltil_avranmst,ltil_objmst, &
          & lambda_til,l_low_til,l_up_til
     CLOSE(12)
     
     open(12,file=TRIM(CDFpath)//'/allMSTs_lamtil.dat',status='replace')
        write(12,14) ltil_allmsts(1:nloop)
     close(12)
     
  END IF

  
  IF (findlamNmed) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_lam'//&
          & trim(adjustl(Nmedstr))//'med.dat',status='replace')
     WRITE(12,13) lNmed_avranmst,lNmed_objmst, &
          & lambda_Nmed,l_low_Nmed,l_up_Nmed
     CLOSE(12)
     
     open(12,file=TRIM(CDFpath)//'/allMSTs_lam_'//&
          & trim(adjustl(Nmedstr))//'med.dat',status='replace')
     write(12,14) lNmed_allmsts(1:nloop)
     close(12)
     
  END IF

  
  IF (findlamstar) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_lamstar.dat',status='replace')
     WRITE(12,13) lstar_avranmst,lstar_objmst, &
          & lambda_star,l_low_star,l_up_star
     CLOSE(12)
     
     open(12,file=TRIM(CDFpath)//'/allMSTs_lamstar.dat',status='replace')
     write(12,14) lstar_allmsts(1:nloop)
     close(12)
     
  END IF

  
  IF (findgam) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_gam.dat',status='replace')
     WRITE(12,13) lgam_avranmst,lgam_objmst, &
          & lambda_gam,l_low_gam,l_up_gam
     CLOSE(12)
     
     open(12,file=TRIM(CDFpath)//'/allMSTs_gam.dat',status='replace')
     write(12,14) lgam_allmsts(1:nloop)
     close(12)
     
  END IF
  
  
  IF (findlamln) THEN
     
     OPEN(12,file=trim(lampath)//'/MST_lamln.dat',status='replace')
     WRITE(12,13) lln_avranmst,lln_objmst, &
          & lambda_ln,l_low_ln,l_up_ln
     CLOSE(12)
     
     open(12,file=TRIM(CDFpath)//'/allMSTs_lamln.dat',status='replace')
     write(12,14) lln_allmsts(1:nloop)
     close(12)

  END IF
  
13 FORMAT(2(2X,F10.5),3(2X,F9.3))
14 format(*(2X,F10.5))


  !===========================================
  
  deallocate(edgelengths)
  IF (findlam) deallocate(l_allmsts)
  IF (findlambar) deallocate(lbar_allmsts)
  IF (findlamrms) deallocate(lrms_allmsts)
  IF (findlamsmr) deallocate(lsmr_allmsts)
  IF (findlamhar) deallocate(lhar_allmsts)
  IF (findlamtil) deallocate(ltil_allmsts)
  IF (findlamNmed)deallocate(lNmed_allmsts)
  IF (findlamstar) deallocate(lstar_allmsts)
  IF (findgam) deallocate(lgam_allmsts)
  IF (findlamln) deallocate(lln_allmsts)

end subroutine lambda_setup


  
SUBROUTINE find_lambda(ncores,param)
  USE directories_module
  USE lambdaparams_module

  IMPLICIT NONE

!-----------
! star data
!-----------
  INTEGER, INTENT(in) :: ncores
! param = parameter called from main program (mass, temperature,etc)
  DOUBLE PRECISION,intent(in) :: param(1:ncores)
  INTEGER :: IDs(1:ncores) !ID numbers for heapsort
  
!-----------
! MST stuff
!-----------
  INTEGER :: nedge !number of edge lengths
  INTEGER :: nint_nedged2, ceint_nedged2 !NINT(... & CEILING(  REAL(nedge)/2.)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: edgeL !length of each MST edge
  INTEGER, DIMENSION(:), ALLOCATABLE :: edgeLlist !IDs of 'edgeL' entries
  double precision, dimension(:,:), allocatable :: connections
  double precision :: totallength  !total length of mst (sum of 'edgeL')
  double precision :: medianlength !median edge length in MST
  double precision :: mst_RA(1:nmst), mst_dec(1:nmst), z(1:nmst)

!total lengths of random MSTs for different lambda:
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: l_ranmst, lbar_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lrms_ranmst, lsmr_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lhar_ranmst, ltil_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lNmed_ranmst, lstar_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lgam_ranmst, lln_ranmst

! done = record which stars have been randomly selected
  LOGICAL :: done(1:ncores)
  REAL :: rand
  
  character(len=6) :: filei !string version of 'i'
  INTEGER :: i,j,k !generic counters

  z=0.  ! z-coordinate
  
!Frequently used expressions:
  nedge=nmst-1
  !for median calculations:
  nint_nedged2=nint(real(nedge)/2.) !nint to round to nearest int from ~x.0
  ceint_nedged2=ceiling(real(nedge)/2.) !ceint to round up from x.5)

!Allocate memory for arrays of length nmst:
  ALLOCATE(edgeL(1:nedge))
  ALLOCATE(edgeLlist(1:nedge))
! Coordinates of the edge connections: (xi,yi,zi,xj,yj,zj for each edge)
  ALLOCATE(connections(1:6,1:nedge))
  

!============================================================
!Sort stars in order of mass & make selection of 'obj' stars
!============================================================
!

  IDs=0

!Assign IDs to stars for heapsort
  DO j=1,ncores
     IDs(j)=j
  END DO

  CALL heapsort(ncores,param,IDs)

!select nmst most massive stars:
!(while checking whether star should be ignored)
  j=0 !j tracks with i but changes with each *attempted* iteration,
      ! even if i does not change
  
  !write positions of object stars:
  open(1,file=trim(lampath)//'/objpositions.dat', status='replace')
  
  DO i = 1,nmst
     j=j+1
     k = IDs(ncores+1-j) !heapsort orders from small to large - invert
     mst_RA(i) = RA(k)  !RA coordinate
     mst_dec(i) = dec(k) !dec coordinate
     
! write object values selected to screen:
     write(6,'(F7.3)') param(IDs(ncores+1-j))
     !write(6,*) param(IDs(ncores+1-j))
     write(1,'(2(2X,F9.5))') mst_RA(i), mst_dec(i)
  END DO
  
  close(1)
  
!          *********************************
!*****************************************************
!      Find the MST length for the OBJECT stars
!*****************************************************
!          *********************************
  
  edgeL = 0.
  
  CALL mst(nmst,mst_RA,mst_dec,z,edgeL,connections)
  
  ! write out positions of nodes between connecting MST edges:
  call system('mkdir -p '//trim(lampath)//'/edgecoords')
  open(1,file=trim(lampath)//'/edgecoords/objconnections.dat'&
       & ,status='replace')
  do i=1,nedge
     write(1,'(6(2X,F9.5))') connections(1:6,i)
  end do
  close(1)
  
!Assign IDs for heapsort to order MST edges
  DO j=1,nedge
     edgeLlist(j)=j
  END DO
  
  CALL heapsort(nedge,edgeL,edgeLlist)

  ! Write out the ordered edge lengths of the MST to plot a CDF:
  OPEN(10,file=trim(CDFPath)//'/MSTedgeL.dat',&
       & status='replace')
  do i=1,nedge
     edgelengths(i)=edgeL(edgeLlist(i))
  end do
  write(10,11) edgelengths(1:nedge)
11 FORMAT(*(2X,F9.5))
  close(10)
  

!######################################
!Average edge lengths for object stars
!######################################

  ! Total length of object MST
  totallength=0.
  totallength=SUM(edgeL)

  !Median edge length:
  medianlength=0.
  DO i = 1, nedge
     edgeLlist(i) = i
  END DO
  CALL heapsort(nedge, edgeL, edgeLlist)

  if (MOD(nedge,2)==0) then !even nedge (odd nmst), take mean of median two
     medianlength=edgeL(edgeLlist(nint_nedged2)) &
          + edgeL(edgeLlist(nint_nedged2 + 1)) !nedge/2=x.0, ensure nearest
     medianlength=medianlength/2.
  else !.not. MOD(nedge,2)==0, odd nedge (even nmst), take median.
     medianlength=edgeL(edgeLlist(ceint_nedged2)) !nedge/2=x.5, round up
  end if

  
!Lambda MST:
! This is just the same as lambar without 1/n
  IF (findlam) l_objmst = totallength


!Lambda bar MST:
!Use the mean MST edge length
  IF (findlambar) lbar_objmst = totallength/(nedge)


!Generalised mean: ((1/n)*SUM(x**p))**(1/p). Arithmetic mean p=1
!rms p=2, smr p=1/2, harmonic p=-1
  
!Lambda rms MST:
  IF (findlamrms) THEN
     DO i=1,nedge
        lrms_objmst = lrms_objmst + edgeL(i)**2.
     END DO
  lrms_objmst = ( (1./REAL(nedge))*lrms_objmst)**0.5
  END IF

  
!Lambda smr MST:
  IF (findlamsmr) THEN
     DO i=1,nedge
        lsmr_objmst = lsmr_objmst + (edgeL(i))**0.5
     END DO
  lsmr_objmst = ( (1./REAL(nedge))*(lsmr_objmst) )**2.
  END IF

  
!Lambda har MST:
  IF (findlamhar) THEN
     DO i=1,nedge
        lhar_objmst = lhar_objmst + (edgeL(i))**(-1.)
     END DO
  lhar_objmst = ( (1./REAL(nedge))*(lhar_objmst) )**(-1.)
  END IF
  
  
!Lambda tilde MST:
  ! Use the median edge length
  IF (findlamtil) ltil_objmst = medianlength
  
  
!Lambda N-median MST:
  IF (findlamNmed) THEN
     if (mod(Nmed,2)==0) then !even edge lengths (odd nmst), round to nint
        lNmed_objmst=medianlength*2 !i=1. *2 for l1+l2 as /2 for median.
        do i=2,Nmed/2
           !take values either side of two median:
           lNmed_objmst=lNmed_objmst &
                + edgeL(edgeLlist(nint_nedged2 - (i-1))) &
                + edgeL(edgeLlist(nint_nedged2 + i))
        end do
     else !odd edge lengths (even nmst), round to ceiling
        lNmed_objmst=medianlength !i=1
        do i=2,(Nmed+1)/2        !take values either side of median:
           lNmed_objmst=lNmed_objmst &
                + edgeL(edgeLlist(ceint_nedged2 - (i-1))) &
                + edgeL(edgeLlist(ceint_nedged2 + (i-1)))
        end do
     end if
     lNmed_objmst=lNmed_objmst/real(Nmed)
  END IF
  
  
!Lambda star MST:
! Find the median length in the tree, & find length of a tree made from these
  IF (findlamstar) THEN
     lstar_objmst = (nedge)*medianlength
     
! Then add on the actual length of the tree
     lstar_objmst = lstar_objmst+totallength
  END IF


!Generalised f-mean: f**(-1) * ((1/n)*SUM(f(x)) where f is some function
!Gamma is a special case where f=ln, so f**(-1)=exp

!Gamma MST:
  IF (findgam) THEN
     DO i = 1,nedge
        lgam_objmst = lgam_objmst + LOG(edgeL(i))  !Add edges to
     END DO                                                  !find total length
!Calculate the geometric mean
     lgam_objmst = EXP( (1./REAL(nedge)) * lgam_objmst)
  END IF

  
!Lambda ln MST:
  IF (findlamln) THEN
     DO i = 1,nedge
        lln_objmst = lln_objmst + EXP(edgeL(i))  !Add edges to
     END DO                                                !find total length
!Calculate the geometric mean
     lln_objmst = LOG( (1./REAL(nedge)) * lln_objmst)
  END IF


!          *********************************
!*****************************************************
!    Find the MST length for the nloop random MSTs
!*****************************************************
!          *********************************

  edgeL = 0.
  
  !set lengths of random MSTs to 0:
  IF (findlam) then
     ALLOCATE(l_ranmst(1:nloop))
     l_ranmst = 0. 
  END IF
  IF (findlambar) THEN
     ALLOCATE(lbar_ranmst(1:nloop))
     lbar_ranmst = 0.
  END IF
  IF (findlamrms) THEN
     ALLOCATE(lrms_ranmst(1:nloop))
     lrms_ranmst = 0.
  END IF
  IF (findlamsmr) THEN
     ALLOCATE(lsmr_ranmst(1:nloop))
     lsmr_ranmst = 0.
  END IF
  IF (findlamhar) THEN
     ALLOCATE(lhar_ranmst(1:nloop))
     lhar_ranmst = 0.
  END IF
  IF (findlamtil) THEN
     ALLOCATE(ltil_ranmst(1:nloop))
     ltil_ranmst = 0.
  END IF
  IF (findlamNmed) THEN
     ALLOCATE(lNmed_ranmst(1:nloop))
     lNmed_ranmst = 0.
  END IF
  IF (findlamstar) THEN
     ALLOCATE(lstar_ranmst(1:nloop))
     lstar_ranmst = 0.
  END IF
  IF (findgam) THEN
     ALLOCATE(lgam_ranmst(1:nloop))
     lgam_ranmst = 0.
  END IF
  IF (findlamln) THEN
     ALLOCATE(lln_ranmst(1:nloop))
     lln_ranmst = 0.
  END IF

  

  DO j = 1,nloop          !Do nloop random MSTs
     mst_RA = 0. ; mst_dec = 0.
     done = .FALSE.

     DO i = 1,nmst        !Select nmst random stars
6       CALL RANDOM_NUMBER(rand) !random number between 0. and 1.
        k = NINT(rand*ncores)
        IF (k == 0) GOTO 6    !There is no star with id=0
        IF (done(k)) GOTO 6    !Don't chose the same star twice
        done(k) = .TRUE.       !You're selected

        mst_RA(i) = RA(k)
        mst_dec(i) = dec(k)
     END DO

     CALL mst(nmst,mst_RA,mst_dec,z,edgeL,connections)
  
!######################################
!Average edge lengths for random stars
!######################################

     totallength=0.
     !DO i = 1,nedge
     !   totallength = totallength + edgeL(i)  !Add the edges of the mst
     !END DO                                    ! to find the total length
     totallength=SUM(edgeL)

     medianlength=0.
     DO i = 1, nedge
        edgeLlist(i) = i
     END DO
     CALL heapsort(nedge, edgeL, edgeLlist)
     
     !CDFs of random stars:
     !if(thisproj=='3D') then
     if (j.le.nCDF) then
        
        ! open files for MSTs of randomly selected stars:
        write(filei,'(I6)') j
        open(20+j,file=trim(CDFPath)//'/MSTedgeL_'&
             & //trim(adjustl(filei))//'.dat',status='replace')
        
        ! Write out the edge lengths of the MST to plot a CDF:
        do k=1,nedge
           edgelengths(k)=edgeL(edgeLlist(k))
        end do
        write(20+j,20) edgelengths(1:nedge)
20      FORMAT(*(2X,F9.5))
        
        !close CDF file:
        close(20+j)
        
     end if
     !end if
     
     !Median edge length of this MST (depends on even/odd nedge):
     if (MOD(nedge,2)==0) then
        ! if even no. of edge lengths, take mean of median 2.
        medianlength=edgeL(edgeLlist(nint_nedged2 +1)) &
             + edgeL(edgeLlist(nint_nedged2))
        medianlength=medianlength/2.
     else !.not. MOD(nedge,2)==0
        !if odd no. of edge lengths (even nmst), take median edge length.
        !(Use CEILING as always need to round up from #.5)
        medianlength=edgeL(edgeLlist(ceint_nedged2))
     end if
     
!Lambda MST:
! This is just the same as lambar without 1/n
     IF (findlam) l_ranmst(j) = totallength


!Lambda bar MST:
!Use the mean MST edge length
     IF (findlambar) lbar_ranmst(j) = totallength/(nedge)


!Generalised mean: ((1/n)*SUM(x**p))**(1/p). Arithmetic mean p=1
!rms p=2, smr p=1/2, harmonic p=-1
  
!Lambda rms MST:
     IF (findlamrms) THEN
        DO i=1,nedge
           lrms_ranmst(j) = lrms_ranmst(j) + (edgeL(i))**2.
        END DO
     lrms_ranmst(j) = ( (1./REAL(nedge))*(lrms_ranmst(j)) )**0.5
     END IF

  
!Lambda smr MST:
     IF (findlamsmr) THEN
        DO i=1,nedge
           lsmr_ranmst(j) = lsmr_ranmst(j) + (edgeL(i))**0.5
        END DO
     lsmr_ranmst(j) = ( (1./REAL(nedge))*(lsmr_ranmst(j)) )**2.
     END IF

  
!Lambda har MST:
     IF (findlamhar) THEN
        DO i=1,nedge
           lhar_ranmst(j) = lhar_ranmst(j) + (edgeL(i))**(-1.)
        END DO
     lhar_ranmst(j) = ( (1./REAL(nedge))*(lhar_ranmst(j)) )**(-1.)
     END IF


!Lambda tilde MST:
! Use the median edge length
     IF (findlamtil) ltil_ranmst(j) = medianlength


!Lambda N-median MST:
!Find the N median points and take the mean of these
     IF (findlamNmed) THEN
        if (mod(Nmed,2)==0) then !even nedge (odd nmst)
           lNmed_ranmst(j) = medianlength*2 !i=1. *2 as /2 earlier.
           do i=2,Nmed/2
              !take values either side of two used for median:
              lNmed_ranmst(j)=lNmed_ranmst(j) &
                   + edgeL(edgeLlist(nint_nedged2 - (i-1))) &
                   + edgeL(edgeLlist(nint_nedged2 + i))
           end do
        else !.not. mod(Nmed,2)==0, odd nedge (even nmst)
           lNmed_ranmst(j) = medianlength !i=1
           do i=2,(Nmed+1)/2
              !take values either side of median:
              lNmed_ranmst(j)=lNmed_ranmst(j) &
                   + edgeL(edgeLlist(ceint_nedged2 - (i-1))) &
                   + edgeL(edgeLlist(ceint_nedged2 + (i-1)))
           end do
        end if
        lNmed_ranmst(j)=lNmed_ranmst(j)/real(Nmed)
     END IF


!Lambda star MST:
! Find the median length in the tree, & find length of a tree made from these
     IF (findlamstar) THEN
        lstar_ranmst(j) = (nedge) * medianlength
        
! Then add on the actual length of the tree
        lstar_ranmst(j) = lstar_ranmst(j) + totallength
     END IF


!Generalised f-mean: f**(-1) * ((1/n)*SUM(f(x)) where f is some function
!Gamma is a special case where f=ln, so f**(-1)=exp

!Gamma MST:
     IF (findgam) THEN
        DO i = 1,nedge
           lgam_ranmst(j) = lgam_ranmst(j) + LOG(edgeL(i))  !Add edges to
        END DO                                               !find total length
        lgam_ranmst(j) = EXP( (1./REAL(nedge)) * lgam_ranmst(j) )
     END IF

  
!Lambda ln MST:
     IF (findlamln) THEN
        DO i = 1,nedge
           lln_ranmst(j) = lln_ranmst(j) + EXP(edgeL(i))  !Add edges to
        END DO                                                !find total length
!Calculate the geometric mean
        lln_ranmst(j) = LOG( (1./REAL(nedge)) * lln_ranmst(j))
     END IF
     
  END DO !end of nloop

!========================================================================
!Finally - find the average value of the random MSTs and the significance
!========================================================================
!
!The average is found using the median value - this means that small
!values don't skew the value.
!The signifincance is found using the 1/6 and 5/6 boundaries


!Calculate lambda:
  IF (findlam) THEN
     CALL calc_lambda(l_objmst, l_ranmst, lambda, &
          & l_allmsts(1:nloop), l_low, l_up, &
          & l_avranmst)
  END IF
  
  
!Calculate lambda bar:
  IF (findlambar) THEN
     CALL calc_lambda(lbar_objmst, lbar_ranmst, lambda_bar, &
          & lbar_allmsts(1:nloop), l_low_bar, l_up_bar, &
          & lbar_avranmst)
  END IF
  
  
!Calculate lambda rms:
  IF (findlamrms) THEN
     CALL calc_lambda(lrms_objmst, lrms_ranmst, lambda_rms, &
          & lrms_allmsts(1:nloop), l_low_rms, l_up_rms, &
          & lrms_avranmst)
  END IF
  

!Calculate lambda smr:
  IF (findlamsmr) THEN
     CALL calc_lambda(lsmr_objmst, lsmr_ranmst, lambda_smr, &
          & lsmr_allmsts(1:nloop), l_low_smr, l_up_smr, &
          & lsmr_avranmst)
  END IF
  
  
!Calculate lambda har:
  IF (findlamhar) THEN
     CALL calc_lambda(lhar_objmst, lhar_ranmst, lambda_har, &
          & lhar_allmsts(1:nloop), l_low_har, l_up_har, &
          & lhar_avranmst)
  END IF
  
  
!Calculate lambda tilde:
  IF (findlamtil) THEN
     CALL calc_lambda(ltil_objmst, ltil_ranmst, lambda_til, &
          & ltil_allmsts(1:nloop), l_low_til, l_up_til, &
          & ltil_avranmst)
  END IF
  
  
!Calculate lambda N median:
  IF (findlamNmed) THEN
     CALL calc_lambda(lNmed_objmst, lNmed_ranmst, lambda_Nmed, &
          & lNmed_allmsts(1:nloop), l_low_Nmed, l_up_Nmed, &
          & lNmed_avranmst)
  END IF
  
  
!Calculate lambda star:
  IF (findlamstar) THEN
     CALL calc_lambda(lstar_objmst, lstar_ranmst, lambda_star, &
          & lstar_allmsts(1:nloop), l_low_star, l_up_star, &
          & lstar_avranmst)
  END IF
  
  
!Calculate gamma:
  IF (findgam) THEN
     CALL calc_lambda(lgam_objmst, lgam_ranmst, lambda_gam, &
          & lgam_allmsts(1:nloop), l_low_gam, l_up_gam, &
          & lgam_avranmst)
  END IF
  
  
!Calculate lambda ln:
  IF (findlamln) THEN
     CALL calc_lambda(lln_objmst, lln_ranmst, lambda_ln, &
          & lln_allmsts(1:nloop), l_low_ln, l_up_ln, &
          & lln_avranmst)
  END IF
  
  
!Deallocate arrays:
  DEALLOCATE(edgeL)
  DEALLOCATE(edgeLlist)
  deallocate(connections)
  IF (findlam) DEALLOCATE(l_ranmst)
  IF (findlambar) DEALLOCATE(lbar_ranmst)
  IF (findlamrms) DEALLOCATE(lrms_ranmst)
  IF (findlamsmr) DEALLOCATE(lsmr_ranmst)
  IF (findlamhar) DEALLOCATE(lhar_ranmst)
  IF (findlamtil) DEALLOCATE(ltil_ranmst)
  IF (findlamNmed) DEALLOCATE(lNmed_ranmst)
  IF (findlamstar) DEALLOCATE(lstar_ranmst)
  IF (findgam) DEALLOCATE(lgam_ranmst)
  IF (findlamln) DEALLOCATE(lln_ranmst)

END SUBROUTINE find_lambda


SUBROUTINE calc_lambda(obj_avl,ran_avl,lam_val,lam_all,minus,plus,avranmstlen)
! A subroutine to calculate final lambda values; various different 
! measures of lambda (e.g. lambda bar, lambda tilde) can be passed 
! in to save the need for repeating this chunk of code every time 
! we want a new lambda.

  USE lambdaparams_module

  implicit none
  DOUBLE PRECISION, intent(in) :: obj_avl ! Average edge length (obj)
  DOUBLE PRECISION, dimension(1:nloop), intent(in) :: ran_avl ! " (ran)
  double precision, intent(out) :: lam_val,minus,plus
  double precision, intent(out) :: avranmstlen ! Average total length of
                                               ! random msts
  double precision, dimension(1:nloop) :: lam_all ! Lambda for each random MST
  integer :: i
  INTEGER, DIMENSION(1:nloop) :: listID        ! IDs of ran_avl list
  REAL :: ranup, ranlow          ! Upper & lower boundaries, 1 sigma
  
  do i = 1,nloop
     listID(i) = i
  END DO

!Sort nloop random MSTs in order of average edge length:
  CALL heapsort(nloop,ran_avl,listID)

! Find lambda for all random MSTs for lambda CDF plots:
  lam_all(1:nloop)=ran_avl(listID(1:nloop))/obj_avl

!Median of the average MST edge length:
  avranmstlen = ran_avl(listID(NINT(REAL(nloop)/2.)))
!1/6 and 5/6 boundaries for significance:
  ranlow = ran_avl(listID(NINT(REAL(nloop)/6.)))
  ranup = ran_avl(listID(NINT(5.*REAL(nloop)/6.)))

  lam_val = avranmstlen/obj_avl
  minus = ranlow/obj_avl
  plus = ranup/obj_avl

! check lambdas are ordered and match median value using this method:  
!  print*,
!  print*,lam_all(1)
!  print*,lam_all(2)
!  print*,lam_all(3)
!  print*,lam_all(4)
!  print*,lam_all(5)
!  print*,'  :'
!  print*,'  : median mst: ',ran_avl(listID(NINT(REAL(nloop)/2.)))*9.
!  print*,'  : object mst: ',obj_avl*9.
!  print*,'  : median lambda: ',ran_avl(listID(NINT(REAL(nloop)/2.)))/obj_avl
!  print*,'  : median lambda from lam_all: ',lam_all(nint(real(nloop/2.)))
!  print*,'  :'
!  print*,lam_all(nloop-4)
!  print*,lam_all(nloop-3)
!  print*,lam_all(nloop-2)
!  print*,lam_all(nloop-1)
!  print*,lam_all(nloop)

end subroutine calc_lambda

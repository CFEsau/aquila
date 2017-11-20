! Read in positions and masses for Aquila stars and calculate lambda
!(debug flags: https://stackoverflow.com/questions/29823099/debugging-with-gdb-and-gfortran-fpes)

program aquila_mst
  use starparams_module
  implicit none

  character (len=30) :: infile(3)       ! input file name for each cluster type
  integer :: corei, corenum             ! counter and number in catalogue
  character (len=15) :: corename        ! core name (15 chars - HGBS_J*)
  character (len=12) :: coretype,sample ! prestellar/starless/protostellar
  integer :: ncores
  logical :: usecores(3)      ! core types to be used in lambda calculations
  integer :: i,j,p            ! generic counters
  character(len=20) :: cdum   ! read non-needed strings
  integer :: io               ! for iostat specifier (end of file)
  LOGICAL :: dirExists        ! does given directory exist?

  infile(1)='../data/starless.dat'
  infile(2)='../data/prestellar.dat'
  infile(3)='../data/protostellar.dat'

  sample="starless" ! Currently all, starless, prestellar, or custom
  
  ! create directory for lambda data:
  outdir = '../cores_'//trim(sample)
  INQUIRE(file = outdir, exist = dirExists)
   IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(outdir))
  END IF
  
  !for starless use starless & prestellar. For prestellar just use 'prestellar'
  !Also set up directory for outputs.
  if (trim(sample)=='custom') then
     usecores(1)=.FALSE.  ! starless
     usecores(2)=.FALSE.  ! prestellar
     usecores(3)=.FALSE.  ! protostellar
  else if (trim(sample)=='all') then
     usecores(1)=.TRUE.   ! starless
     usecores(2)=.TRUE.   ! prestellar
     usecores(3)=.TRUE.   ! protostellar
  else if (trim(sample)=='starless') then
     usecores(1)=.TRUE.   ! starless
     usecores(2)=.TRUE.   ! prestellar
     usecores(3)=.FALSE.   ! protostellar
  else if (trim(sample)=='prestellar') then
     usecores(1)=.FALSE.   ! starless
     usecores(2)=.TRUE.   ! prestellar
     usecores(3)=.FALSE.   ! protostellar
  else
     write(6,*) "Core type not recognised."; stop
  end if
  
  
  ncores=0
  do i=1,3
     if(usecores(i)) then       ! If we want this core, open end of file
        open(unit=i,file=trim(infile(i)),status='old',position='append')
        backspace(i) !backspace to beginning of last line in file
        read(i,*) j
        rewind(i) !rewind to beginning of file for read-in below
        ncores=ncores+j !add number of cores in this file to total number
     end if
  end do
  
  
  allocate(RA(1:ncores))
  allocate(dec(1:ncores))
  allocate(rcore(1:ncores))
  allocate(robs(1:ncores))
  allocate(m(1:ncores))
  allocate(merr(1:ncores))
  allocate(T(1:ncores))
  allocate(Terr(1:ncores))
  allocate(ncolPeak(1:ncores))
  allocate(ncolCore(1:ncores))
  allocate(ncolObs(1:ncores))
  allocate(nvolPeak(1:ncores))
  allocate(nvolCore(1:ncores))
  allocate(nvolObs(1:ncores))
  allocate(mBE(1:ncores))
 
  
  ! Read in data:
  j=1 ! start with j=1 for first entry. Now j is array element, not counter
  do i=1,3
     if(usecores(i)) then
        do
           read(i,*,iostat=io) corei, corenum, corename, RA(j), dec(j), &
                & rcore(j), robs(j), m(j), merr(j), T(j), Terr(j), &
                & ncolPeak(j), ncolCore(j), ncolObs(j), &
                & nvolPeak(j), nvolCore(j), nvolObs(j), mBE(j), coretype
           if (io < 0) goto 100
           j=j+1  ! Successful read: increment j
        end do
100     close(i)
     end if
  end do
  
  

  paramdir = trim(outdir)//'/mass'
  !create output directory
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  !see whether two masses are equal. If so, shift one
  call shift(ncores,m)
  !and find lambda
  call lambda_setup(ncores,m)

  
  ! And do for all other paramaters:
  
  !create directory for parameter
  paramdir = trim(outdir)//'/rcore'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,rcore)
  call lambda_setup(ncores,rcore)

  
  paramdir = trim(outdir)//'/robs'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,robs)
  call lambda_setup(ncores,robs)

 
  paramdir = trim(outdir)//'/T'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,T)
  call lambda_setup(ncores,T)

  
  print*,ncores,maxval(ncolPeak)
  paramdir = trim(outdir)//'/ncolPeak'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,ncolPeak)
  call lambda_setup(ncores,ncolPeak)
  
  
  paramdir = trim(outdir)//'/ncolCore'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,ncolCore)
  call lambda_setup(ncores,ncolCore)

  
  paramdir = trim(outdir)//'/ncolObs'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,ncolObs)
  call lambda_setup(ncores,ncolObs)

  
  paramdir = trim(outdir)//'/nvolPeak'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,nvolPeak)
  call lambda_setup(ncores,nvolPeak)

  
  paramdir = trim(outdir)//'/nvolCore'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,nvolCore)
  call lambda_setup(ncores,nvolCore)

  
  paramdir = trim(outdir)//'/nvolObs'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,nvolObs)
  call lambda_setup(ncores,nvolObs)

  
  paramdir = trim(outdir)//'/mBE'
  INQUIRE(file = trim(paramdir), exist = dirExists)
  IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(paramdir))
  END IF
  call shift(ncores,mBE)
  call lambda_setup(ncores,mBE)
  
  
  deallocate(RA)
  deallocate(dec)
  deallocate(rcore)
  deallocate(robs)
  deallocate(m)
  deallocate(merr)
  deallocate(T)
  deallocate(Terr)
  deallocate(ncolPeak)
  deallocate(ncolCore)
  deallocate(nvolPeak)
  deallocate(nvolCore)
  deallocate(mBE)

end program aquila_mst





subroutine shift(ncores,param)

  implicit none
  integer, intent(in) :: ncores
  double precision, intent(inout) :: param(1:ncores)
  integer :: i,j,p

  
  ! Want to know if any stars have equal masses - that will mess up ordering.
  ! To avoid any biases in the catalogue, add a small random number to the mass
  ! of stars with equal masses so that any sorting bias in the catalogue is
  ! erased and the stars are sorted reliably.
  do
     p = 0
     do j = 1, ncores - 1
        do i = j+1, ncores
           if (param(i) == param(j)) then
              p=p+1
              !write(6,*) "Two values equal", j, i, param(j), param(i), p
              param(i) = param(i) + rand(0)*10.e-5
              !write(6,*) "Shouldn't be now...", j, i, param(j), param(i)
           end if
        end do
     end do
     if (p==0) exit
  end do
  
end subroutine shift

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
  character(len=20) :: cdum   ! read un-needed strings
  integer :: io               ! for iostat specifier (end of file)
  LOGICAL :: dirExists        ! does given directory exist?

  infile(1)='../data/starless.dat'
  infile(2)='../data/prestellar.dat'
  infile(3)='../data/protostellar.dat'

  sample="all" ! Currently all, starless, prestellar, or custom
  
  !If 'starless', use starless & prestellar.
  !If 'prestellar', just use prestellar
  
  !Set up output directory
  outdir = '../cores_'//trim(sample)
  INQUIRE(file = outdir, exist = dirExists)
   IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//trim(outdir))
  END IF
  
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
  
  
  j=0
  do i=1,3
     if(usecores(i)) then       ! If we want this core, open file
        open(i,file=trim(infile(i)),status='old')
        !count number of cores to allocate arrays:
        do
           read(i,*,iostat=io) cdum
           if (io < 0) goto 100
           if (len_trim(cdum)==0) then
              cycle
           else
              j=j+1  !Successful read: increment j
           end if
        end do
100     rewind(i)
     end if
  end do
  ncores=j

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
           if (io < 0) goto 200
           j=j+1  ! Successful read: increment j
        end do
200     close(i)
     end if
  end do

  ! Want to know if any stars have equal masses - that will mess up ordering.
  ! To avoid any biases in the catalogue, add a small random number to the mass
  ! of stars with equal masses so that any sorting bias in the catalogue is
  ! erased and the stars are sorted reliably.
  do
     p = 0
     do j = 1, ncores - 1
        do i = j+1, ncores
           !if (m(i) == m(j)) then
           if (T(i) == T(j)) then
              p=p+1
              !write(6,*) "Two masses equal", j, i, m(j), m(i), p
              !m(i) = m(i) + rand(0)*10.e-5
              !write(6,*) "Shouldn't be now...", j, i, m(j), m(i)
              !write(6,*) "Two temperatures equal", j, i, T(j), T(i), p
              T(i) = T(i) + rand(0)*10.e-5
              !write(6,*) "Shouldn't be now...", j, i, T(j), T(i)
              !write (6,*) ""
           end if
        end do
     end do
     if (p==0) exit
  end do
  
  !rcore, robs, m, T, ncolPeak, ncolCore, ncolObs,
  !nvolPeak, nvolCore, nvolObs, or mBE as input
  call lambda_setup(ncores,T)

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

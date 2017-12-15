!Read all AquilaM2 data and write selected data to output files
!Columns from original data:
! 1-2:  Core running number, Core name
! 3-4:  Core RA, Core dec
! 5-6:  Deconvolved radius, Observed radius
! 7-8:  Core mass (Msun), Uncertainty in m
! 9-10: Tdust, Uncertainty in Tdust
! 11:   Peak H2 column density
! 12:   Average column density of core using observed radius
! 13:   Average column density of core using deconvolved radius
! 14:   Peak volume density
! 15:   Average volume density using observed radius
! 16:   Average volume density using deconvolved radius
! 17:   Bonnor-Ebert mass ratio of core M_BE (crit/Mobs)
! 18:   Core type (starless, prestellar, protostellar)
! 19:   Comments. To be removed for outputs.

program readdata
  
  implicit none

  ! File I/O
  character (len=100) :: infile,outfile1,outfile2,outfile3
  character :: char
  integer :: corenum
  character (len=15) :: corename      ! Core name (15 chars)
  ! Core Right ascension and declination (J2000)
  character (len=11) :: RAstr, decstr ! RA, dec all 11 chars
  real :: rcore, robs                 ! Deconvolved & observerd radius
  real :: mcore, merr                 ! Mass of core & uncertainty (Msun)
  real :: Tdust, Terr                 ! Dust temperature & uncertainty
  real :: ncolPeak, ncolObs, ncolCore ! peak & average column densities
  real :: nvolPeak, nvolObs, nvolCore ! peak & average volume densities
  real :: mBE                         ! Bonnor-Ebert mass
  character (len=20) :: coretype      ! prestellar, starless, etc.
  character (len=500) :: wholeline
  character (len=13) :: ignore
  logical :: commentInLine
  real :: RAdeg,decdeg
! counters
  integer :: i, nstarless, nprestellar, nproto, n_all,nskip
  integer :: ofile,corei  !output file index and core number for output
  
  infile = 'HGBS_aquilaM2_derived_core_catalog.txt'
  outfile1 = 'starless.dat'
  outfile2 = 'prestellar.dat'
  outfile3 = 'protostellar.dat'
  write(6,*) ""
  write(6,*) "Input file: ",infile
  write(6,*) ""
  ignore='CO high-V_LSR' !ignore any objects that have this comment

  !skip past header info
  open(1,file=infile,status='old')
  read(1,*)char
  do while (char=='|')
     read(1,*)char
     !print *, char
  end do
  
  backspace(1)
  nstarless=0
  nprestellar=0
  nproto=0
  n_all=0
  nskip=0

  
  open(10,file=outfile1)
  open(11,file=outfile2)
  open(12,file=outfile3)
  
  do
     n_all=n_all+1
     read (1,'(a)',end=100) wholeline     !read line of data
     
     ! Check whether object has 'ignore' string in & cycle if so
     call hasstring(ignore,wholeline,commentInLine)
     if (commentInLine) then
        nskip=nskip+1
        !write (6,*) "SKIPPING ", corenum
        !print *, wholeline
        cycle ! go to next iteration
     end if
     
     backspace(1)      ! Go back to previous line
     read (1,*) corenum, corename, RAstr, decstr, rcore, robs, &
          & mcore, merr, Tdust, Terr, ncolPeak, ncolObs, ncolCore, &
          &  nvolPeak, nvolObs, nvolCore, mBE, coretype

     ! Find core type for output file direction
     if (coretype=='starless') then
        nstarless=nstarless+1
        ofile=10
        corei=nstarless
     else if (coretype=='prestellar') then
        nprestellar=nprestellar+1
        ofile=11
        corei=nprestellar
     else if (coretype=='protostellar') then
        nproto=nproto+1
        ofile=12
        corei=nproto
     end if
     
     ! convert from deg:arcmin:arcsec to decimal degrees:
     call degrees(RAstr,decstr,RAdeg,decdeg)
     write(ofile,20) corei, corenum, corename, RAdeg, decdeg, rcore, robs, &
             & mcore, merr, Tdust, Terr, ncolPeak, ncolObs, ncolCore, &
             & nvolPeak, nvolObs, nvolCore, mBE, coretype
  end do
  
  ! Character formatting: ID, ncore, name, RA, dec, rc, ro
20 format(2(2X,I4),2X,A,2(2X,F9.5),2(2X,(1P E7.2E1)) &
        & 0P, 2X,F6.2,2X,F5.2,2X,F4.1,2X,F3.1, & !Mc, Merr, T, Terr
        & 6(2X,F5.1),2X,F4.1,2X,A) !N, Nc, No, n, nc, no, mBE, type
 
100 CONTINUE
  n_all=n_all-1
  
  write(6,*) "Number of objects: ", n_all
  write(6,*) "Number of starless cores: ", nstarless
  write(6,*) "Number of prestellar cores: ", nprestellar
  write(6,*) "Number of protostellar cores:", nproto
  write(6,*) "Number of useable cores: ", nstarless+nprestellar+nproto,&
       & "    (Skipped",nskip,")"
  write(6,*) ""
  
  close(10)
  close(11)
  close(12)
  if (.not. nstarless+nprestellar+nproto+nskip==n_all) then
     write(6,*) "!!!!!!!!"
     write(6,*) "WARNING: Sum of objects doesn't equal number in file"
     write(6,*) "!!!!!!!!"
  end if
  
end program readdata



SUBROUTINE HASSTRING(A,B,AinB)	!Text B appears somewhere in text A?
  CHARACTER*(*) :: A,B
  LOGICAL :: AinB
  INTEGER :: L

  AinB = .FALSE.
  L = INDEX(B,A)		!The first position in A where B matches.
  IF (L.LE.0) THEN
     AinB = .FALSE.
  ELSE
     AinB = .TRUE.
  END IF
  !write(6,*) "String A: ", A,  " String B: ", B
  !write(6, *) AinB, index(B,A)
  RETURN

END SUBROUTINE HASSTRING



subroutine degrees(RAstr,decstr,radeg,decdeg)
  character(len=11), intent(in) :: RAstr, decstr
  real, intent(out) :: RAdeg, decdeg
  real :: deg, min, sec


  !RA always positive, dec always negative here:
  !RA: hh:mm:ss.ss (colons at index 3,6)
  read(RAstr(1:2),*) deg
  read(RAstr(4:5),*) min
  read(RAstr(7:11),*) sec
  !print *, deg, min, sec
  radeg = deg + min/60 + sec/3600
  !print *,radeg
  
  !decination: -hh:mm:ss.s (colons at index 4,7)
  read(decstr(1:3),*) deg
  read(decstr(5:6),*) min
  read(decstr(8:11),*) sec
  !print *, deg, min, sec
  decdeg = deg - min/60 - sec/3600
  !print *,decdeg
  return
  
end subroutine degrees

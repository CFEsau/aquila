!******************************************************************************!
!
! Modules for analysis code
!
!******************************************************************************!

MODULE directories_module
  character(len=100) :: paramdir  ! Destination directory (e.g. 'outputs/m')
  character(len=100) :: lampath, CDFpath
  !RA & dec in decimal format
  double precision, dimension(:), allocatable :: RA, dec
END MODULE directories_module


MODULE lambdaparams_module
  
  integer :: nmst            ! number of stars in MST
  INTEGER :: nloop         ! number of random MSTs calculated in loop
  
! Lengths of 'object' edges in MST for CDFs
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: edgelengths
! Number of median values to include in l_Nmed (e.g. 2 or 3 for even/odd nedge)
  INTEGER :: Nmed
  INTEGER :: nCDF       ! Number of CDF plots of random MSTs
! lambda = measure of mass segregation & errors
  DOUBLE PRECISION :: lambda, l_up, l_low
! mean MST length of nloop random MSTs, object MST length
  DOUBLE PRECISION :: l_avranmst, l_objmst
  double precision, dimension(:), allocatable :: l_allmsts
! lambda_bar uses mean length of MST (basically same as lambda)
  DOUBLE PRECISION :: lambda_bar, l_up_bar, l_low_bar
  DOUBLE PRECISION :: lbar_avranmst, lbar_objmst
  double precision, dimension(:), allocatable :: lbar_allmsts
! lambda_rms uses root mean square length of MST (generalised mean with p=2)
  DOUBLE PRECISION :: lambda_rms, l_up_rms, l_low_rms
  DOUBLE PRECISION :: lrms_avranmst, lrms_objmst
  double precision, dimension(:), allocatable :: lrms_allmsts
! lambda_smr uses square mean root length of MST (generalised mean with p=1/2)
  DOUBLE PRECISION :: lambda_smr, l_up_smr, l_low_smr
  DOUBLE PRECISION :: lsmr_avranmst, lsmr_objmst
  double precision, dimension(:), allocatable :: lsmr_allmsts
! lambda_har uses harmonic mean (generalised mean with p=-1)
  DOUBLE PRECISION :: lambda_har, l_up_har, l_low_har
  DOUBLE PRECISION :: lhar_avranmst, lhar_objmst
  double precision, dimension(:), allocatable :: lhar_allmsts
! lambda_tilde uses median length of MST
  DOUBLE PRECISION :: lambda_til, l_up_til, l_low_til
  DOUBLE PRECISION :: ltil_avranmst, ltil_objmst
  double precision, dimension(:), allocatable :: ltil_allmsts
! lambda_Nmed uses the mean of the N median lengths of MST (Nmed)
  DOUBLE PRECISION :: lambda_Nmed, l_up_Nmed, l_low_Nmed
  DOUBLE PRECISION :: lNmed_avranmst, lNmed_objmst
  double precision, dimension(:), allocatable :: lNmed_allmsts
! lambda_star uses median length of MST and adds this to the actual length of the MST
  DOUBLE PRECISION :: lambda_star, l_up_star, l_low_star
  DOUBLE PRECISION :: lstar_avranmst, lstar_objmst
  double precision, dimension(:), allocatable :: lstar_allmsts
! lambda_gam uses the geometric mean
  DOUBLE PRECISION :: lambda_gam, l_up_gam, l_low_gam
  DOUBLE PRECISION :: lgam_avranmst, lgam_objmst
  double precision, dimension(:), allocatable :: lgam_allmsts
! lamda_ln: sum the exponents and take the log (sort of inverse of geometric mean)
  DOUBLE PRECISION :: lambda_ln, l_up_ln, l_low_ln
  DOUBLE PRECISION :: lln_avranmst, lln_objmst
  double precision, dimension(:), allocatable :: lln_allmsts
  
! state which types of lambda you want to find, to save doing all every time
  LOGICAL :: findlam, findlambar, findlamrms, findlamsmr, findlamhar
  LOGICAL :: findlamtil, findlamNmed, findlamstar, findgam, findlamln
  
END MODULE lambdaparams_module

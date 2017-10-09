Aquila study

Data from http://gouldbelt-herschel.cea.fr/archives
File HGBS_aquilaM2_derived_core_catalog.txt

Mass segregation checks: starless (unbound+bound)  vs  prestellar (bound) subsamples.

From Vera:
To get the STARLESS cores for the tests you may choose the "starless" AND "prestellar"
cores together based on their Core_type (col. 18), excluding the ones with "CO high-V_LSR"
comments (col. 19). This should give 651 unbound/bound starless cores.

For the PRESTELLAR cores, please use the subsample of "prestellar" cores (Core_type, col.
18), excluding the "CO high-V_LSR" commented ones (col. 19), I think there's 446 of them.

Again, what I could see with Allison+2009's Lambda_MST script, is that among the starless
sample there is positive mass segregation for the highest mass (N=10-20-50) cores, AND also
for the smallest mass cores, while these latters are not part of the bound prestellar sample,
so there only the highest-mass bins had >0 Lambda_MST.

This all shows up for the whole field, but maybe even more significantly for the central W40 -
Sepens South filaments region. I presented both the part vs whole in my talk. In order to get a
subsample only around W40 - Sepens S, I basically selected the cores in that part of the map
above an A_V (visual extinction) of ~ 7 mag from the corresponding column density map (also
available in the same line of the archive: HGBS_aquilaM2_column_density_map.fits.gz).

-------------------------------------------------------------

'data' directory:

Data file read in and written to alternative files using readdata.f90

Columns in prestellar.dat and starless.dat:
ID    core#    corename    RA(deg)    dec(deg)    Rdeconvolved    Robs
mass    merr    Tdust    Terr    npeak(cm^-2)    n(cm^-2)    nobs(cm^-2)
Npeak(cm^-3)    N(cm^-3)    Nobs(cm^-3)    M_BE    coretype

Outputs:
prestellar.dat
starless.dat

--------------------------------------------------------------

'src' directory:

runcode.txt: compile code/modules/subroutines & run findlambda

aquila_mst.f90
Input file: .dat files from 'data'
- Set which core types are required (starless, prestellar, protostellar)
- Cycle through core type files and count objects to allocate arrays
- Allocate variables desired for MSTs
- Cycle through core type files again reading in data
- Call lambda_setup to allocate, zero, and deallocate arrays used in lambda calculations

  lambda_setup:
  - Set number of msts
  - Set up CDF & lambda file units
  - Allocate array for MST edge lengths
  - Set types of lambda to calculate
  - Create CDF directory if needed
  - Open MSTedgeL.dat - MST lengths for object stars (1 row of nedge entries, as one snapshot)
  - Open STedgeL_#.dat - MST lengths for random stars, where suffix is one of nCDF
  - call find_lambda to calculate various lambda values
  - Close CDF files
  - Open lambda data files (MST_[lamtype].dat) and write MST lengths for object & random MSTs,
    lambda values, and lambda errors (again one row per file, as one snapshot)
  - Open file and write lambda values for each random MST in allMSTs_[lamtype].dat
  - Deallocate lambda arrays

    find_lambda:
    - 


Summary of outputs:
lambda directory, MST_[type].dat for obj & medianrandom MST lengths, lambda values, & errors
CDFdata directory, allMSTs_[type].dat for lambda calculated using all nloop MSTs
		   MSTedgeL_[#] for edge lengths in various MSTs (one MST per file, 1 snapshot)

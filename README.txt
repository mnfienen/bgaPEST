bgaPEST -- Version 1.0
Mike Fienen (mnfienen@usgs.gov)
and
Marco D'Oria (marco.doria@unipr.it)
November 19, 2012

This folder contains the code, compiled executables, and examples to accompany:
Approaches in Highly Parameterized Inversion: bgaPEST, a Bayesian Geostatistical Approach 
Implementation With PEST—Documentation and Instructions
By Michael N. Fienen, Marco D’Oria, John E. Doherty, and Randall J. Hunt
Techniques and Methods, Book 7, Section C9
Published 2012
http://pubs.usgs.gov/tm/tm7c9

GENERAL INSTRUCTIONS:
1) Create an input file, following the documentation, with the .bgp extension
2) Copy bgaPEST.exe to the directory in which the problem is to run
2a) If using derivmode=0, also copy pest.exe and pst_generator.exe into the run directory
3) Type bgaPEST.exe <casename>.bgp at the command line to run


The subdirectories in this directory contain the following:
1_EXECUTABLES
This folder has subdirectory containing a 32-bit and 64-bit version of bgaPEST.exe, compiled
for Windows using the Intel Visual Fortran Compiler in Microsoft Visual Studio.

2_EXAMPLES
This subdirectory contains the input and executable files needed to run the examples presented
in the report cited above. Furthermore, an example using the parallel facility that uses
HTCondor is included and a fully 3-D version of the 3 layer problem is also included.

3_SOURCE
This subdirectory includes a file bgaPEST.sln for use as a solution file with Visual Studio 10.
The src subdirectory contains all Fortran source code used for bgaPEST. The inschek subdirectory 
is a double-precision version of the PEST utility inschek for use with the parallel option of
bgaPEST. Finally, in the src subdirectory is a python subdirectory containing the source code 
for pst_generator and the parallel CHTCondor utilities.

DISCLAIMER and NOTICE
Please refer to the USGS Software User Rights Notice (http://water.usgs.gov/software/help/notice/) for complete use, copyright, and distribution information. The USGS provides no warranty, expressed or implied, as to the correctness of the furnished software or the suitability for any purpose. The software has been tested, but as with any complex software, there could be undetected errors. Users who find errors are requested to report them to the USGS.

References to non-USGS products, trade names, and (or) services are provided for information purposes only and do not constitute endorsement or warranty, express or implied, by the USGS, U.S. Department of Interior, or U.S. Government, as to their suitability, content, usefulness, functioning, completeness, or accuracy.
# Info (EOF MultiRegion (MR) Modification)

Modification of the EOF Library for the use of with multiple regions

Compilation of the solver:


The following changes have to be made bevor using this:

modifiy - SOLVER.KEYWORDS in /opt/elmerfem/elmer-mpi-install/share/elmersolver/lib 


addiational entries are: 

$ for(i=1:nexp) "Solver:Integer: 'Body " i2str(i) " Use Target Variable'"
Solver:Integer:     'Bodies'


There seems to be some identical fuctionallity also introduced in development 
branch of the EOF library, but I have not tested it by now.

https://github.com/jvencels/EOF-Library/tree/devel

(master branch commit 82f3bdb3c509f5bd243671ef3e89bbc172557506, 27.05.2019)
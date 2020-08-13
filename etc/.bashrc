PATH_OLD=$PATH
LD_LIBRARY_PATH_OLD=$LD_LIBRARY_PATH

# --- Modify paths:
MYPATH="/home/user/OpenFOAM/user-v1912/applications"
FOAM_PATH="/opt/openfoam/OpenFOAM-v1912"
ELMER_PATH="/opt/elmerfem/elmer-mpi-install"
# --- End modify

EOF_OF_SRC="$MYPATH/libs/eof/coupleElmer"
EOF_SOLVER_SRC="$MYPATH/libs/eof/elmer_solver"
EOF_MR_SOLVER_SRC="$MYPATH/libs/elmer_mr_solver"

source_of1912="source $FOAM_PATH/bashrc"
path_elmerfem="export PATH=$ELMER_PATH/bin:\$PATH"

ld_user="export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH"
ld_foam_user="export LD_LIBRARY_PATH=\$FOAM_USER_LIBBIN:\$LD_LIBRARY_PATH"
set_eofsrc="export EOF_SRC=$EOF_SOLVER_SRC"
set_eofsrc_mrmod="export EOF_SRC=$EOF_MR_SOLVER_SRC"
ld_eof="export LD_LIBRARY_PATH=\$EOF_OF_SRC:\$LD_LIBRARY_PATH"

reset_path="export PATH=$PATH_OLD"
reset_ld="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_OLD"

alias elmerfem="$source_of1912;$path_elmerfem;$ld_user "
alias of1912="$source_of1912"
alias eof1912="elmerfem;$set_eofsrc;$ld_eof$"
alias eof1912_mrmod="elmerfem;$set_eofsrc_mrmod;$ld_eof"

eofCompile() {
    wmake $EOF_OF_SRC/coupleElmer
    elmerf90 -o $EOF_SRC/Elmer2OpenFOAM.so -J$EOF_SRC $EOF_SRC/Elmer2OpenFOAM.F90
    elmerf90 -o $EOF_SRC/OpenFOAM2Elmer.so -J$EOF_SRC $EOF_SRC/OpenFOAM2Elmer.F90
}




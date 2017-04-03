#!/bin/bash -x
#
# Script: SvarDiag.sh
# Purpose: Produce aux images to help diagnose the SVAR spikes 
# Author: T. Nichols & S. Afyouni
# Version: $Id$
#


###############################################################################
#
# Environment set up
#
###############################################################################

shopt -s nullglob # No-match globbing expands to null
TmpDir=/tmp
Tmp=$TmpDir/`basename $0`-${$}-
trap CleanUp INT
TMPOUTPUTTYPE="$FSLOUTPUTTYPE"
export FSLOUTPUTTYPE="NIFTI"


###############################################################################
#
# Functions
#
###############################################################################

Usage() {
cat <<EOF
Usage: `basename $0` [options] fMRI_4Dtimeseries BND1 BND2 DVARSout
INPUTS:
fMRI_4Dtimeseries : 4D Resting-state image
BND1              : Upper bound where the Svar should start from
BND2              : Lowe bound where the Svar image should stop
DVARSout          : Name and directory for the output image (Svar image) 

SA & TEN, UoW, 2017
_________________________________________________________________________
\$Id$
EOF
exit
}

CleanUp () {
    /bin/rm -f ${Tmp}*
    exit 0
}


###############################################################################
#
# Parse arguments
#
###############################################################################
while (( $# > 1 )) ; do
    case "$1" in
        "-help")
            Usage
            ;;
        -*)
            echo "ERROR: Unknown option '$1'"
            exit 1
            break
            ;;
        *)
            break
            ;;
    esac
done

if (( $# < 2 )) ; then
    Usage
fi

FUNC="$1"
BND1="$2"
BND2="$3"
OUT="$4"
###############################################################################
#
# Script Body
#
###############################################################################
echo -n "."

Nvol=$(fslnvols "$FUNC")

# Find mean over time
fslmaths "$FUNC" -Tmean $Tmp-Mean
echo $Tmp-Mean
# Find the brain
bet $Tmp-Mean  $Tmp-MeanBrain

# Demean
fslmaths "$FUNC" -sub $Tmp-Mean -mas $Tmp-MeanBrain $Tmp-Demean

#std
fslmaths $Tmp-Demean -Tstd $OUT-DemeanStd

BNDR=$(($BND2 - $BND1))
echo "Designated bounds are length of : ${BNDR}"
echo "  From ${BND1} to ${BND2}"

fslroi $Tmp-Demean $Tmp-BND0  $((BND1 - 1))     "$BNDR"
fslroi $Tmp-Demean $Tmp-BND1  "$BND1"           "$BNDR"

fslmaths $Tmp-BND0 -add $Tmp-BND1 -div 2 -sqr $OUT-Svar


###############################################################################
#
# Exit & Clean up
#
###############################################################################

CleanUp



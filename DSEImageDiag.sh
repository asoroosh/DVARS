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
export FSLOUTPUTTYPE="NIFTI_GZ" 
# you might wanna change this to NIFTI, but note on 7T resting-state data, a single file will exceed 5GB!


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
#BND1="$2"
#BND2="$3"
OUT="$2"

#Create a new directory, if doesn't already exists, as variable OUT
Dir2Save=`dirname "$OUT"`
PreFix=`basename "$OUT"`

mkdir -p $Dir2Save
echo "Created: $Dir2Save"

###############################################################################
#
# Script Body
#
###############################################################################
#echo -n "."

Nvol=$(fslnvols "$FUNC")
echo "The input image has $Nvol volumes."

# Find mean over time
fslmaths "$FUNC" -Tmean $Tmp-Mean
echo $Tmp-Mean

# Form a mask, later will be used to get the Dvar and Svar time series
fslmaths $Tmp-Mean -bin $Tmp-Mean-mask 

# Demean
#fslmaths "$FUNC" -sub $Tmp-Mean -mas $Tmp-MeanBrain $Tmp-Demean
fslmaths "$FUNC" -sub $Tmp-Mean $Tmp-Demean

#std
fslmaths $Tmp-Demean -Tstd $Dir2Save/$PreFix-DemeanStd

fslroi $Tmp-Demean $Tmp-BND0  0     $((Nvol - 1))
fslroi $Tmp-Demean $Tmp-BND1  1     "$Nvol"

#Sanity check
#fslinfo $Tmp-BND0
#fslinfo $Tmp-BND1

echo "Generating Svar and Dvar 4D data..."
fslmaths $Tmp-Demean -sqr $Dir2Save/$PreFix-Avar
fslmaths $Tmp-BND0   -add $Tmp-BND1 -div 2 -sqr $Dir2Save/$PreFix-Svar
fslmaths $Tmp-BND0   -sub $Tmp-BND1 -div 2 -sqr $Dir2Save/$PreFix-Dvar

fslmeants -i $Dir2Save/$PreFix-Avar -m $Tmp-Mean-mask -o $Dir2Save/$PreFix-Avar-meants.txt
fslmeants -i $Dir2Save/$PreFix-Svar -m $Tmp-Mean-mask -o $Dir2Save/$PreFix-Svar-meants.txt
fslmeants -i $Dir2Save/$PreFix-Dvar -m $Tmp-Mean-mask -o $Dir2Save/$PreFix-Dvar-meants.txt

echo "Generating Svar and Dvar 3D images..."
fslmaths $Dir2Save/$PreFix-Svar -Tmean $Dir2Save/$PreFix-mSvar
fslmaths $Dir2Save/$PreFix-Dvar -Tmean $Dir2Save/$PreFix-mDvar

echo "Done!"
###############################################################################
#
# Exit & Clean up
#
###############################################################################

CleanUp

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
Usage: `basename $0` [options] fMRI_4Dtimeseries OUT 4DFlag
==INPUTS:
fMRI_4Dtimeseries : 4D Resting-state image
OUT               : Directory and a prefix for the outputs, e.g. path/to/dir/SUB001
4DFlag [optional] : if set to 1, both 3D and 4D variability images are saved, if set to 0 then only 3D images are saved [default 1]. 

==OUTPUTS:
    4D images:
        Svar, Dvar & Avar. pSvar & pDvar. DeltapDvar & DeltapSvar
    3D images:
        Svar, Dvar & Avar. pSvar & pDvar. DeltapDvar & DeltapSvar
    TXT:
        mean Svar, mean Dvar & mean Avar. mean pSvar & mean pDvar. mean DeltapSvar & mean DeltapDvar

==DEPENDENCIES:
FSL should have been already installed. 

==REFERENCE:
Afyouni, Soroosh, and Thomas E. Nichols. "Insight and Inference for DVARS." NeuroImage (2018).

SA & TEN, Ox, 2018
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
OUT="$2"

if [ -z "$3" ]; 
then 
	echo "4D flag is unset so 4D images will be saved..."
	need4D=1 
else 
	need4D="$3"
	echo "4D flag is set to $need4D so 4D images will not be saved." 
fi

#Create a new directory, if doesn't already exists, as variable OUT
Dir2Save=`dirname "$OUT"`
PreFix=`basename "$OUT"`

mkdir -p $Dir2Save/DSE
mkdir -p $Dir2Save/pDSE
mkdir -p $Dir2Save/DeltapDSE

echo "Created: $Dir2Save"

###############################################################################
#
# Script Body
#
###############################################################################
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

echo "Generating Avar, Svar and Dvar 4D data + non-normalised time series."
fslmaths $Tmp-Demean -sqr $Dir2Save/DSE/$PreFix-Avar
fslmaths $Tmp-BND0   -add $Tmp-BND1 -div 2 -sqr $Dir2Save/DSE/$PreFix-Svar
fslmaths $Tmp-BND0   -sub $Tmp-BND1 -div 2 -sqr $Dir2Save/DSE/$PreFix-Dvar

fslmeants -i $Dir2Save/DSE/$PreFix-Avar -m $Tmp-Mean-mask -o $Dir2Save/DSE/$PreFix-Avar-meants.txt
fslmeants -i $Dir2Save/DSE/$PreFix-Svar -m $Tmp-Mean-mask -o $Dir2Save/DSE/$PreFix-Svar-meants.txt
fslmeants -i $Dir2Save/DSE/$PreFix-Dvar -m $Tmp-Mean-mask -o $Dir2Save/DSE/$PreFix-Dvar-meants.txt

awk '{ total += $1; count++ } END { print total/count }' $Dir2Save/DSE/$PreFix-Avar-meants.txt >  $Dir2Save/DSE/$PreFix-Whole-DSE.txt
awk '{ total += $1; count++ } END { print total/count }' $Dir2Save/DSE/$PreFix-Svar-meants.txt >> $Dir2Save/DSE/$PreFix-Whole-DSE.txt
awk '{ total += $1; count++ } END { print total/count }' $Dir2Save/DSE/$PreFix-Dvar-meants.txt >> $Dir2Save/DSE/$PreFix-Whole-DSE.txt

echo "Generating Avar 3D image..."
fslmaths $Dir2Save/DSE/$PreFix-Avar -Tmean $Dir2Save/DSE/$PreFix-mAvar

#Generate the median images for calculating the DeltapDvar and DeltapSvar...
fslmaths $Dir2Save/DSE/$PreFix-Svar -Tmedian $Dir2Save/$PreFix-mdSvar
fslmaths $Dir2Save/DSE/$PreFix-Dvar -Tmedian $Dir2Save/$PreFix-mdDvar

echo "Generating %Svar and %Dvar 4D images..."
#Normalise Dvar and Svar by Avar i.e. %Svar and %Dvar
fslmaths $Dir2Save/DSE/$PreFix-Svar -div $Dir2Save/DSE/$PreFix-mAvar $Dir2Save/pDSE/$PreFix-pSvar
fslmaths $Dir2Save/DSE/$PreFix-Dvar -div $Dir2Save/DSE/$PreFix-mAvar $Dir2Save/pDSE/$PreFix-pDvar

echo "Generating Delta%Svar and Delta%Dvar 4D images..."
fslmaths $Dir2Save/DSE/$PreFix-Svar -sub $Dir2Save/$PreFix-mdSvar -div $Dir2Save/DSE/$PreFix-mAvar $Dir2Save/DeltapDSE/$PreFix-DeltapSvar
fslmaths $Dir2Save/DSE/$PreFix-Dvar -sub $Dir2Save/$PreFix-mdDvar -div $Dir2Save/DSE/$PreFix-mAvar $Dir2Save/DeltapDSE/$PreFix-DeltapDvar

#remove the median images
rm $Dir2Save/$PreFix-mdSvar.nii.gz $Dir2Save/$PreFix-mdDvar.nii.gz

echo "Generating %Svar and %Dvar time series..."
fslmeants -i $Dir2Save/pDSE/$PreFix-pSvar -m $Tmp-Mean-mask -o $Dir2Save/pDSE/$PreFix-pSvar-meants.txt
fslmeants -i $Dir2Save/pDSE/$PreFix-pDvar -m $Tmp-Mean-mask -o $Dir2Save/pDSE/$PreFix-pDvar-meants.txt

awk '{ total += $1; count++ } END { print total/count }' $Dir2Save/pDSE/$PreFix-pSvar-meants.txt >  $Dir2Save/pDSE/$PreFix-Whole-p-DSE.txt
awk '{ total += $1; count++ } END { print total/count }' $Dir2Save/pDSE/$PreFix-pDvar-meants.txt >> $Dir2Save/pDSE/$PreFix-Whole-p-DSE.txt

echo "Generating Delta%Svar and Delta%Dvar time series..."
fslmeants -i $Dir2Save/DeltapDSE/$PreFix-DeltapSvar -m $Tmp-Mean-mask -o $Dir2Save/DeltapDSE/$PreFix-DeltapSvar-meants.txt
fslmeants -i $Dir2Save/DeltapDSE/$PreFix-DeltapDvar -m $Tmp-Mean-mask -o $Dir2Save/DeltapDSE/$PreFix-DeltapDvar-meants.txt

echo "Generating %Svar and %Dvar 3D images..."
fslmaths $Dir2Save/pDSE/$PreFix-pSvar -Tmean $Dir2Save/pDSE/$PreFix-mpSvar
fslmaths $Dir2Save/pDSE/$PreFix-pDvar -Tmean $Dir2Save/pDSE/$PreFix-mpDvar

echo "Generating Delta%Svar and Delta%Dvar 3D images..."
fslmaths $Dir2Save/DeltapDSE/$PreFix-DeltapSvar -Tmean $Dir2Save/DeltapDSE/$PreFix-mDeltapSvar
fslmaths $Dir2Save/DeltapDSE/$PreFix-DeltapDvar -Tmean $Dir2Save/DeltapDSE/$PreFix-mDeltapDvar

if [ $need4D == 0 ]
then
	echo "Now deleting 4D images to free up some space..."
        # remove the raw Svar and Dvar
        rm $Dir2Save/DSE/$PreFix-Dvar.nii.gz $Dir2Save/DSE/$PreFix-Svar.nii.gz $Dir2Save/DSE/$PreFix-Avar.nii.gz
        # remove the pSvar and pDvar images
	    rm $Dir2Save/pDSE/$PreFix-pDvar.nii.gz $Dir2Save/pDSE/$PreFix-pSvar.nii.gz
        # remove DeltapSvar and DeltapDvar
        rm $Dir2Save/DeltapDSE/$PreFix-DeltapDvar.nii.gz $Dir2Save/DeltapDSE/$PreFix-DeltapSvar.nii.gz
fi

echo "Done!"
###############################################################################
#
# Exit & Clean up
#
###############################################################################

CleanUp

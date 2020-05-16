#bin/bash

set -e

origimage=$1   # original image
maskimage=$2   # mask of the oirginal image
smoothparMM=$3 # smoothing parameter in mm
finalresult=$4 # Name of the ouput parameter

smoothparSig=$(bc -l <<< "${smoothparMM}/2.3548")

finalresultdir=$(dirname $finalresult)

if [ -z "$maskimage" ] || [ "$maskimage" == " " ]; then
	echo "No mask was set, so I will make one. "
	maskimage=${finalresultdir}/mask_internal
	${FSLDIR}/bin/fslmaths $origimage -Tmean -bin ${maskimage}
fi

echo "===SMOOTHING:"
echo "-IN:   ${origimage}"
echo "-MASK: ${maskimage}"
echo "-OUT:  ${finalresult}"
echo "-KERNEL(mm): ${smoothparMM} , sigma: ${smoothparSig}"

${FSLDIR}/bin/fslmaths $origimage -s $smoothparSig -mas $maskimage ${finalresultdir}/tmp_result1_tmp
${FSLDIR}/bin/fslmaths $maskimage -s $smoothparSig -mas $maskimage $finalresult
${FSLDIR}/bin/fslmaths ${finalresultdir}/tmp_result1_tmp -div $finalresult $finalresult

${FSLDIR}/bin/imrm ${finalresultdir}/tmp_result1_tmp

#bin/bash

set -e


funcdata=$1
smoothparMM=$3
mask=$2
smoothdfuncdata=$4

finalresultdir=$(dirname $smoothdfuncdata)

if [ -z "$mask" ] || [ "$mask" == " " ]; then
        echo "No mask was set, so I will make one. "
        mask=${finalresultdir}/mask_internal
        ${FSLDIR}/bin/fslmaths $funcdata -Tmean -bin ${mask}
fi

int2=$(${FSLDIR}/bin/fslstats $funcdata -p 2 -p 98 |awk {'print int($1)'})
int98=$(${FSLDIR}/bin/fslstats $funcdata -p 2 -p 98 |awk {'print $2'})

echo ${int2}

#line 5349, featlib.tcl
median_intensity=$(${FSLDIR}/bin/fslstats $funcdata -k $mask -p 50)
smoothsigma=$(bc -l <<< "${smoothparMM}/2.3548")
susan_int=$(bc -l <<< "($median_intensity-$int2) * 0.75")

echo "===fslsusan:"
echo "-IN:   ${funcdata}"
echo "-MASK: ${mask}"
echo "-OUT:  ${smoothdfuncdata}"
echo "fslsusan:: susan bt:${susan_int}, sigma:${smoothsigma}, median_intensity:${median_intensity},"

${FSLDIR}/bin/fslmaths $funcdata -Tmean $finalresultdir/mean_func_tmp
${FSLDIR}/bin/susan $funcdata $susan_int $smoothsigma 3 1 1 $finalresultdir/mean_func_tmp $susan_int $smoothdfuncdata
${FSLDIR}/bin/fslmaths $smoothdfuncdata -mas $mask $smoothdfuncdata

${FSLDIR}/bin/imrm mean_func_tmp

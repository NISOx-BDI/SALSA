#bin/bash

# SA, Ox, 2020

set -e

COPEDIR=$1
VCOPEDIR=$2

tmpdir=$(mktemp -d)

NUMS=$(${FSLDIR}/bin/fslinfo ${VCOPEDIR} | sed '5q;d' | awk '{print $2}')
echo "Number of subjects: $NUMS"

imgdir=$(dirname "${VCOPEDIR}")
#mkdir -p ${imgdir}/FLAME1 # falmeo makes one

## mat & grp -------------------------
grp=${tmpdir}/grp.txt
mat=${tmpdir}/mat.txt
cont=${tmpdir}/cont.mat

echo $grp
echo $mat
echo $cont

for (( ll=1; ll<=$NUMS; ll++ ))
do
	echo "1" >> ${grp}
	echo "1" >> ${mat}
done
## design.grp ------------------------
${FSLDIR}/bin/Text2Vest ${grp} ${tmpdir}/design.grp
${FSLDIR}/bin/Text2Vest ${mat} ${tmpdir}/design.mat

## cont.mat --------------------------
echo "/NumWaves       1" >> $cont
echo "/NumContrasts   1" >> $cont
echo "" >> $cont
echo "" >> $cont
echo "/Matrix" >> $cont
echo "1" >> $cont


## Make a mask:
${FSLDIR}/bin/fslmaths ${VCOPEDIR} -Tmean ${tmpdir}/mask

${FSLDIR}/bin/flameo --cope=${COPEDIR} \
--varcope=${VCOPEDIR} \
--runmode=flame1 \
--dm=${tmpdir}/design.mat \
--tc=${tmpdir}/cont.mat \
--logdir=${imgdir}/FLAME1 \
--mask=${tmpdir}/mask \
--cs=${tmpdir}/design.grp

echo "Results: ${imgdir}/FLAME1"
echo "DONE"

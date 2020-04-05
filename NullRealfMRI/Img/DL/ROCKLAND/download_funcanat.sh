#/bin/bash

module load Python/3.6.6-foss-2018b

CohortID=ROCKLAND
outp="/well/nichols/users/scf915/${CohortID}"

python3 download_rockland_raw_bids.py -o $outp -e REST645 -t func anat -v DS2


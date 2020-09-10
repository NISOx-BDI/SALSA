#/bin/bash

module load Python/3.6.6-foss-2018b

CohortID=ROCKLAND
outp="/well/nichols/users/scf915/${CohortID}/raw"

mkdir -p $outp

ImgTyp=CHECKERBOARD645 #REST1400 #RESTCAP #REST645

python3 download_rockland_raw_bids.py -o $outp -e ${ImgTyp} -t func -v DS2


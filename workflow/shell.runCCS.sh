#!/bin/bash
#source /ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Denovo/PacBio/SMRTAnalysis/smrtanalysis/install/smrtanalysis_2.3.0.140936/etc/setup.sh
source /home/zhouze/team/local/smrtanalysis/current/etc/setup.sh
fofnToSmrtpipeInput.py my_inputs.fofn > my_inputs.xml
smrtpipe.py -D TMP=./output/TMP -D SHARED_DIR=./output/SHARE --output=./output --params=./my_ReadOfInsert.xml  xml:my_inputs.xml

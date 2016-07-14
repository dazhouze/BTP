#!/bin/bash
source /ifs1/ST_IM/USER/zhouze/local/smrtanalysis/current/etc/setup.sh
referenceUploader -c -p/ifs1/ST_IM/USER/zhouze/dataBase/smrt -n ucsc_hg19 -f/ifs1/ST_IM/USER/database/Resequence/database/Ref/hg19/ucsc.hg19.fa --saw='sawriter -welter' --samIdx='samtools faidx'
referenceUploader -c -p/ifs1/ST_IM/USER/zhouze/dataBase/smrt -n hg38_chr6 -f/ifs1/ST_IM/USER/database/Resequence/database/Ref/hg38/ucsc.chr6.fa --saw='sawriter -welter' --samIdx='samtools faidx'

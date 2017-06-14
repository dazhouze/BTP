#!/bin/bash
dir=MHC
ref=/ifs1/ST_IM/USER/database/Resequence/database/Ref/hg19/ucsc.hg19.fa
ccs_bam=/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/merged5YH.best.ccs.sort.bam
ccs_fastq="/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.1.ccs.fastq  /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.2.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.3.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.1.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.2.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.3.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.1.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.2.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.3.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.1.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.2.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.3.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.1.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.2.ccs.fastq /ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.3.ccs.fastq"
python3 phasing.py -o $dir 
mkdir -p $dir/Pos
mkdir -p $dir/Fastq
mkdir -p $dir/Qname
mv $dir/*fastq $dir/Fastq
mv $dir/*txt $dir/Qname
python3 example/fastq_fetch.py $dir $ccs_fastq
python3 example/canu_fetch.py $dir $ccs_bam
(perl example/canu_assembly_0.pl $dir &); (perl example/canu_assembly_1.pl $dir &)

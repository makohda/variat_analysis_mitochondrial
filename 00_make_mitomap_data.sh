#!/bin/sh

current_dir=`pwd`
mitomap_dir="./MitoMap"
day=`date "+%y%m%d"`
outfile=MITOMAP.cfrm_report_all.${day}

if [ -e ${mitomap_dir} ]; then
  :
else
  mkdir ${mitomap_dir}
fi

cd ${mitomap_dir}

if [ -e ${day} ]; then
  cd ${day}
else
  mkdir ${day}
  cd ${day}
fi

wget -c https://www.mitomap.org/foswiki/bin/view/MITOMAP/MutationsRNA
#=> MutationsRNA
wget -c https://www.mitomap.org/foswiki/bin/view/MITOMAP/MutationsCodingControl
#=> MutationsCodingControl

# Parse MutationsRNA and MutationsCodingControl
echo "Position\tLocus\tDisease\tAllele\tRNA\tHomoplasmy\tHeteroplasmy\tStatus\ttRNA Score\tGB Frequency\tReferences"  >  ${outfile}.tRNA.txt
grep GenomeLoci MutationsRNA | perl -pe 's/\<.*?\>//g; s/\[|\],//g' | perl -pe 's/\",\"/\t/g; s/^\"|\"$//g; s/\"\]//' >> ${outfile}.tRNA.txt

echo "Position\tLocus\tDisease\tAllele\tNucleotide Change\tAmino Acid Change\tHomoplasmy\tHeteroplasmy\tStatus\tGB Frequency\tReferences" >  ${outfile}.Coding.txt
grep GenomeLoci MutationsCodingControl | perl -pe 's/\<.*?\>//g; s/\[|\],//g' | perl -pe 's/\",\"/\t/g; s/^\"|\"$//g; s/\"\]//'           >> ${outfile}.Coding.txt

perl -F"\t" -lane 'next if /Allele/ && $. == 1; print join("\t", "MT", @F[3,4,1,8,2])' MITOMAP.cfrm_report_all.${day}.Coding.txt > ${outfile}.genericdb.Coding.tmp
#perl -i -ne 'print unless /15498del24|15649-15666del/' ${outfile}.genericdb.Coding.tmp
#perl -i -pe 's/\t\w+?(\d+)\w+\t/\t$1\t/; s/\-/\t/; s/\tdel\t/\t-\t/; s/Confirmed/Cfrm/' ${outfile}.genericdb.Coding.tmp
perl -i -pe 's/\t[ACGT](\d+)[ACGT]\t/\t$1\t/; s/\-/\t/' ${outfile}.genericdb.Coding.tmp
perl -i -pe 's/Confirmed/Cfrm/'                         ${outfile}.genericdb.Coding.tmp

perl -F"\t" -lane 'next if /Allele/ && $. == 1; print join("\t", "MT", @F[3,1,7,2])' MITOMAP.cfrm_report_all.${day}.tRNA.txt > ${outfile}.genericdb.tRNA.tmp
perl -i -pe 's/\t([ACGT]+)(\d+)([ACGT]+)\t/\t$2\t$1\t$3\t/' ${outfile}.genericdb.tRNA.tmp
perl -i -pe 's/Confirmed/Cfrm/'                             ${outfile}.genericdb.tRNA.tmp

cat ${outfile}.genericdb.Coding.tmp ${outfile}.genericdb.tRNA.tmp > ${outfile}.genericdb
perl -i -pe 's/A302ACC\tA\tACC/302\t-\tCC/'                                         ${outfile}.genericdb
perl -i -pe 's/C309CC\tC\tCC/309\t-\tC/'                                            ${outfile}.genericdb
perl -i -pe 's/T310TC\tT\tTC/310\t-\tC/'                                            ${outfile}.genericdb
perl -i -pe 's/C960del/960\tC\t-/'                                                  ${outfile}.genericdb
perl -i -pe 's/960\tC\tCC/960\t-\tC/'                                               ${outfile}.genericdb
perl -i -pe 's/961\tT\tTC/961\t-\tC/'                                               ${outfile}.genericdb
perl -i -pe 's/T3271del/3271\tT\t-/'                                                ${outfile}.genericdb
perl -i -pe 's/4322\tC\tCC/4322\t-\tC/'                                             ${outfile}.genericdb
perl -i -pe 's/4369\tA\tAA/4369\t-\tA/'                                             ${outfile}.genericdb
perl -i -pe 's/A5001AA\tA\tAA/5001\t-\tA/'                                          ${outfile}.genericdb
perl -i -pe 's/AA5134d\tAAA\tA/5134\tAA\t-/'                                        ${outfile}.genericdb
perl -i -pe 's/A5537insT/5537\t-\tT/'                                               ${outfile}.genericdb
perl -i -pe 's/CGAGC6020d\tCGAGC\tdel/6020\tCGAGC\t-/'                              ${outfile}.genericdb
perl -i -pe 's/A6698del\tA\tdel/6698\tA\t-/'                                        ${outfile}.genericdb
perl -i -pe 's/C7402del\tC\tdel/7402\tC\t-/'                                        ${outfile}.genericdb
perl -i -pe 's/C7471insC/7471\t-\tC/'                                               ${outfile}.genericdb
perl -i -pe 's/A7472insC/7472\t-\tC/'                                               ${outfile}.genericdb
perl -i -pe 's/8042delAT\tAT\tdel/8042\tAT\t-/'                                     ${outfile}.genericdb
perl -i -pe 's/C8611CC\tC\tCC/8611\t-\tC/'                                          ${outfile}.genericdb
perl -i -pe 's/G8156del\tG\tdel/8156\tG\t-/'                                        ${outfile}.genericdb
perl -i -pe 's/T8618TT\tT\tTT/8618\t-\tT/'                                          ${outfile}.genericdb
perl -i -pe 's/9127\t9128delAT\tAT-del/9127\tAT\t-/'                                ${outfile}.genericdb
perl -i -pe 's/9205\t9206delTA\tTA-del/9205\tTA\t-/'                                ${outfile}.genericdb
perl -i -pe 's/9480del15\tTTTTTCTTCGCAGGA\tdel/9480\tTTTTTCTTCGCAGGA\t-/'           ${outfile}.genericdb
perl -i -pe 's/C9537insC\tC\tCC/9537\t-\tC/'                                        ${outfile}.genericdb
perl -i -pe 's/C9559del\tC\tdel/9559\tC\t-/'                                        ${outfile}.genericdb
perl -i -pe 's/11621delTA\tTA\tdel/11621\tTA\t-/'                                   ${outfile}.genericdb
perl -i -pe 's/A12425del\tA\tdel/12425\tA\t-/'                                      ${outfile}.genericdb
perl -i -pe 's/14787delTTAA\tTTAA\tdel/14787\tTTAA\t-/'                             ${outfile}.genericdb
perl -i -pe 's/15498del24\t24bpdeletion\tMT\tCYB/15498del24\t24bpdeletion\tMT-CYB/' ${outfile}.genericdb
perl -i -pe 's/T15944del/15944\tT\t-/'                                              ${outfile}.genericdb
perl -i -pe 's/16018\tT\tTTCTCTGTTCTTTCAT/16018\t-\tTCTCTGTTCTTTCAT/'               ${outfile}.genericdb
perl -i -pe 's/16021_16022delCT/11621\tTA\t-/'                                      ${outfile}.genericdb
perl -i -pe 's/16032\tT\tTTCTCTGTTCTTTCAT/16032\t-\tTCTCTGTTCTTTCAT/'               ${outfile}.genericdb
perl -i -pe 's/16033\tG\tTCTCTGTTCTTTCATG/16033\t-\tTCTCTGTTCTTTCAT/'               ${outfile}.genericdb

perl -i".with_removed_entries" -ne 'print unless /MT\t3902_3908invACCTTGC|MT\t15498del24|MT\tT961delT\+|MT\t15649\t15666del/' ${outfile}.genericdb

ruby -i -F"\t" -lane 'puts [$F[0], $F[1], $F[1].to_i + $F[2].size - 1, $F[2], $F[3], $F[5] + ": " + $F[6]].join("\t")'         ${outfile}.genericdb

cd ${current_dir}


echo "Currently, these are removed, due to difficult to format"
echo "MT\t3902_3908invACCTTGC"
echo "MT\t15498del24"
echo "MT\tT961delT\+"
echo "MT\t15649\t15666del"
echo "They should be formatted for Annovar or grep"

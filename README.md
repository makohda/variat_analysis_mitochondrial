# variat_analysis_mitochondrial

##  パスなど
- cd /lustre/data/mako/tmp/
- module load matchclips2/20140409
- module load samtools/1.6
- ln -s /usr/local/matchclips2-20140409/matchclips ./ # module load だとパスが通らない感じだったので


## mtDNA deletion の検出
数kbp欠失していることがあるが、GATK variant callでは検出されないので別の手段が必要になる  
この時  
1. 何らかの手法で欠失を検出
2. IGV で目視で確認 (この時ペアエンド表示にしていると一目でわかる)

あたりが簡単な方法と思われる  
サンプル数が増えると2から1になっていく  
ここではペアエンド間の距離、スプリッドリード、シークエンスデプスを用いた解析を行える matchclips2 での解析 (1) を行う  

### Pre-require
- rCRS.fasta の用意  
```
wget https://www.mitomap.org/foswiki/pub/MITOWIKI/ArchivedProjectDataCuration/rCRS.fasta
cp rCRS.fasta rCRS_MT.fasta
emacs rCRS_MT.fasta # fastaのヘッダー行を MT に書き換える
samtools faidx rCRS_MT.fasta
```

### MatchClips2 を用いた mtDNA deletion の検出
yhwu/matchclips2: paired end distance, split reads matching, read depth are used to detect and precisely locate structure variations https://github.com/yhwu/matchclips2  

テストしてみる  
```
./matchclips -t 4 -L 10000 -f rCRS_MT.fasta -b Pt1812.MT.bam -o Pt1812.matchclips.out
./matchclips -t 4 -L 10000 -f rCRS_MT.fasta -b Pt1088*.MT.bam -o Pt1088.matchclips.out
./matchclips -t 4 -L 10000 -f rCRS_MT.fasta -b Pt0798*.MT.bam -o Pt0798.matchclips.out
./matchclips -t 4 -L 10000 -f rCRS_MT.fasta -b Pt1422*.MT.bam -o Pt1422.matchclips.out
```

全サンプルで  
```
foreach i in ./Pt*MT.bam
  echo $i
  ./matchclips -t 4 -L 10000 -f rCRS.fasta -b ${i} -o ${i:r}.matchclips.out
end
```
(検出としてはいいが、微妙にポジション違う)

### 課題
- 目視と比較して数塩基ズレてるように見える場合がある (仕方ないが)
- 感度がどの程度か (これは欠失あり検体データに、健常データを混ぜて、欠失の割合をコントロールしたサンプルデータをつくれば簡単に検証できそう)


### 正解データ
Pt0372  
Pt0798  
Pt1088  
Pt1422  
Pt1812  
(内部データベースから引いてください)


----

## mtDNA mutation の検出 (Annovar + MITOMAP)
ミトコンドリアゲノム上にある変異は dbSNP, ClinVar, OMIM などなどいろいろ記載があるが MITOMAP のものが最もキュレーションされているようなので、これをアノテーションとして付加していく  

### 既知の変異例 (この後の実行例の正解例)
内部DBから引いてください

### Annovar でのミトコンドリアバリアントのアノテーション
combined_genotyped.vcf は GATK でバリアントコールして生成した .vcf から MT だけ抜き出したもの  
この時、１つ注意として他の染色体のバリアントのようにバリアントフィルタリングをする前の combined_genotyped.vcf から作ること  
GATK の様々なステップは常染色体を前提としているので、既知の変異もバリアントフィルタリングで除かれる場合があるため  

combined_genotyped.MT.vcf の作成
```
grep -e "^#" -e "^MT" combined_genotyped.vcf > combined_genotyped.MT.vcf
./annovar/convert2annovar.pl -format vcf4 --includeinfo --withzyg \
                             --allsample combined_genotyped.MT.vcf \
                             --outfile combined_genotyped.MT
```

Annovar解析の実行
```
id=Pt0871BL_F_C1P_MC__.SSv5
./annovar/table_annovar.pl combined_genotyped.MT.${id}.avinput annovar/humandb/ \
                           -buildver GRCh37_MT \
                           --protocol ensGene,generic,generic,generic,generic,generic \
                           --genericdbfile MITOMAP.cfrm_report_all.180809.genericdb,tommo-3.5kjpnv2-20180625-af_snvall.MAF.genericdb,tommo-3.5kjpnv2-20180625-af_snvall.INFO.genericdb,clinvar_20170905.OMIM.genericdb,clinvar_20170905.genericdb \
                           --operation g,f,f,f,f,f \
                           --argument '--hgvs --exonicsplicing --splicing_threshold 5',,,,, \
                           --nastring NA \
                           --otherinfo \
                           --remove \
                           --outfile ${id}.MT.avoutput #=> PtXXXX.MT.avoutput.GRCh37_MT_multianno.txt

perl -i -pe 's/generic/MITOMAP.180809/'           ${id}.MT.avoutput.GRCh37_MT_multianno.txt
perl -i -pe 's/generic2/3.5KJPNv2.20180625.MAF/'  ${id}.MT.avoutput.GRCh37_MT_multianno.txt
perl -i -pe 's/generic3/3.5KJPNv2.20180625.INFO/' ${id}.MT.avoutput.GRCh37_MT_multianno.txt
perl -i -pe 's/generic4/clinvar_20170905.OMIM/'   ${id}.MT.avoutput.GRCh37_MT_multianno.txt
perl -i -pe 's/generic5/clinvar_20170905/'        ${id}.MT.avoutput.GRCh37_MT_multianno.txt

cat Pt0871*.MT.avoutput.GRCh37_MT_multianno.txt B grep Cfrm CL
```

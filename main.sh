#!/usr/bin/env bash

TMP_DIR=tmp

if [ ! -e "$TMP_DIR" ]; then
	mkdir $TMP_DIR
fi

if [ ! -e "$TMP_DIR/img" ]; then
	mkdir $TMP_DIR/img
fi

if [ ! -e "$TMP_DIR/img/site" ]; then
	mkdir $TMP_DIR/img/site
fi

if [ ! -e "$TMP_DIR/img/gene" ]; then
	mkdir $TMP_DIR/img/gene
fi


for fam in $(seq 1 1 6)
#for fam in $(seq 1 1 1)
	do
	for t in SI ST
#	for t in SI 
		do 
		echo "processing..: $fam, $t"

		./get_ASE_applicables.pl all_family.SNP.genes.vcf_gff id.tsv $fam $t > $TMP_DIR/fam_$fam.$t.tsv
		./filter_ASE_sites.pl $TMP_DIR/fam_$fam.$t.tsv $fam $t $TMP_DIR/fam_$fam.$t.filtered.tsv $TMP_DIR/fam_$fam.$t.filtered.gene.tsv 

# Per SNP
		cat $TMP_DIR/fam_$fam.$t.filtered.tsv | perl -lanF/"\t"/ -e 'if (m/^#/) {tr/#//d; print} elsif ($F[4]!=0 && $F[5]!=0 && $F[6]!=0 && $F[7]!=0) {print};' > $TMP_DIR/fam_$fam.$t.filtered.non-zero.tsv
		cat $TMP_DIR/fam_$fam.$t.filtered.tsv | perl -lanF/"\t"/ -e 'if (m/^#/) {print} elsif ($F[4]==0 || $F[5]==0 || $F[6]==0 || $F[7]==0) {print};' > $TMP_DIR/fam_$fam.$t.filtered.zero.tsv

# Per Gene, simple average
		cat $TMP_DIR/fam_$fam.$t.filtered.gene.tsv | perl -lanF/"\t"/ -e 'if (m/^#/) {tr/#//d; print} elsif ($F[6]>=0.5 && $F[7]>=0.5 && $F[8]>=0.5 && $F[9]>=0.5) {print};' > $TMP_DIR/fam_$fam.$t.filtered.gene.non-zero.tsv
		cat $TMP_DIR/fam_$fam.$t.filtered.gene.tsv | perl -lanF/"\t"/ -e 'if (m/^#/) {print} elsif ($F[6]<0.5 || $F[7]<0.5 || $F[8]<0.5 || $F[9]<0.5) {print};' > $TMP_DIR/fam_$fam.$t.filtered.gene.zero.tsv


		if [[ $fam -eq 4 || $fam -eq 6 ]]; then
			PL_P1=2
			PL_P2=4
		else
			PL_P1=2
			PL_P2=2
		fi
		
		echo "...statistical testing [site]"
		R CMD BATCH "--args $TMP_DIR/fam_$fam.$t.filtered.non-zero.tsv $TMP_DIR/fam_$fam.$t.filtered.non-zero.stat.tsv $PL_P1 $PL_P2" stat_regulatory.r $TMP_DIR/stat_regulatory.r.fam_$fam.$t.Rout
		echo "...plotting [site]"
		R CMD BATCH "--args $TMP_DIR/fam_$fam.$t.filtered.non-zero.stat.tsv $TMP_DIR/img/site/fam_$fam.$t.relative $TMP_DIR/fam_$fam.$t.filtered.non-zero.stat.regulatory.tsv $PL_P1 $PL_P2" scatter_plot.r $TMP_DIR/scatter_plot.r.fam_$fam.$t.Rout

		echo "...statistical testing [gene]"
		R CMD BATCH "--args $TMP_DIR/fam_$fam.$t.filtered.gene.non-zero.tsv $TMP_DIR/fam_$fam.$t.filtered.gene.non-zero.stat.tsv $PL_P1 $PL_P2" stat_regulatory.r $TMP_DIR/stat_regulatory.r.fam_$fam.$t.gene.Rout
		echo "...plotting [gene]"
		R CMD BATCH "--args $TMP_DIR/fam_$fam.$t.filtered.gene.non-zero.stat.tsv $TMP_DIR/img/gene/fam_$fam.$t.relative $TMP_DIR/fam_$fam.$t.filtered.gene.non-zero.stat.regulatory.tsv $PL_P1 $PL_P2" scatter_plot.r $TMP_DIR/scatter_plot.r.fam_$fam.$t.gene.Rout
	done
done




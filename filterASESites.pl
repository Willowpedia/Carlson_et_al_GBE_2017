#!/usr/bin/env perl

use strict;
use List::Util qw (sum);


# FILE_EXPRESSION
#GENE_ID	CHR	POS	REF	ALT	A	C	G	T	DNA_P1_A	DNA_P1_C	DNA_P1_G	DNA_P1_T	DNA_READ_DEPTH_P1	DNA_P2_A	DNA_P2_C	DNA_P2_G	DNA_P2_TDNA_READ_DEPTH_P2	RNA_P1_A	RNA_P1_C	RNA_P1_G	RNA_P1_T	RNA_P1_TOTAL	RNA_P2_A	RNA_P2_C	RNA_P2_G	RNA_P2_T	RNA_P2_TOTAL	RNA_H_A	RNA_H_C	RNA_H_G	RNA_H_T	RNA_H_TOTAL
#SapurV1A.0646s0100	chr16	40182	C	T	0	P1	0	P2	0	77	0	0	77	0	0	0	26	26	0	0	0	0	0	0	0	0
#SapurV1A.0646s0090	chr16	63437	A	G	P1	0	P2	0	68	0	0	0	68	0	0	8	0	8	7	0	0	0	7	0	0	6

# FILE_READ_COUNT
#P1	26424111
#P2	29581541
#H	22489012

my $file_expression = $ARGV[0];
my $family_num = $ARGV[1];
my $tissue = $ARGV[2];

my $f_out_filtered = $ARGV[3];
my $f_out_filtered_gene = $ARGV[4];


my $file_read_count = "num_mapped_reads/fam_$family_num.$tissue.mapped_reads_count.tsv";

# config
# per site
my $min_dna_depth = 5; #3;
my $min_dna_match_pct = 95;
my $min_cpm_total = 1; #0.3; #1; # P1+P2+H, matching alleles only
my $min_rna_match_pct = 95;

# per gene
my $min_sites_per_gene = 1;
####

my $MIL = 1000000;

my @libs;
my %mapped_read_counts = ();
open(my $fh,"<$file_read_count");
while(<$fh>) {
	chomp;
	my ($lib, $total_count) = split/\t/;
	$mapped_read_counts{$lib} = $total_count;

	push(@libs, $lib);
}
close($fh);

#header
#print join("\t", qw(#CHR POS ALLELE_P1 ALLELE_P2 GENE_ID DNA_MATCHED_PCT_P1 DNA_MATCHED_P1 DNA_UNMATCHED_P1 DNA_MATCHED_PCT_P2 DNA_MATCHED_P2 DNA_UNMATCHED_P2));
#foreach my $lib (@libs) {
#	print "\tCPM_$lib\tCOUNT_$lib\tCPM_P1_$lib\tCPM_P2_$lib";
#}

#print "\n";

# per gene information
my %num_sites = {};
my %sum_rna_reads = {};

open(my $fh_out_site, ">$f_out_filtered") or die();


open(my $fh,"<$file_expression");

while (<$fh>) {
	chomp;
	if (m/^#/) {
		tr/#//d;
		printf $fh_out_site join("\t", qw( #NUM_DNA_READS_P1 NUM_DNA_READS_P2 NUM_RNA_READS_P1_IN_P1 NUM_RNA_READS_P2_IN_P2 NUM_RNA_READS_P1_IN_P1_NORM NUM_RNA_READS_P2_IN_P2_NORM NUM_RNA_READS_P1_IN_H NUM_RNA_READS_P2_IN_H CPM_P1_IN_P1 CPM_P2_IN_P2 CPM_P1_IN_H CPM_P2_IN_H ), $_) . "\n";
		next;
	}
	my @F = split/\t/;
	my @alleles = @F[5..8];
	my @DNA_P1 = @F[9..13];
	my @DNA_P2 = @F[14..18];
	my @RNA_P1 = @F[19..23];
	my @RNA_P2 = @F[24..28];
	my @RNA_H = @F[29..33];

	my $dna_reads_p1 = 0;
	my $dna_reads_p2 = 0;
	my $dna_reads_p1_unmatched = 0;
	my $dna_reads_p2_unmatched = 0;
	my $exp_p1 = 0;
	my $exp_p2 = 0;
	my $exp_p1h = 0;
	my $exp_p2h = 0;

	foreach my $i (0..3) {
		if ($alleles[$i] eq "P1") {
			$dna_reads_p1 += $DNA_P1[$i];
			$exp_p1 += $RNA_P1[$i];
			$exp_p1h += $RNA_H[$i];
			$dna_reads_p2_unmatched += $DNA_P2[$i];
		} elsif ($alleles[$i] eq "P2") {
			$dna_reads_p2 += $DNA_P2[$i];
			$exp_p2 += $RNA_P2[$i];
			$exp_p2h += $RNA_H[$i];
			$dna_reads_p1_unmatched += $DNA_P1[$i];
		} else {
			$dna_reads_p1_unmatched += $DNA_P1[$i];
			$dna_reads_p2_unmatched += $DNA_P2[$i];
		}
	}

	# minimum DNA reads count
	if ($dna_reads_p1 < $min_dna_depth || $dna_reads_p2 < $min_dna_depth) {
		next;
	}

	# proportion of matched DNA reads
	if ($dna_reads_p1*100/($dna_reads_p1+$dna_reads_p1_unmatched) < $min_dna_match_pct) {
		next;
	}
	if ($dna_reads_p2*100/($dna_reads_p2+$dna_reads_p2_unmatched) < $min_dna_match_pct) {
		next;
	}

	my $exp_p1_cpm = sprintf("%.4f", $exp_p1*$MIL/$mapped_read_counts{"P1"});
	my $exp_p2_cpm = sprintf("%.4f", $exp_p2*$MIL/$mapped_read_counts{"P2"});
	my $exp_p1h_cpm = sprintf("%.4f", $exp_p1h*$MIL/$mapped_read_counts{"H"});
	my $exp_p2h_cpm = sprintf("%.4f", $exp_p2h*$MIL/$mapped_read_counts{"H"});

	# expression level
	if (sum($exp_p1_cpm, $exp_p2_cpm, $exp_p1h_cpm, $exp_p2h_cpm) < $min_cpm_total) {
	#if (sum($exp_p1_cpm, $exp_p2_cpm) < $min_cpm_total) {
		next;
	}

	# proportion of matched RNA reads
	if ($RNA_P1[4] != 0 && $exp_p1*100/$RNA_P1[4] < $min_rna_match_pct) {
		next;
	}
	if ($RNA_P2[4] != 0 && $exp_p2*100/$RNA_P2[4] < $min_rna_match_pct) {
		next;
	}
	if ($RNA_H[4] != 0 && ($exp_p1h+$exp_p2h)*100/$RNA_H[4] < $min_rna_match_pct) {
		next;
	}

	my $exp_p1_norm;
	my $exp_p2_norm;
	if ($mapped_read_counts{"P1"} > $mapped_read_counts{"P2"}) {
		$exp_p1_norm = sprintf("%.0f", $exp_p1*$mapped_read_counts{"P2"}/$mapped_read_counts{"P1"});
		$exp_p2_norm = $exp_p2;
	} else {
		$exp_p1_norm = $exp_p1;
		$exp_p2_norm = sprintf("%.0f", $exp_p2*$mapped_read_counts{"P1"}/$mapped_read_counts{"P2"});
	}
	
	printf $fh_out_site join("\t", $dna_reads_p1, $dna_reads_p2, $exp_p1, $exp_p2, $exp_p1_norm, $exp_p2_norm, $exp_p1h, $exp_p2h, $exp_p1_cpm, $exp_p2_cpm, $exp_p1h_cpm, $exp_p2h_cpm, @F) ."\n";


	# store information per gene
	my $gene_id = $F[0];
	if (exists($num_sites{$gene_id})) {
		$num_sites{$gene_id}++;
		$sum_rna_reads{$gene_id}->[0] += $exp_p1;
		$sum_rna_reads{$gene_id}->[1] += $exp_p2;
		$sum_rna_reads{$gene_id}->[2] += $exp_p1h;
		$sum_rna_reads{$gene_id}->[3] += $exp_p2h;
	} else {
		$num_sites{$gene_id} = 1;
		$sum_rna_reads{$gene_id} = [ ($exp_p1, $exp_p2, $exp_p1h, $exp_p2h) ];
	}
}

close($fh);
close($fh_out_site);

open(my $fh_out_gene, ">$f_out_filtered_gene") or die();

printf $fh_out_gene join("\t", qw( #GENE_ID SUM_RNA_READS_P1_IN_P1 SUM_RNA_READS_P2_IN_P2 SUM_RNA_READS_P1_IN_H SUM_RNA_READS_P2_IN_H NUM_SITES NUM_RNA_READS_P1_IN_P1_NORM NUM_RNA_READS_P2_IN_P2_NORM NUM_RNA_READS_P1_IN_H NUM_RNA_READS_P2_IN_H CPM_P1_IN_P1 CPM_P2_IN_P2 CPM_P1_IN_H CPM_P2_IN_H )) . "\n";

foreach my $k (sort keys %num_sites) {
	if ($num_sites{$k} < $min_sites_per_gene) {
		next;
	}
	my $cpm_p1_avg = $sum_rna_reads{$k}->[0]*$MIL/$mapped_read_counts{"P1"}/$num_sites{$k};
	my $cpm_p2_avg = $sum_rna_reads{$k}->[1]*$MIL/$mapped_read_counts{"P2"}/$num_sites{$k};
	my $cpm_p1h_avg = $sum_rna_reads{$k}->[2]*$MIL/$mapped_read_counts{"H"}/$num_sites{$k};
	my $cpm_p2h_avg = $sum_rna_reads{$k}->[3]*$MIL/$mapped_read_counts{"H"}/$num_sites{$k};

	my @rna_reads_avg = ();
	for my $i (0..3) {
		$rna_reads_avg[$i] = $sum_rna_reads{$k}->[$i]/$num_sites{$k};
	}
	if ($mapped_read_counts{"P1"} > $mapped_read_counts{"P2"}) {
		$rna_reads_avg[0] *= $mapped_read_counts{"P2"}/$mapped_read_counts{"P1"};
	} else {
		$rna_reads_avg[1] *= $mapped_read_counts{"P1"}/$mapped_read_counts{"P2"};
	}
	
	printf $fh_out_gene join("\t", $k, @{$sum_rna_reads{$k}}, $num_sites{$k}, @rna_reads_avg, $cpm_p1_avg, $cpm_p2_avg, $cpm_p1h_avg, $cpm_p2h_avg) . "\n";
}


close($fh_out_gene);



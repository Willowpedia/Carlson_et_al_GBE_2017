#!/usr/bin/env perl

use strict;
use List::Util qw (sum);
use List::MoreUtils qw(uniq);

# INPUT FILE
#1	#CHROM
#2	POS
#3	ID
#4	REF
#5	ALT
#6	QUAL
#7	FILTER
#8	INFO
#9	FORMAT
#10	01-200-003
#11	07-MBG-5027
#12	94001
#13	94006
#14	9882-34
#15	9882-41
#16	Jorr
#17	seqname, CHROM (GFF part starts)
#18	source, JGI_JCVI
#19	feature, gene
#20	start
#21	end
#22	score, .
#23	strand
#24	frame
#25	attribute (e.g. ID=SapurV1A.0035s0010;Name=SapurV1A.0035s0010)

my %column_num = (
	"01-200-003" => 10,
	"07-MBG-5027" => 11,
	"94001" => 12,
	"94006" => 13,
	"9882-34" => 14,
	"9882-41" => 15,
	"Jorr" => 16,
);

my $file = $ARGV[0];
my $id_file = $ARGV[1];
my $family_num = $ARGV[2];
my $tissue = $ARGV[3];

my $idx_p1;
my $idx_p2;
my $ID_P1;
my $ID_P2;
my $ID_H;

open(my $fh,"<$id_file") or die "cannot open: $!";
while(<$fh>) {
	chomp;
	my @F = split/\t/;
	if ($F[0] == $family_num) {
		$idx_p1 = $column_num{$F[1]} -1;
		$idx_p2 = $column_num{$F[2]} -1;
		$ID_P1 = $F[3];
		$ID_P2 = $F[4];
		$ID_H = $F[5];
		last;
	};
}
close($fh);

my $RNA_map_file_P1 = "mapping_coverage/$ID_P1-$tissue.SNP_only.tsv";
my $RNA_map_file_P2 = "mapping_coverage/$ID_P2-$tissue.SNP_only.tsv";
my $RNA_map_file_H = "mapping_coverage/$ID_H-$tissue.SNP_only.tsv";
my $total_read_count_file = "num_mapped_reads/fam_$family_num.$tissue.mapped_reads_count.tsv";


#print "RNA mapping file for P1: $RNA_map_file_P1\n";


my @libs;
my %mapped_read_counts = ();
open(my $fh,"<$total_read_count_file") or die "cannot open: $!";
while(<$fh>) {
	chomp;
	my ($lib, $total_count) = split/\t/;
	$mapped_read_counts{$lib} = $total_count;

	push(@libs, $lib);
}
close($fh);

my %RNA_data_P1 = read_RNA_mapping_file($RNA_map_file_P1);
my %RNA_data_P2 = read_RNA_mapping_file($RNA_map_file_P2);
my %RNA_data_H = read_RNA_mapping_file($RNA_map_file_H);

#header
print join("\t", qw(#GENE_ID CHR POS REF ALT A C G T
			 DNA_P1_A DNA_P1_C DNA_P1_G DNA_P1_T DNA_READ_DEPTH_P1
			 DNA_P2_A DNA_P2_C DNA_P2_G DNA_P2_T DNA_READ_DEPTH_P2
			 RNA_P1_A RNA_P1_C RNA_P1_G RNA_P1_T RNA_P1_TOTAL
			 RNA_P2_A RNA_P2_C RNA_P2_G RNA_P2_T RNA_P2_TOTAL
			 RNA_H_A RNA_H_C RNA_H_G RNA_H_T RNA_H_TOTAL)) ."\n";
#foreach my $lib (@libs) {
#	print "\tCPM_$lib\tCOUNT_$lib\tCPM_P1_$lib\tCPM_P2_$lib";
#}

#print "\n";


open(my $fh,"<$file") or die "cannot open: $!";

while(<$fh>) {
	chomp;
	my @F = split/\t/;

	my @P1 = split/:/, $F[$idx_p1];
	my @P2 = split/:/, $F[$idx_p2];

	my @allele_p1 = uniq split/\//,$P1[0];
	my @allele_p2 = uniq split/\//,$P2[0];

	# check overlap between P1 and P2 alleles
	my $ASE = 1;
	foreach my $a1 (@allele_p1) {
		if ($a1 eq ".") {
			$ASE = 0;
			last;
		}
		foreach my $a2 (@allele_p2) {
			if ($a2 eq "." || $a1 eq $a2) {
				$ASE = 0;
				last;
			}
		}
	}
	next if ($ASE == 0);

	my @gt;
	foreach my $a1 (@allele_p1) {
		$gt[$a1] = "P1";
	}
	foreach my $a2 (@allele_p2) {
		$gt[$a2] = "P2";
	}

	$F[24] =~ m/ID=([^;]+)/;
	my $ID = $1;

	#print join("\t", @F[0..8], $F[$idx_p1], $F[$idx_p2], @F[19..20], $F[22], $F[24]) ."\n";
	print join("\t", $ID, @F[0..1], $F[3], $F[4],
				get_acgt_order($F[3], $F[4], join(",",@gt)),
				get_acgt_order($F[3], $F[4], $P1[1]), $P1[2],
				get_acgt_order($F[3], $F[4], $P2[1]), $P2[2],
				$RNA_data_P1{join("\t",@F[0..1])},
				$RNA_data_P2{join("\t",@F[0..1])},
				$RNA_data_H{join("\t",@F[0..1])},
				) ."\n";
}
close($fh);

sub get_acgt_order {
	my ($ref_allele, $alt_alleles, $AD_values) =  @_;

	my @NT_order = ("A", "C", "G", "T");

	my @NT = ($ref_allele);
	push(@NT, split(',',$alt_alleles) );
	my %DNA_count;
	foreach my $a (@NT_order) {
		$DNA_count{$a} = 0;
	}
	my @count = split/,/,$AD_values;
	for (my $i=0; $i<=$#count; $i++) {
		$DNA_count{$NT[$i]} = $count[$i];
	}

	my @ret;
	foreach my $a (@NT_order) {
		push(@ret, $DNA_count{$a});
	}

	return @ret;
}


sub read_RNA_mapping_file {
	my ($file) = @_;
	# RNA_map_file
	#chr16	6295	-	A	3	0	0	0	0	0	3
	#chr16	6311	-	G	0	0	4	0	0	0	4
	#chr16	6312	-	T	0	0	0	4	0	0	4
	my %RNA_data = ();
	open(my $fh,"<$file") or die "cannot open: $!";
	while(<$fh>) {
		chomp;
		my @F = split/\t/;
		$RNA_data{join("\t",@F[0..1])} = join("\t", @F[4..7], $F[10]);
	}
	close($fh);

	return %RNA_data;
}



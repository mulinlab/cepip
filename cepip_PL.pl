#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

=head1 Name

cepip_PL.pl  --   Perl version of cepip for prediction cell type-specific regulatory variant

=head1 Version

  Author: Mulin Jun Li, mulin0424.li@gmail.com
  Version: 0.12,  Date: 2016-12-12
  Note:

=head1 Usage

  -h			output help information to screen
  -i			the input variants file
  -f			the format of input file, supporting VCF and ANNOVAR format; default: vcf
  -t			the path of executable tabix program; default: /user/bin/tabix
  -r			the path of dbNCFP reference database; default: ftp://147.8.193.36/PRVCS/v1.1/dbNCFP_whole_genome_SNVs.bgz
  -s			the selected tool scores; default: 1,3,4,5,6,8,10,11
  				1	CADD_cscore
  				2	CADD_PHRED
  				3	DANN_score
  				4	FunSeq_score
  				5	FunSeq2_score
  				6	GWAS3D_score
  				7	GWAVA_region_score
  				8	GWAVA_TSS_score
  				9	GWAVA_unmatched_score
  				10	SuRFR_score
  				11	Fathmm_MKL_score
  -p			the casual distribution folder; default: resources/all_causal_distribution
  -n			the neutral distribution folder; default: resources/all_neutral_distribution
  -a			the path of tissue/cell type reference epigenome; default: ftp://147.8.193.36/cepip/cell_signal/
  -c			the dependent tissue/cell type; default: E116
  -o			the path of output file; default: cepip1.flt.txt
  
=head1 Exmple

perl cepip_PL.pl -i query.vcf -f vcf -t /user/bin/tabix -r ftp://147.8.193.36/PRVCS/v1.1/dbNCFP_whole_genome_SNVs.bgz -s 1,3,4,5,6,8,10,11 -a ftp://147.8.193.36/cepip/cell_signal/ -c E116 -o cepip1.flt.txt
	
=cut


my $TABIX = "/g/jun/Tools/tabix-0.2.6/tabix";
my $input_file = '';
my $input_format = "vcf";
my $output_file = "cepip1.flt.txt";
my $reference_file = "dbNCFP_whole_genome_SNVs.bgz";
my $selected_tools = "1,3,4,5,6,8,10,11";
my $positive_dis = "resources/all_causal_distribution";
my $negative_dis = "resources/all_neutral_distribution";
my $cell_type_annotation = "resources/cell_type_annotation/";
my $CELL_ID = "E116";
my $cell_type_annotation_list = "resources/cell_signal";

open(SL,$cell_type_annotation_list);
my %cellsg2hash =();
while(<SL>){
	chomp;
	$_=~ s/[\r\n]+$//;
	$cellsg2hash{$_} = 1;
}

my $DNASE = "DNase.macs2.narrowPeak.sorted.gz";
my $DNASE_NE = "DNase.narrowPeak.sorted.gz";
my $H2AZ = "H2A.Z.narrowPeak.sorted.gz";
my $H3K27ac = "H3K27ac.narrowPeak.sorted.gz";
my $H3K27me3 = "H3K27me3.narrowPeak.sorted.gz";
my $H3K36me3 = "H3K36me3.narrowPeak.sorted.gz";
my $H3K4me1 = "H3K4me1.narrowPeak.sorted.gz";
my $H3K4me2 = "H3K4me2.narrowPeak.sorted.gz";
my $H3K4me3 = "H3K4me3.narrowPeak.sorted.gz";
my $H3K79me2 = "H3K79me2.narrowPeak.sorted.gz";
my $H3K9ac = "H3K9ac.narrowPeak.sorted.gz";
my $H3K9me3 = "H3K9me3.narrowPeak.sorted.gz";
my $H4K20me1 = "H4K20me1.narrowPeak.sorted.gz";
my $help;

my %dis_index = ("1"=>"CADD_cscore","2"=>"CADD_PHRED","3"=>"DANN_score","4"=>"FunSeq_score","5"=>"FunSeq2_score","6"=>"GWAS3D_score","7"=>"GWAVA_region_score","8"=>"GWAVA_TSS_score","9"=>"GWAVA_unmatched_score","10"=>"SuRFR_score","11"=>"Fathmm_MKL_score");
my %mean_index = ("1"=>-0.0784622,"2"=>3.68016,"3"=>0.291601,"4"=>1,"5"=>1E-50,"6"=>2.80186,"7"=>0.153309,"8"=>0.205698,"9"=>0.33259,"10"=>9.44965,"11"=>0.00403947);

my @ARGVR = @ARGV;
while (defined(my $arg=shift(@ARGVR))){
	if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { die `pod2text $0`; }
}

my $result = GetOptions(
	'help|?' => \$help,
    "i=s" => \$input_file,
    "f=s" => \$input_format,
    "t=s" => \$TABIX,
    "s=s" => \$selected_tools,
    "p=s" => \$positive_dis,
    "n=s" => \$negative_dis,
    "a=s" => \$cell_type_annotation,
    "c=s" => \$CELL_ID,
    "r=s" => \$reference_file,
    "o=s" => \$output_file
);

if(!-f $input_file){
	die "the input file is not found!\n";
}

#if(!-f $reference_file){
#	die "the reference database file is not found!\n";
#}else{
#	my $reference_file_index = $reference_file . ".tbi";
#	if(!-f $reference_file_index){
#		die "the index file for reference database is required!\n";
#	}
#}

if(!-X $TABIX){
	die "the TABIX is not righlty set!\n";
}

#init selected tools
my @selected_tools_list = split(/,/,$selected_tools);
my %selected_dis_hash = ();
my $print_header = "CHR\tPOS\tREF\tALT\t";
foreach my $sd (@selected_tools_list){
	my $dis_p_file = $positive_dis."/".$dis_index{$sd}.".dis";
	my @dis_p_array = ();
	open(PDIS,$dis_p_file);
	while(<PDIS>){
		chomp;
		push(@dis_p_array,$_);		
	}
	close(PDIS);
	my $p_variance = pop(@dis_p_array) + 0;
	my $p_mean = pop(@dis_p_array) + 0;
#	print $dis_index{$sd}."\t".$p_mean."\t".$p_variance."\n";
	
	my $dis_n_file = $negative_dis."/".$dis_index{$sd}.".dis";
	my @dis_n_array = ();
	open(NDIS,$dis_n_file);
	while(<NDIS>){
		chomp;
		push(@dis_n_array,$_);		
	}
	close(NDIS);
	my $n_variance = pop(@dis_n_array) + 0;
	my $n_mean = pop(@dis_n_array) + 0;
	
	my %c_dis_obj = ("p_dis"=>\@dis_p_array,"p_mean"=>$p_mean,"p_variance"=>$p_variance, "n_dis"=>\@dis_n_array,"n_mean"=>$n_mean,"n_variance"=>$n_variance); 
	$selected_dis_hash{$sd} = \%c_dis_obj;
#	print $dis_index{$sd}."\t".$n_mean."\t".$n_variance."\n\n";
	$print_header .= $dis_index{$sd}."\t";
}
$print_header .= "BF\tComposite_score\tCell_P\tCombine_P\n";
open(OUT,">".$output_file);
print OUT $print_header;

#print OUT "Chr\tPos\tRef\tAlt\tCADD_Cscore\tCADD_PHRED\tDANN_score\tFunSeq_score\tFunSeq2_score\tGWAS3D_score\tGWAVA_region_score\tGWAVA_TSS_score\tGWAVA_unmatched_score\tSuRFR_score\tFathmm_MKL_score\n";
open(IN,$input_file);
while(<IN>){
	chomp;
	$_ =~ s/[\s\r\n]+$//g;
	next if($_ eq '');
	next if($_ =~ /^#/);
	my ($chr, $position, $id_or_position2, $ref, $alt);
	if($input_format =~ /annovar/i ){
		($chr, $position, $id_or_position2, $ref, $alt) = split(/\t/,$_);
	}else{
		($chr, $position, $id_or_position2, $ref, $alt) = split(/\t/,$_);
	}
	
	if($ref !~ /[ATCG]+/ || $alt !~ /[ATCG]+/){
		die "make sure your input is VCF or ANNOVAR format, and each variant should contain ref and at least one alt allele!\n";
	}
	
	my $isIndel = 0;
	if(length($ref) > 1){
		$isIndel = 1;
	}
	my @alts = split (/,/, $alt);
	foreach my $alt_one (@alts) {
		if (length($alt_one) > 1) {
			$isIndel = 1;
		}
	}
	next if($isIndel ==1);
	
	my $query = $chr.':'.$position.'-'.$position;
#	print "$TABIX $reference_file $query\n";
	my $T_ALL = `$TABIX $reference_file $query`;
	my $curr_snp_info = '';
	my @T_ALL_array = split(/\n/,$T_ALL);
	my %allelehash_funseq2_cs = (); #need fix funseq2 score multiple line
	foreach my $one_site (@T_ALL_array){
		my ($chr_o, $position_o) = split(/\t/,$one_site);
		if($position_o == $position){
			$curr_snp_info = $one_site;
		}
	}
	
	if($curr_snp_info ne ''){
		#print $curr_snp_info."\n";
		my @site_array = split(/\t/,$curr_snp_info);
		
		#CADD 
		my %allelehash_cadd_cs = ();
		my %allelehash_cadd_ph = ();
		my @cadd_alleles = split(/,/,(split(/\|/,$site_array[2]))[1]);
		if($site_array[2] ne '.'){
			my @cadd_score_cs = split(/,/,$site_array[3]);
			my @cadd_score_ph = split(/,/,$site_array[4]);
			for (my $i=0;$i<scalar(@cadd_alleles);$i++){
				$allelehash_cadd_cs{$ref.":".$cadd_alleles[$i]} = $cadd_score_cs[$i];
				$allelehash_cadd_ph{$ref.":".$cadd_alleles[$i]} = $cadd_score_ph[$i];
			}
		}
		
		
		#fathmm-MKL 
		my %allelehash_fathmm_cs = ();
		my %allelehash_fathmm_ns = ();
		if($site_array[5] ne '.'){
			my @fathmm_alleles = split(/,/,(split(/\|/,$site_array[5]))[1]);
			my @fathmm_score_ns = split(/,/,$site_array[6]);
			my @fathmm_score_cs = split(/,/,$site_array[7]);
			for (my $i=0;$i<scalar(@fathmm_alleles);$i++){
				$allelehash_fathmm_cs{$ref.":".$fathmm_alleles[$i]} = $fathmm_score_cs[$i];
				$allelehash_fathmm_ns{$ref.":".$fathmm_alleles[$i]} = $fathmm_score_ns[$i];
			}
		}
		
		#DANN 
		my %allelehash_dann_cs = ();
		if($site_array[8] ne '.'){
			my @dann_alleles = split(/,/,(split(/\|/,$site_array[8]))[1]);
			my @dann_score_cs = split(/,/,$site_array[9]);
			for (my $i=0;$i<scalar(@dann_alleles);$i++){
				$allelehash_dann_cs{$ref.":".$dann_alleles[$i]} = $dann_score_cs[$i];
			}
		}
		
		#Funseq2
		if($site_array[10] ne '.'){
			my @funseq2_alleles = split(/,/,(split(/\|/,$site_array[10]))[1]);
			my $funseq2_score_cs = split(/,/,$site_array[11]);
			if(scalar(@funseq2_alleles) > 1){
				my @funseq2_score_cs = split(/,/,$site_array[11]);
				for (my $i=0;$i<scalar(@funseq2_alleles);$i++){
					if(!defined($allelehash_funseq2_cs{$ref.":".$funseq2_alleles[$i]})){
						#$allelehash_funseq2_cs{$ref.":".$funseq2_alleles[$i]} = log10($funseq2_score_cs[$i]);
						$allelehash_funseq2_cs{$ref.":".$funseq2_alleles[$i]} = $funseq2_score_cs[$i];
					}
				}
			}else{
				my $funseq2_score_cs = $site_array[11];
				for (my $i=0;$i<scalar(@cadd_alleles);$i++){
					if(!defined($allelehash_funseq2_cs{$ref.":".$cadd_alleles[$i]})){
						#$allelehash_funseq2_cs{$ref.":".$cadd_alleles[$i]} = log10($funseq2_score_cs);
						$allelehash_funseq2_cs{$ref.":".$cadd_alleles[$i]} = $funseq2_score_cs;
					}
				}
			}
		}
		
		#GWAS3D
		my %allelehash_gwas3d_cs = ();
		if($site_array[12] ne '.'){
			my @gwas3d_alleles = split(/,/,(split(/\|/,$site_array[12]))[1]);
			my $gwas3d_score_cs = $site_array[13];
			for (my $i=0;$i<scalar(@gwas3d_alleles);$i++){
				$allelehash_gwas3d_cs{$ref.":".$gwas3d_alleles[$i]} = $gwas3d_score_cs;
			}
		}
		
		#Funseq
		my %allelehash_funseq_cs = ();
		if($site_array[14] ne '.'){
			my @funseq_alleles = split(/,/,(split(/\|/,$site_array[14]))[1]);
			my $funseq_score_cs = $site_array[15];
			for (my $i=0;$i<scalar(@funseq_alleles);$i++){
				$allelehash_funseq_cs{$ref.":".$funseq_alleles[$i]} = $funseq_score_cs;
			}
		}
		
		#GWAVA
		my %allelehash_gwava_re = ();
		my %allelehash_gwava_tss = ();
		my %allelehash_gwava_um = ();
		if($site_array[16] ne '.'){
			my @gwava_alleles = split(/,/,(split(/\|/,$site_array[16]))[1]);
			for (my $i=0;$i<scalar(@gwava_alleles);$i++){
				$allelehash_gwava_re{$ref.":".$gwava_alleles[$i]} = $site_array[17];
				$allelehash_gwava_tss{$ref.":".$gwava_alleles[$i]} = $site_array[18];
				$allelehash_gwava_um{$ref.":".$gwava_alleles[$i]} = $site_array[19];
			}
		}
		
		#SuRFR
		my %allelehash_surfr_cs = ();
		if($site_array[20] ne '.'){
			my @sirfr_alleles = split(/,/,(split(/\|/,$site_array[20]))[1]);
			my $surfr_score_cs = $site_array[21];
			for (my $i=0;$i<scalar(@sirfr_alleles);$i++){
				$allelehash_surfr_cs{$ref.":".$sirfr_alleles[$i]} = $surfr_score_cs;
			}
		}
		
		
		for (my $i=0;$i<scalar(@alts);$i++){
			my $CADD_Cscore = defined($allelehash_cadd_cs{$ref.":".$alts[$i]}) ? $allelehash_cadd_cs{$ref.":".$alts[$i]} : ".";
			my $CADD_PHRED = defined($allelehash_cadd_ph{$ref.":".$alts[$i]}) ? $allelehash_cadd_ph{$ref.":".$alts[$i]} : ".";
			my $FunSeq_score = defined($allelehash_funseq_cs{$ref.":".$alts[$i]}) ? $allelehash_funseq_cs{$ref.":".$alts[$i]} : ".";
			my $GWAS3D_score = defined($allelehash_gwas3d_cs{$ref.":".$alts[$i]}) ? $allelehash_gwas3d_cs{$ref.":".$alts[$i]} : ".";
			my $SuRFR_score = defined($allelehash_surfr_cs{$ref.":".$alts[$i]}) ? $allelehash_surfr_cs{$ref.":".$alts[$i]} : ".";
			my $GWAVA_region_score = defined($allelehash_gwava_re{$ref.":".$alts[$i]}) ? $allelehash_gwava_re{$ref.":".$alts[$i]} : ".";
			my $GWAVA_TSS_score = defined($allelehash_gwava_tss{$ref.":".$alts[$i]}) ? $allelehash_gwava_tss{$ref.":".$alts[$i]} : ".";
			my $GWAVA_unmatched_score = defined($allelehash_gwava_um{$ref.":".$alts[$i]}) ? $allelehash_gwava_um{$ref.":".$alts[$i]} : ".";
			my $FunSeq2_score = defined($allelehash_funseq2_cs{$ref.":".$alts[$i]}) ? $allelehash_funseq2_cs{$ref.":".$alts[$i]} : ".";
			my $DANN_score = defined($allelehash_dann_cs{$ref.":".$alts[$i]}) ? $allelehash_dann_cs{$ref.":".$alts[$i]} : ".";
			my $Fathmm_MKL_score = defined($allelehash_fathmm_ns{$ref.":".$alts[$i]}) ? $allelehash_fathmm_ns{$ref.":".$alts[$i]} : ".";
			
			#print $chr."\t".$position."\t".$ref."\t".$alts[$i]."\t".$CADD_Cscore."\t".$FunSeq_score."\t".$GWAS3D_score."\t".$SuRFR_score."\t".$GWAVA_region_score."\t".$GWAVA_TSS_score."\t".$GWAVA_unmatched_score."\t".$FunSeq2_score."\t".$DANN_score."\t".$Fathmm_MKL_score."\n";
			#print OUT $chr."\t".$position."\t".$ref."\t".$alts[$i]."\t".$CADD_Cscore."\t".$CADD_PHRED."\t".$DANN_score."\t".$FunSeq_score."\t".$FunSeq2_score."\t".$GWAS3D_score."\t".$GWAVA_region_score."\t".$GWAVA_TSS_score."\t".$GWAVA_unmatched_score."\t".$SuRFR_score."\t".$Fathmm_MKL_score."\n";
			my $score_string = $CADD_Cscore."\t".$CADD_PHRED."\t".$DANN_score."\t".$FunSeq_score."\t".$FunSeq2_score."\t".$GWAS3D_score."\t".$GWAVA_region_score."\t".$GWAVA_TSS_score."\t".$GWAVA_unmatched_score."\t".$SuRFR_score."\t".$Fathmm_MKL_score;
			my @score_array = split(/\t/,$score_string);
			
			#calculate cs
			my $bf = 1;
			my $combined_p = 1;
			my $crvs = 1;
			my $prior_p = 0.5;
			my $cell_sig_p = 0.3696304;
			
			my $print_line = $chr."\t".$position."\t".$ref."\t".$alts[$i]."\t";
			foreach my $sd (@selected_tools_list){
				my $score = $score_array[$sd-1];
				$print_line .=	$score."\t"; 
				if($score eq "." || $score eq "NA"){
					$score = $mean_index{$sd};
				}
				if($sd eq "5"){
					$score = log10($score);
				}
				my $p_distribution = $selected_dis_hash{$sd}{'p_dis'};
				my $p_prob = check_distribution_probability($score,$p_distribution,'pos');
				my $n_distribution = $selected_dis_hash{$sd}{'n_dis'};
				my $n_prob = check_distribution_probability($score,$n_distribution,'neg');
				if($p_prob == 0){
					$p_prob = 0.01;
				}	
		        if($n_prob == 0){
					$n_prob = 0.01;
				}
				#print $score."\t".$p_prob."|".$n_prob."\t";
				$bf *=($p_prob/$n_prob);
				$combined_p *= ($prior_p*$p_prob)/($prior_p*$p_prob + (1-$prior_p)*$n_prob);
			}
			if(defined($CELL_ID) && $CELL_ID ne ""){
				$cell_sig_p = check_cell_type_prior($chr,$position,$CELL_ID,\%cellsg2hash);
			}
			$crvs = ($cell_sig_p * $combined_p)/$prior_p;
			$print_line .= $bf."\t".$combined_p."\t".$cell_sig_p."\t".$crvs."\n";
			print OUT $print_line;
		}
	}
}
close OUT;

sub log10 {
      my $n = shift;
      return log($n)/log(10);
}

sub check_distribution_probability{
	my ($score,$dis,$type) = @_;
	my @distribution = @{$dis};
#	my $prob = 0.1;
#	if($type eq 'pos'){
#		$prob = 0.1;
#	}
#	if($type eq 'neg'){
#		$prob = 0.9;
#	}
#	if($score eq '.' || $score eq 'NA'){
#		return $prob;
#	}
	my $start = 0;
	my $end = scalar(@distribution) -1;
	while($end - $start > 1){
		my $mid = $start + int(($end-$start)/2);
#		print $start."\t".$mid."\t".$end."\n";
		my $value = $distribution[$mid]; 
		my ($a,$b,$p) = split(/\s+/,$value);
		if($score >= $a){
			$start = $mid;
		}else{
			$end = $mid;
		}
	}
	
	my $value_s = $distribution[$start];
	my ($a_s,$b_s,$p_s) = split(/\s+/,$value_s);
	
	my $value_e = $distribution[$end];
	my ($a_e,$b_e,$p_e) = split(/\s+/,$value_e);
	
	my $distance_to_s = abs($score - $a_s);
	my $distance_to_e = abs($score - $a_e);
	if($distance_to_s >= $distance_to_e){
		return $p_e;
	}else{
		return $p_s;
	}
}

sub check_cell_type_prior{
	my ($chr,$pos,$CELL_ID,$sghash) = @_;
	my %cellsg2hash = %{$sghash};
	$chr =~ s/^chr//;
	$chr = "chr".$chr;
	my $query = $chr.":".$pos."-".$pos;
	
	my @hit_array = ();
	
	if($CELL_ID !~ /^E/){
		$DNASE = $DNASE_NE;
	}
	my $DNASE_FILE = $cell_type_annotation.$CELL_ID."-".$DNASE;
	my $DNASE_hit = 0;
	my $DNASE_score = 0;
	my $DNASE_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$DNASE})){
		my $T_DNASE = `$TABIX $DNASE_FILE $query |head -1`;
#		print "$TABIX $DNASE_FILE $query |head -1\n";
		if($T_DNASE ne '' && defined($T_DNASE)){
			$DNASE_hit = 1;
			my @array = split(/\s+/,$T_DNASE);
			my $start = $array[1];
			my $end = $array[2];
			$DNASE_score = $array[4];
			my $peak = $array[9];
			$DNASE_centrality = abs($peak-($pos-$start));
			push(@hit_array,$DNASE_hit);
		}
	}

	my $H2AZ_FILE = $cell_type_annotation.$CELL_ID."-".$H2AZ;
	my $H2AZ_hit = 0;
	my $H2AZ_score = 0;
	my $H2AZ_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H2AZ})){
		my $T_H2AZ = `$TABIX $H2AZ_FILE $query |head -1`;
		if($T_H2AZ ne '' && defined($T_H2AZ)){
			$H2AZ_hit = 1;
			my @array = split(/\s+/,$T_H2AZ);
			my $start = $array[1];
			my $end = $array[2];
			$H2AZ_score = $array[4];
			my $peak = $array[9];
			$H2AZ_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H2AZ_hit);
		}
	}
	
	my $H3K27ac_FILE = $cell_type_annotation.$CELL_ID."-".$H3K27ac;
	my $H3K27ac_hit = 0;
	my $H3K27ac_score = 0;
	my $H3K27ac_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K27ac})){
		my $T_H3K27ac = `$TABIX $H3K27ac_FILE $query |head -1`;
		if($T_H3K27ac ne '' && defined($T_H3K27ac)){
			$H3K27ac_hit = 1;
			my @array = split(/\s+/,$T_H3K27ac);
			my $start = $array[1];
			my $end = $array[2];
			$H3K27ac_score = $array[4];
			my $peak = $array[9];
			$H3K27ac_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K27ac_hit);
		}
	}
	
	my $H3K27me3_FILE = $cell_type_annotation.$CELL_ID."-".$H3K27me3;
	my $H3K27me3_hit = 0;
	my $H3K27me3_score = 0;
	my $H3K27me3_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K27me3})){
		my $T_H3K27me3 = `$TABIX $H3K27me3_FILE $query |head -1`;
		if($T_H3K27me3 ne '' && defined($T_H3K27me3)){
			$H3K27me3_hit = 1;
			my @array = split(/\s+/,$T_H3K27me3);
			my $start = $array[1];
			my $end = $array[2];
			$H3K27me3_score = $array[4];
			my $peak = $array[9];
			$H3K27me3_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K27me3_hit);
		}
	}	
	
	my $H3K36me3_FILE = $cell_type_annotation.$CELL_ID."-".$H3K36me3;
	my $H3K36me3_hit = 0;
	my $H3K36me3_score = 0;
	my $H3K36me3_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K36me3})){
		my $T_H3K36me3 = `$TABIX $H3K36me3_FILE $query |head -1`;
		if($T_H3K36me3 ne '' && defined($T_H3K36me3)){
			$H3K36me3_hit = 1;
			my @array = split(/\s+/,$T_H3K36me3);
			my $start = $array[1];
			my $end = $array[2];
			$H3K36me3_score = $array[4];
			my $peak = $array[9];
			$H3K36me3_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K36me3_hit);
		}
	}
	
	my $H3K4me1_FILE = $cell_type_annotation.$CELL_ID."-".$H3K4me1;
	my $H3K4me1_hit = 0;
	my $H3K4me1_score = 0;
	my $H3K4me1_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K4me1})){
		my $T_H3K4me1 = `$TABIX $H3K4me1_FILE $query |head -1`;
		if($T_H3K4me1 ne '' && defined($T_H3K4me1)){
			$H3K4me1_hit = 1;
			my @array = split(/\s+/,$T_H3K4me1);
			my $start = $array[1];
			my $end = $array[2];
			$H3K4me1_score = $array[4];
			my $peak = $array[9];
			$H3K4me1_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K4me1_hit);
		}
	}
	
	my $H3K4me2_FILE = $cell_type_annotation.$CELL_ID."-".$H3K4me2;
	my $H3K4me2_hit = 0;
	my $H3K4me2_score = 0;
	my $H3K4me2_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K4me2})){
		my $T_H3K4me2 = `$TABIX $H3K4me2_FILE $query |head -1`;
		if($T_H3K4me2 ne '' && defined($T_H3K4me2)){
			$H3K4me2_hit = 1;
			my @array = split(/\s+/,$T_H3K4me2);
			my $start = $array[1];
			my $end = $array[2];
			$H3K4me2_score = $array[4];
			my $peak = $array[9];
			$H3K4me2_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K4me2_hit);
		}
	}
	
	my $H3K4me3_FILE = $cell_type_annotation.$CELL_ID."-".$H3K4me3;
	my $H3K4me3_hit = 0;
	my $H3K4me3_score = 0;
	my $H3K4me3_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K4me3})){
		my $T_H3K4me3 = `$TABIX $H3K4me3_FILE $query |head -1`;
		if($T_H3K4me3 ne '' && defined($T_H3K4me3)){
			$H3K4me3_hit = 1;
			my @array = split(/\s+/,$T_H3K4me3);
			my $start = $array[1];
			my $end = $array[2];
			$H3K4me3_score = $array[4];
			my $peak = $array[9];
			$H3K4me3_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K4me3_hit);
		}
	}		

	my $H3K79me2_FILE = $cell_type_annotation.$CELL_ID."-".$H3K79me2;
	my $H3K79me2_hit = 0;
	my $H3K79me2_score = 0;
	my $H3K79me2_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K79me2})){
		my $T_H3K79me2 = `$TABIX $H3K79me2_FILE $query |head -1`;
		if($T_H3K79me2 ne '' && defined($T_H3K79me2)){
			$H3K79me2_hit = 1;
			my @array = split(/\s+/,$T_H3K79me2);
			my $start = $array[1];
			my $end = $array[2];
			$H3K79me2_score = $array[4];
			my $peak = $array[9];
			$H3K79me2_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K79me2_hit);
		}
	}
	
	my $H3K9ac_FILE = $cell_type_annotation.$CELL_ID."-".$H3K9ac;
	my $H3K9ac_hit = 0;
	my $H3K9ac_score = 0;
	my $H3K9ac_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K9ac})){
		my $T_H3K9ac = `$TABIX $H3K9ac_FILE $query |head -1`;
		if($T_H3K9ac ne '' && defined($T_H3K9ac)){
			$H3K9ac_hit = 1;
			my @array = split(/\s+/,$T_H3K9ac);
			my $start = $array[1];
			my $end = $array[2];
			$H3K9ac_score = $array[4];
			my $peak = $array[9];
			$H3K9ac_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K9ac_hit);
		}
	}
	
	my $H3K9me3_FILE = $cell_type_annotation.$CELL_ID."-".$H3K9me3;
	my $H3K9me3_hit = 0;
	my $H3K9me3_score = 0;
	my $H3K9me3_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H3K9me3})){
		my $T_H3K9me3 = `$TABIX $H3K9me3_FILE $query |head -1`;
		if($T_H3K9me3 ne '' && defined($T_H3K9me3)){
			$H3K9me3_hit = 1;
			my @array = split(/\s+/,$T_H3K9me3);
			my $start = $array[1];
			my $end = $array[2];
			$H3K9me3_score = $array[4];
			my $peak = $array[9];
			$H3K9me3_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H3K9me3_hit);
		}
	}

	my $H4K20me1_FILE = $cell_type_annotation.$CELL_ID."-".$H4K20me1;
	my $H4K20me1_hit = 0;
	my $H4K20me1_score = 0;
	my $H4K20me1_centrality = -1;
	if(defined($cellsg2hash{$CELL_ID."-".$H4K20me1})){
		my $T_H4K20me1 = `$TABIX $H4K20me1_FILE $query |head -1`;
		if($T_H4K20me1 ne '' && defined($T_H4K20me1)){
			$H4K20me1_hit = 1;
			my @array = split(/\s+/,$T_H4K20me1);
			my $start = $array[1];
			my $end = $array[2];
			$H4K20me1_score = $array[4];
			my $peak = $array[9];
			$H4K20me1_centrality = abs($peak-($pos-$start));
			push(@hit_array,$H4K20me1_hit);
		}
	}
	
	
	my $prior = 1/(1+exp(-(-0.5339527052 + 1.0513562209 * $H3K4me1_hit + 1.5659681399 * $H3K36me3_hit + 1.2131942069 * $DNASE_hit + 0.9750312605 * $H3K79me2_hit + -0.4843821400 * $H3K9me3_hit + 1.5150317212 * $H3K27me3_hit + 0.0008691201 * $H3K4me2_score + 0.0003089830 * $H3K4me3_score + 0.0043517819 * $H3K36me3_score + -0.0001497833 * $H3K79me2_centrality)));
	if($prior < 0.3696304){
		$prior = 0.3696304;
	}
	return $prior;
}

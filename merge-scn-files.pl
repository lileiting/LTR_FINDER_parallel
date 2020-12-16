#!/usr/bin/env perl

use warnings;
use strict;

my $coord_adj = 0;
my $harvest_format = 1;
my $annotator = "LTR_FINDER_parallel";
my $finder_para = "-w 2 -C -D 15000 -d 1000 -L 7000 -l 100 -p 20 -M 0.85";

open GFF, "> out.finder.combine.gff3" or die $!;
print GFF "##gff-version   3\n";
open my $out, "> out.finder.combine.scn" or die $!;

#print out headers
if ($harvest_format == 0){
    print $out "index	SeqID	Location	LTR len	Inserted element len	TSR	PBS	PPT	RT	IN (core)	IN (c-term)	RH	Strand	Score	Sharpness	Similarity\n";
} else {
print $out "#LTR_FINDER_parallel -seq \$seq_file -size \$size_of_each_piece -time \$timeout -try1 \$try1 -harvest_out -threads \$threads -cut \$cut -finder \$ltr_finder
# LTR_FINDER args=$finder_para
# LTR_FINDER_parallel version=\$version
# predictions are reported in the following way
# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr
# where:
# s = starting position
# e = ending position
# l = length
# ret = LTR-retrotransposon
# lLTR = left LTR
# rLTR = right LTR
# sim = similarity
# seq-nr = sequence order\n";
}

my $count = 1;
open my $list_scn_fh, "< list-of-scn.txt" or die $!;
while( my $scn_file = <$list_scn_fh> ){
	chomp $scn_file;
	open my $scn_fh, $scn_file or die $!;
	while (<$scn_fh>){
		next unless /^\[/;
		s/^\[\s+[0-9]+\]/\[NA\]/;
		my ($index, $id, $loc, $ltr_len, $ele_len, $TSR, $PBS, $PPT, $RT, $IN_core, $IN_cterm, $RH, $strand, $score, $sharpness, $sim) = (split);
                my $base = $id;
                my $order = $base;
                $order =~ s/^\D+//;
                $order--; 

		#convert coordinates back to the genome scale
		my @coord = ($loc, $PBS, $PPT, $RT, $IN_core, $IN_cterm, $RH);
		my $i = -1;
		foreach (@coord) {
			$i++;
			next if /^N-N$/;
			my ($start, $end) = ($1+$coord_adj, $2+$coord_adj) if /([0-9]+)\-([0-9]+)/;
			$coord[$i] = "$start-$end";
			}
		my ($start, $end) = ($1, $2) if $coord[0] =~ /([0-9]+)\-([0-9]+)/;
		my ($lltr, $rltr) = ($1, $2) if $ltr_len=~/([0-9]+),([0-9]+)/; 
		my ($lltr_e, $rltr_s) = ($start+$lltr-1, $end-$rltr+1);

		#output LTRharvest or LTR_FINDER (-w 2) format
		if ($harvest_format == 0){
			print $out "[NA]\t$base\t$coord[0]\t$ltr_len\t$ele_len\t$TSR\t$coord[1]\t$coord[2]\t$coord[3]\t$coord[4]\t$coord[5]\t$coord[6]\t$strand\t$score\t$sharpness\t$sim\n";
			} else {
			$sim*=100;
			print $out "$start $end $ele_len $start $lltr_e $lltr $rltr_s $end $rltr $sim $order $base\n";
			}

		#print GFF format
		my $chr = $base;
		print GFF "$chr\t$annotator\trepeat_region\t$start\t$end\t.\t$strand\t.\tID=repeat_region$count\n";
		#print GFF "$chr\t$annotator\ttarget_site_duplication\t$lTSD\t.\t$strand\t.\tParent=repeat_region$count\n" unless $TSD eq "NA";
		print GFF "$chr\t$annotator\tLTR_retrotransposon\t$start\t$end\t.\t$strand\t.\tID=LTR_retrotransposon$count;Parent=repeat_region$count;tsd=$TSR;ltr_identity=$sim;seq_number=$order\n";
		print GFF "$chr\t$annotator\tlong_terminal_repeat\t$start\t$lltr_e\t.\t$strand\t.\tParent=LTR_retrotransposon$count\n";
		print GFF "$chr\t$annotator\tlong_terminal_repeat\t$rltr_s\t$end\t.\t$strand\t.\tParent=LTR_retrotransposon$count\n";
		#print GFF "$chr\t$annotator\ttarget_site_duplication\t$rTSD\t.\t$strand\t.\tParent=repeat_region$count\n" unless $TSD eq "NA";
		print GFF "###\n";

		}
	$count++;
        close $scn_fh;
	}
close $list_scn_fh;
close GFF;

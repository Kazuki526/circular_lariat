#!/usr/bin/perl
use warnings;
use strict;

my $pwd = `pwd`;chomp $pwd;
if($pwd ne "$ENV{HOME}/Dropbox/work/segdup"){die "ERROR:: please work on /Dropbox/work/segdup\n";}

my $genome = "/Volumes/areca42TB/GRCh38.d1.vd1.fa";
(-e $genome) or die "ERROR::there are not exist $genome\n";

if(scalar(@ARGV) != 2){die "ERROR::designate input file & output dir\n";}
my ($in_file, $out_dir) = @ARGV;

open(IN,"nkf -L $in_file|");
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		if($line[0] eq ""){next;}
		my $group_dir = "group$line[0]";
		if(!-e "$group_dir"){mkdir "$out_dir/$group_dir";}
		if($line[2] ==0){next;}
		my ($chr,$start,$end)=("",0,0);
		if($line[3] =~/^(chr.+):(\d+)-(\d+)$/){($chr,$start,$end)=($1,$2,$3);}else{die "ERROR::what region $line[3]\n";}
		my $width= $end -$start+1;
		my $reserve=int($width*0.2);
		$start-=$reserve;$end+=$reserve;
		`samtools faidx $genome $chr:$start-$end >$out_dir/$group_dir/$line[1]_$line[3].fa`;
		&pairwise_search($line[3],"$out_dir/$group_dir/$line[1]_$line[3]_pairwise.tsv");
}
close IN;
exit;

sub pairwise_search ( $ $ ){
		my($region,$out_file)=@_;
		my($focal_chr,$focal_start,$focal_end);
		if($region =~ /^(chr.*):(\d+)-(\d+)$/){
				$focal_chr = $1;
				$focal_start = $2;
				$focal_end = $3;
		}else{die "ERROR::what region => $region ??\n";}
		open(PAIRW,"nkf -L Pairwise_SDs_Human_hg38.txt|");
		open(OUT,">$out_file");
		<PAIRW>;
		while(my $pairw_line = <PAIRW>){
				chomp $pairw_line;
				my @pairw_line = split(/  /,$pairw_line);
				if((($focal_chr eq $pairw_line[1]) && 
						((($focal_start <=$pairw_line[2])&&($focal_end >=$pairw_line[2])) || (($focal_start <=$pairw_line[3])&&($focal_end >=$pairw_line[3])))) ||
				   (($focal_chr eq $pairw_line[5]) &&
						((($focal_start <=$pairw_line[6])&&($focal_end >=$pairw_line[6])) || (($focal_start <=$pairw_line[7])&&($focal_end >=$pairw_line[7]))))){
						print OUT join("\t",@pairw_line) . "\n";
				}
		}
		close PAIRW;
		close OUT;
}





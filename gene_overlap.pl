#!/usr/bin/perl
use strict;
use warnings;

#read focal region
my %focal_list = ();
my %region2id = ();
open(IN,"nkf -L group/faidx_region.tsv|") or die "ERROR:problem in faidx_region.tsv\n";
my @id=();
while(<IN>){
		chomp;
		if($_ !~ /\S/){next;}
		my @line=split(/\t/,);
		push(@id,$line[0]);
		for(my$i=0;$i<@line;$i++){$line[$i] =~ s/\s//g;} #delete space
		$region2id{$line[1]}{$line[2]}{end}=$line[3];
		$region2id{$line[1]}{$line[2]}{id}=$line[0];
		$region2id{$line[1]}{$line[2]}{former}=0;
		$focal_list{$line[0]}{region}="$line[1]\t$line[2]\t$line[3]";
}
close IN;

#read GFF
open(GFF,"gunzip -c /Volumes/areca42TB/Homo_sapiens.GRCh38.90.gff3.gz|") or die "ERROR::cannot open GFF file\n";
my $bef_gene="";
my %gene=();
while(<GFF>){
		if($_ =~ /^#/){next;}
		chomp;
		my @line = split(/\t/,);
		if(($line[2] ne "gene")||($line[8] !~ /;biotype=protein_coding;/)){next;}
		my $chr="chr$line[0]";
		my $gene;
		if($line[8] =~/;Name=([^;]+);/){$gene=$1;}
		my($gene_start,$gene_end)=($line[3],$line[4]);
		while($gene_start =~ s/^(-?\d+)(\d\d\d)/$1,$2/){;}
		while($gene_end =~ s/^(-?\d+)(\d\d\d)/$1,$2/){;}
		$gene{$gene}="$chr:$gene_start-$gene_end";
		foreach my $start(sort {$a <=> $b} keys%{$region2id{$chr}}){
				my ($id,$end)=($region2id{$chr}{$start}{id},$region2id{$chr}{$start}{end});
				if((($start>=$line[3])&&($start<=$line[4]))||
						(($end>=$line[3])&&($end<=$line[4]))){
						$focal_list{$id}{overlap}.="$gene;";
						if(!defined$focal_list{$id}{former}){
								$focal_list{$id}{former}=$bef_gene;
						}
				}elsif(($region2id{$chr}{$start}{former}==0)&&($start<$line[3])){
						$focal_list{$id}{back}=$gene;
						if(!defined$focal_list{$id}{former}){
								$focal_list{$id}{former}=$bef_gene;
						}
						$region2id{$chr}{$start}{former}++;
				}
		}
		$bef_gene=$gene;
}
close GFF;

open(OUT,">gene_overlap.tsv");
my $bef_group="";
foreach my $id(@id){
		my$group;
		if($id=~/^([LC]\d+)/){$group=$1;}
		if($bef_group ne $group){print OUT "\n";$bef_group=$group;}
		my $overlap = "\t\t";
		if(defined$focal_list{$id}{overlap}){
				my@overlap=split(/;/,$focal_list{$id}{overlap});
				my@overlap_posi=();
				foreach my $gene(@overlap){push(@overlap_posi,$gene{$gene});}
				$overlap=join(",",@overlap)."\t".join(",",@overlap_posi)."\t";
		}
		my($former,$back)=($focal_list{$id}{former},$focal_list{$id}{back});
		print OUT "$id\t$focal_list{$id}{region}\t\t$overlap\t$former\t$gene{$former}\t\t$back\t$gene{$back}\n";
}
close OUT;

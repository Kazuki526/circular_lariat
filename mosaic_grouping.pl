#!/usr/bin/perl
use warnings;
use strict;

my %index = (); #=>group_index
my %group_index = (); #=>SD
my %SD = (); #=>index

my $mosaicSD_file ="MosaicSDs_Human_hg38.txt";
open(IN,"$mosaicSD_file") or die "ERROR:cannot open $mosaicSD_file\n";
<IN>; #header line
my $new_group_index=1;
while(<IN>){
		chomp;
		my @line=split(/  /,);
		if(scalar(@line) <=4){next;}
		$line[2] =~ s/ :$//;
		my $sd="$line[0]:$line[1]-$line[2]";
		my @group_index_list=();
		for(my $i=3; $i < scalar(@line); $i++){
				my $id=$line[$i];$id=~s/^-//;
				if(defined $index{$id}){
						if(grep{$_ == $index{$id}}@group_index_list){next;
					}else{push(@group_index_list, $index{$id});}
				}
		}
		if(scalar(@group_index_list)==0){
				for(my $i=3; $i < scalar(@line); $i++){
						my $id=$line[$i];$id=~s/^-//;
						$index{$id}=$new_group_index;
				}
				$group_index{$new_group_index}="$sd";
				$SD{$sd}=join("\t",@line[3..$#line]);
				$new_group_index++;
		}elsif(scalar(@group_index_list)==1){
				my $merge_gi=$group_index_list[0];
				for(my $i=3; $i < scalar(@line); $i++){
						my $id=$line[$i];$id=~s/^-//;
						$index{$id}=$merge_gi;
				}
				$group_index{$merge_gi}.="\t$sd";
				$SD{$sd}=join("\t",@line[3..$#line]);
		}else{
				@group_index_list = sort{$a<=>$b}@group_index_list;
				my $merge_gi=shift(@group_index_list);
				foreach my $gi(@group_index_list){
						my @sd_list = split(/\t/,$group_index{$gi});
						foreach my $merge_sd(@sd_list){
								my @index_list=split(/\t/,$SD{$merge_sd});
								foreach my $id(@index_list){
										$id =~s/^-//;
										$index{$id}=$merge_gi;
								}
						}
						$group_index{$merge_gi}.="\t$group_index{$gi}";
						$group_index{$gi}="";
				}
				for(my $i=3; $i < scalar(@line); $i++){
						my $index=$line[$i];$index=~s/^-//;
						$index{$index}=$merge_gi;
				}
				$group_index{$merge_gi}.="\t$sd";
				$SD{$sd}=join("\t",@line[3..$#line]);
		}
}
close IN;


open(OUT,">mosaic_group.tsv");
my$l=0;
foreach my $gi(keys %group_index){
		my @sd_list=split(/\t/,$group_index{$gi});
		if(scalar(@sd_list)<=1){next;}
		$l++;
		if(scalar(@sd_list)>10){print scalar(@sd_list)."\n";}
		foreach my $sd(@sd_list){
				print OUT "$sd\t$SD{$sd}\n";
		}
		print OUT "\n";
}
close OUT;
print "$l\n";

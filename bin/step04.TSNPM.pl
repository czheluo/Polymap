#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin1,$fin2,$fout,$min);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"int:s"=>\$fin1,
	"out:s"=>\$fout,
			) or &USAGE;
&USAGE unless ($fout);

open IN,$fin1;
my %idg;
my $nid;
while (<IN>) {
	chomp;
	if (/^ID/){
		my @nid=split(/\,/,$_);
		$nid= (scalar @nid) -3;
		#print Dumper $nid;die;
	}else{
		my ($idg1,$geno)=split(/\,/,$_,2);
		$idg{$idg1}=$geno;
	}
	#print Dumper $geno;die;
}
close IN;

my @tsnpm=glob("$out/step4.qtlanalysis/TetraploidSNPMap_QTLfiles_pp/*.txt");
foreach my $tsnp (@tsnpm) {
	my $fileID=split(/\./,basename($tsnp));
	open In,$tsnp;
	open Out,">$out/step4.qtlanalysis/TetraploidSNPMap_QTLfiles_pp/$fileID.SNPloc";
	my $n1=<In>;
	chomp($n1);
	#print Dumper $n1;
	#my %id;
	print OUT "$nid\t$n1\n";
	while (<In>) {
	chomp;
	my($id1,undef)=split/\s+/,$_;
	my @geno=split/\,/,$idg{$id1};
	my @out;
	for (my $i=0;$i<@geno ;$i++) {
		if ($geno[$i] eq "NA") {
			$geno[$i]=~s/NA/9/g;
			push @out,join("\t",$geno[$i]);
		}else{
			push @out,join("\t",$geno[$i]);
		}
	}
	print OUT "$id1\n @out\n";
	}
	close In;
	close Out;
}

my @tsnpm=glob("$out/step4.qtlanalysis/TetraploidSNPMap_QTLfiles_P1/*.txt");
foreach my $tsnp (@tsnpm) {
	my $fileID=split(/\./,basename($tsnp));
	open In,$tsnp;
	open Out,">$out/step4.qtlanalysis/TetraploidSNPMap_QTLfiles_P1/$fileID.SNPloc";
	my $n1=<In>;
	chomp($n1);
	#print Dumper $n1;
	#my %id;
	print OUT "$nid\t$n1\n";
	while (<IN>) {
	chomp;
	my($id1,undef)=split/\s+/,$_;
	my @geno=split/\,/,$idg{$id1};
	my @out;
	for (my $i=0;$i<@geno ;$i++) {
		if ($geno[$i] eq "NA") {
			$geno[$i]=~s/NA/9/g;
			push @out,join("\t",$geno[$i]);
		}else{
			push @out,join("\t",$geno[$i]);
		}
	}
	print OUT "$id1\n @out\n";
	}
	close In;
	close Out;
}

my @tsnpm=glob("$out/step4.qtlanalysis/TetraploidSNPMap_QTLfiles_P2/*.txt");
foreach my $tsnp (@tsnpm) {
	my $fileID=split(/\./,basename($tsnp));
	open In,$tsnp;
	open Out,">$out/step4.qtlanalysis/TetraploidSNPMap_QTLfiles_P2/$fileID.SNPloc";
	my $n1=<In>;
	chomp($n1);
	#print Dumper $n1;
	#my %id;
	print OUT "$nid\t$n1\n";
	while (<In>) {
	chomp;
	my($id1,undef)=split/\s+/,$_;
	my @geno=split/\,/,$idg{$id1};
	my @out;
	for (my $i=0;$i<@geno ;$i++) {
		if ($geno[$i] eq "NA") {
			$geno[$i]=~s/NA/9/g;
			push @out,join("\t",$geno[$i]);
		}else{
			push @out,join("\t",$geno[$i]);
		}
	}
	print OUT "$id1\n @out\n";
	}
	close In;
	close Out;
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:

	eg: perl -int filename -out filename 
	

Usage:
  Options:
	-int input genotype file name
	-out ouput file name SNPloc
	-h         Help

USAGE
        print $usage;
        exit;
}

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$pdep,$gdep,$ms,$fdat,$fdep,$pid,$mid);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$fin,
	"out:s"=>\$fout,
	"pid:s"=>\$pid,
	"mid:s"=>\$mid,
	"pdep:s"=>\$pdep,
	"gdep:s"=>\$gdep,
	"dat:s"=>\$fdat,
	"dep:s"=>\$fdep,
	"map"s"=>\$map,
	"ms:s"=>\$ms,
			) or &USAGE;
&USAGE unless ($fout);

mkdir $fout if (!-d $fout);
my $updog="$out/updog";
mkdir $updog if (!-d $updog);
my $RData="$out/RData";
mkdir $RData if (!-d $RData);
$fout=ABSOLUTE_DIR($fout);
$pdep||=10;
$gdep||=10;
$ms||=0;
my @idr;
my $npid=0;
my $nmid=0;
open IN,$fin;
open Out,">$fout/dosage.matrix";
open Dat,">$fdat";
open Dep,">$fdep";
open Map,">$map";
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	if (/^#/) {
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$inform)=split/\s+/,$_,10;
		@idr=split/\s+/,$inform;
		for (my $i=0;$i<@idr ;$i++) {
			if ($idr[$i] eq $pid) {
				$npid=$i+1;
			}elsif ($idr[$i] eq $mid) {
				$nmid=$i+1;
			}else{
				next;
				}
		}
		print Out "$id\t$inform\n";
		print Dat "id\tsnp\tcounts\tsize\n";
		print Dep "$id\t$inform\n";
		print Map "Chr\tPos\n";
	}else {
		next if (/^sca/);
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$inform)=split/\s+/,$_,10;
		my @inform=split/\s+/,$inform;
		my @f=split/\:/,$inform[$npid];
		my @m=split/\:/,$inform[$nmid];
		if ($f[0] eq "./." || $f[0] eq "./././." || $f[0] eq "./././././." ) {
				next;
			}elsif ($m[0] eq "./." ||  $m[0] eq "./././." || $m[0] eq "./././././.") {
				next;
			}
		next if ( $f[2] < $pdep || $m[2] < $pdep);
		my $k=1;
		my $sm=0;
		my $total=scalar(@inform);
		my @out;
		my @dep;
		#print Dumper $total;die;
		for (my $i=0;$i<@inform ;$i++) {
			my @ad=split(/\:/,$inform[$i]);
			if ($ad[0] eq "./." || $ad[0] eq "./././." || $ad[0] eq "./././././.") {
				$sm++;
				next;
			}
			if ($ad[2]<$gdep) {
				$k=2;
				#last;
			}
		}
		my $miss=$sm/$total;
		next if ($k==2);
		next if ($miss>$ms);
		for (my $i=0;$i<@inform ;$i++) {
			my @ad=split(/\:/,$inform[$i]);
			if ($ad[0] eq "./././.") {
				push @out,join("\t",qw(NA));
				push @dep,join("\t",qw(NA));
				}else{
					my @geno=split/\//,$ad[0];
					my $gm=0;
					foreach (@geno) {
						if ($_ eq "0") {
							next;
						}else{
							$gm++;
						}
						}
				#print Dumper @geno;die;
				push @out,join("\t",$gm);
				push @dep,join("\t",$ad[2]);
			}
		}
		#print Dumper @out;die;
		print Out join("\t",$id,@out),"\n";
		print Dep join("\t",$id,@dep),"\n";
		print Map "$chrom\t$pos\n";
		for (my $i=0;$i<@idr ;$i++) {
			if ($inform[$i] =~ /\.\/./){
				print Dat "$idr[$i]\t$id\tNA\tNA\n";
				}else{
					my @format=split/\:/,$inform[$i];
					my ($ad,undef,undef)=split/\,/,$format[1];
					print Dat "$idr[$i]\t$id\t$ad\t$format[2]\n";
			}
		}
	}
}
close IN;
close Out;
close Dat;
close Dep;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:

	perl vcf2filter.pl -int pop.final.vcf -out pop.vcf -dep 5 -ms 0.2

Usage:
  Options:
	-int input filename
	-out output filename
	-pid	paternal ID
	-mid	maternale ID
	-pdep parents depth 
	-gdep generation depth 
	-ms missing rate default was 0
	-dat output snpdat filename
	-dep output all sample depth filename
	-map output doasge matrix map {file names};
 	-h         Help

USAGE
        print $usage;
        exit;
}


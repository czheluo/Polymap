#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$dOut,$adjust);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"out:s"=>\$dOut,
	"adjust:s"=>\$adjust,
			) or &USAGE;
&USAGE unless ($dmap and $dOut);

my @map=glob("$dmap/*.sexAver.map");
open Out,">$dOut/total.sexAver.map";
#open Draw,">$dOut/total.sexAver.map.draw";
my %info;
foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	my $nloc=`wc -l $map`;
	$nloc=(split(/\s+/,$nloc))[0];
	next if ($nloc <= 2);
	print Out "group\t$lgID\n";
	open In,$map;
	my $max=0;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		my ($id,$pos)=split(/\s+/,$_);
		$info{$id}{lgID}=$lgID;
		$info{$id}{pos}=$pos;
		$max=$pos if ($pos > $max);
	}
	close In;
	my $newdis=$max;
	if ($adjust) {
		$newdis=rand(80)+120;
	}
	my $n=0;
	my $pos0;
	foreach my $id (sort {$info{$a}{pos}<=>$info{$b}{pos}} keys %info) {
		next if ($info{$id}{lgID} ne $lgID);
		$info{$id}{pos}=$info{$id}{pos}*$newdis/$max;
		if ($n == 0) {
			$pos0=$info{$id}{pos};
			$n++;
		}
		$info{$id}{pos}=$info{$id}{pos}-$pos0;
		print Out $id,"\t",$info{$id}{pos},"\n";
	}

}
close Out;
@map=glob("$dmap/*.male.map");
open Out,">$dOut/total.male.map";
my %males;

foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	my $nloc=`wc -l $map`;
	$nloc=(split(/\s+/,$nloc))[0];
	next if ($nloc <= 2);
	print Out "group\t$lgID\n";
	open In,$map;
	my $max=0;
	my %male;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		my ($id,$pos)=split(/\s+/,$_);
		$male{$id}{lgID}=$lgID;
		$male{$id}{pos}=$pos;
		$max=$pos if ($pos > $max);
	}
	close In;
	my $newdis=$max;
	if ($adjust) {
		$newdis=rand(80)+120;
	}
	my $n=0;
	my $pos0;
	foreach my $id (sort {$male{$a}{pos}<=>$male{$b}{pos}} keys %male) {
		next if ($info{$id}{lgID} ne $lgID);
		$male{$id}{pos}=$male{$id}{pos}*$newdis/$max;
		if ($n == 0) {
			$pos0=$male{$id}{pos};
			$n++;
		}
		$male{$id}{pos}=$male{$id}{pos}-$pos0;
		$males{$id}{pos}=$male{$id}{pos};
		print Out $id,"\t",$male{$id}{pos},"\n";
	}

}
close Out;
@map=glob("$dmap/*.female.map");
open Out,">$dOut/total.female.map";
	my %females;

foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	my $nloc=`wc -l $map`;
	$nloc=(split(/\s+/,$nloc))[0];
	next if ($nloc <= 2);
	print Out "group\t$lgID\n";
	open In,$map;
	my $max=0;
	my %female;

	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		my ($id,$pos)=split(/\s+/,$_);
		$female{$id}{lgID}=$lgID;
		$female{$id}{pos}=$pos;
		$max=$pos if ($pos > $max);
	}
	close In;
	my $newdis=$max;
	if ($adjust) {
		$newdis=rand(100)+120;
	}
	my $n=0;
	my $pos0;
	foreach my $id (sort {$female{$a}{pos}<=>$female{$b}{pos}} keys %female) {
		next if ($info{$id}{lgID} ne $lgID);
		$female{$id}{pos}=$female{$id}{pos}*$newdis/$max;
		if ($n == 0) {
			$pos0=$female{$id}{pos};
			$n++;
		}
		$female{$id}{pos}=$female{$id}{pos}-$pos0;
		$females{$id}{pos}=$female{$id}{pos};
		print Out $id,"\t",$female{$id}{pos},"\n";
	}
	
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -dmap	<file>	input map dir name 
  -out	<file>	output dir 
  -h         Help

USAGE
        print $usage;
        exit;
}

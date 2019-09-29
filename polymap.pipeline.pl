#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$popt,$out,$bin,$pid,$mid,$nchr,$ref,$step);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"pid:s"=>\$pid,
	"mid:s"=>\$mid,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($vcf and $popt and $out and $pid and $mid and $nchr);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
my $dsh="$out/work_sh";
mkdir $dsh if (!-d $dsh);
open Log,">$out/work_sh/polymap.$BEGIN_TIME.log";
$step||=1;
if ($step==1) {
	open SH1,">$dsh/01.vcf-convert.sh";
	print SH1 "perl $Bin/bin/vcf2dosage.pl -vcf $vcf -out $out/01.vcf-convert -pid $pid -mid $mid -dat $out/01.vcf-convert/snpdata.dat -dep $out/01.vcf-convert/dosage.matrix.dep -map $out/01.vcf-convert/dosage.matrix.map &&";
	print SH1 "Rscript $Bin/bin/step01.filter.R --out $out --pid $pid --mid $mid "
	close SH1;
	my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=10G --CPU 3  $dsh/01.vcf-convert.sh";
	`$job`;
	$step++ if ($step ne $stop);
}
if ($step == 2) {
	open SH2,">$dsh/02.linkage.sh";
	print SH2 "Rscript $Bin/bin/step02.linkage.R --pid $pid --mid $mid ";
	print SH2 "perl "
	close SH2;
	my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=10G --CPU 3  $dsh/02.linkage.sh";
	`$job`;
	$step++ if ($step ne $stop);
}
if ($step == 3) {
	open SH3,">$dsh/03.evaluating.map.sh";
	print SH3 "Rscript $Bin/bin/step03.evalu.map.R --out $out && ";
	print SH3 "perl $Bin/bin/collinear.pl -dmap $out/step3.evaluating.map -out $out/step3.evaluating.map -adjust && "
	print SH3 "perl $Bin/bin/drawAligmentRalationMap.pl -m $out/step3.evaluating.map/total.female.map -k $out/step3.evaluating.map/total.female -o $out/step3.evaluating.map && "
	print SH3 "perl $Bin/bin/drawAligmentRalationMap.pl -m $out/step3.evaluating.map/total.male.map -k $out/step3.evaluating.map/total.male -o $out/step3.evaluating.map && "
	print SH3 "perl $Bin/bin/drawAligmentRalationMap.pl -m $out/step3.evaluating.map/total.sexAver.map -k $out/step3.evaluating.map/total.sexAver -o $out/step3.evaluating.map"
	close SH3;
	my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=3G --CPU 1  $dsh/03.evaluating.map.sh";
	`$job`;
	$step++ if ($step ne $stop);
}
if ($step == 4) {
	open SH4,">$dsh/04.qtlanalysis.sh";
	print SH4 "Rscript $Bin/bin/step04.qtlanalysis.R --out $out && ";
	print SH4 "perl $Bin/bin/step04.TSNPM.pl -int $out/01.vcf-convert/dosage.matrix -out $out/step4.qtlanalysis "
	close SH4;
	my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=3G --CPU 1  $dsh/04.qtlanalysis.sh";
	`$job`;
	$step++ if ($step ne $stop);
}
if ($step == 5) {
	open SH5,">$dsh/05.output.stat.sh";
	print SH5 "Rscript $Bin/bin/step05.output.stat.R --out $out &&";
	close SH5;
	my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=3G --CPU 1  $dsh/05.output.stat.sh";
	`$job`;
	$step++ if ($step ne $stop);
}
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
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -vcf	<file>	input file name
  -out	<dir>	output dir of filename
  -pid	<str>	paternal id
  -mid	<str>	maternal id
	-step which step you want 
	-stop control the steps 

  -h         Help

USAGE
        print $usage;
        exit;
}

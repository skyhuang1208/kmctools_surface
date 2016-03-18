#! /usr/bin/perl -w
#############################################################################################
# make a .ltcp or .xyz for only solute and vacancies
# By Sky, 2015
#############################################################################################

use strict;
$|=1;  # don't buffer output

if (scalar(@ARGV) != 1) 
{
        print "Make a conf file with no solvents\n";
        print "ABORT:: parameters required!\a\n";
        print "Usage: $0 [<t0.xyz>||<t0.ltcp>]\n";
	die;
}

my ($history_file) = @ARGV;

if (! -e $history_file)	{ die("ABORT::history file <$history_file> doesn't exist!\n");}
else			{ open(IN, $history_file)   || die ("Cannot open t0.xyz file for reading\n");}

my $line=0;
my $n= 0;
my @data;
while (my $buff=<IN>){
	$line ++;

	next if($line<3);

	my ($type, $xi, $yi, $zi) = split(/\s+/, $buff);

	if($type == -1){
		$n ++;
		$data[$n][0]= $type;
		$data[$n][1]= $xi;
		$data[$n][2]= $yi;
		$data[$n][3]= $zi;
	}
}
$n ++;
print "$n\n\n";
for(my $i= 1; $i<= $n-1; $i ++){
	my $el;

	if(-1==$data[$i][0]){
		$el= 'B';
	}
	elsif(0==$data[$i][0]){
		$el= 'V';
	}
	else{
		die("strange element $data[$i][0]\n");
	}
	print "B $data[$i][1] $data[$i][2] $data[$i][3]\n";
}
print "V -1 -1 -1\n";

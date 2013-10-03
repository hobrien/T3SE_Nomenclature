#!/usr/bin/perl -w

=head1 NAME

MakeDistPlot.pl
9 Nov 2012

=head1 SYNOPSIS

MakeDistPlot.pl -m <mlst_distances> -t <t3se_distances> [ -f <family_name> -o <out_file> ]


=head1 DESCRIPTION

Matches distance data output from Mega (columns, csv) for a T3SE and for MLST genes
and produces a scatterplot of distances, with a line indicating unity (equal distance
between T3SE and MLST genes).

If T3SE family name is not specified, it is guessed from the file name
If outfile is not specified, default filename is created ( <t3se_name>_dist.pdf )

=head2 NOTES

Determinations are also made if strains are in the same phylogroup and if sequences
have been assigned to the same subfamily. This information is not currently included
in the plot, but would be easy to add

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use autodie qw(open close);

my $mlst_file;
my $t3se_file;
my $t3se_fam;
my $outfile;

GetOptions( 'mlst_file=s' => \$mlst_file,
	        't3se_file=s' => \$t3se_file,
	        'family_name=s' => \$t3se_fam,
	        'outfile=s' => \$outfile
	      )
or die "failure to communicate\n";

unless ( $mlst_file and -e $mlst_file ) {
  die "please select valid mlst distance file file\n";
}
unless ( $t3se_file and -e $t3se_file ) {
  die "please select valid T3SE distance file file\n";
}

#extract file base name for use in default filenames
(my $name, my $path, my $ext) = fileparse($t3se_file, qr/\.[^.]*/);

unless ( $t3se_fam ) {
  if ( $name =~ /((hop[A-Z]{1,3})\d?|(avr[A-Z]{1,3})\d?|(hrp[A-Z])\d?)/ ) {
    $t3se_fam =$1;
  }
  else {
    die "T3SE family could not be parsed from file name. Please specify\n";
  }
}

unless (  $outfile ) {
  $outfile = $t3se_fam . '_dist.pdf';
}

#Read MLST distances from file
my %mlst_dist;

open(my $mlst_fh, $mlst_file);
while (<$mlst_fh>) {
  chomp;
  if ($_ =~/^Species/ or $_ =~ /^Table/ or $_ =~ /^The number/ or $_ =~ /^Disclaimer/ or $_ =~ /^\d\. / or $_ =~ /^\s*$/) { next;}
  my @fields = split(/,/, $_);
  if ( abs($fields[2]) > 5 ) { $fields[2] = 5; }
  $mlst_dist{join(",", sort(@fields[0,1]))} = $fields[2];
}
close($mlst_fh);

#Read T3SE distances, combine with corresponding MLST distances and write to temp file
open(my $dist_fh, ">distances.txt");
print $dist_fh "Species1,Species2,$t3se_fam,MLST,Phylo,SubFam\n";
open(my $t3se_fh, $t3se_file);
while (<$t3se_fh>) {
  chomp($_);
  if ($_ =~ /^Species/ or $_ =~ /^Table/ or $_ =~ /^The number/ or $_ =~ /^Disclaimer/ or $_ =~ /^\d\. / or $_ =~ /^\s*$/) { next;}
  if ($_ =~ / FS/ or $_ =~ / N/ or $_ =~ / C/ or $_ =~ / Y/ or $_ =~ /partial/ ) { next; }
  my @fields = split(/,/, $_);
  my $strain1 = $fields[0];
  my $strain2 = $fields[1];
  $strain1 =~ s/ (.*)//;
  my $fam1 = $1;
  $strain2 =~ s/ (.*)//;
  my $fam2 = $1;
  my $key = join(",", sort($strain1, $strain2));
  if ( $strain1 eq $strain2 ) { $mlst_dist{$key} = 0; }
  unless ( exists($mlst_dist{$key}) ) { next; } 
  if ( abs($fields[2]) > 5 ) { $fields[2] = 5; }
  my $phylocomp = GetPhylogroup($key);
  unless ( $phylocomp ) { next; }
  pop(@fields);
  print $dist_fh join(",", (@fields, $mlst_dist{$key},$phylocomp,GetSubFamily($fam1,$fam2))), "\n";
}
close($t3se_fh);
close($dist_fh);

#Write R script to file and execute
open(my $r_script, ">make_plot.R");
print $r_script "x<-read.csv(file='distances.txt', header=T)\n";
print $r_script "library(ggplot2)\n";
print $r_script "pdf(file='$outfile')\n";
print $r_script "ggplot(x, aes(MLST, $t3se_fam))+geom_point()+geom_abline()+theme(
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.background = element_blank(),
     axis.line=element_line())\n";
print $r_script "dev.off()\n";
close($r_script);
print `R CMD BATCH make_plot.R`;

#clean up
print `rm make_plot.R distances.txt make_plot.Rout`;

#Compare strain names and determine if comparison is within a phylogroup or between
#This info is written to the infile for R, but is not currently plotted. This can be easily
#added to the R script
sub GetPhylogroup {
  my %phylogroups = (
    "Pac302273" => 2,
    "Pan302091" => 1,
    "Pan7286" => 1,
    "Pan3739" => 1,
    "Pae3681" => 3,
    "Pae2250" => 3,
    "PgyB076" => 3,
    "PgyR4Q" => 3,
    "Pja301072" => 2,
    "Pla302278" => 1,
    "Pmo301020" => 3,
    "Pmp302280" => 1,
    "Por16" => 4,
    "Ppa2367" => 2,
    "Pph1448A" => 3,
    "Pph1644" => 3,
    "PphNPS3121" => 3,
    "PphHB10Y" => 3,
    "Ppi1704B" => 2,
    "PseHC1" => 3,
    "Psv3335" => 3,
    "Psy1" => 2,
    "Psy108" => 1,
    "Psy11" => 5,
    "Psy117" => 2,
    "Psy119" => 2,
    "Psy13" => 3,
    "Psy135" => 1,
    "Psy14" => 3,
    "Psy15" => 3,
    "Psy16" => 3,
    "Psy27" => 3,
    "Psy30" => 3,
    "Psy34" => 5,
    "Psy37" => 1,
    "Psy38" => 1,
    "Psy39" => 1,
    "Psy40" => 5,
    "Psy41" => 3,
    "Psy44" => 3,
    "Psy47" => 3,
    "Psy49" => 3,
    "Psy52" => 4,
    "Psy55" => 3,
    "Psy56" => 3,
    "Psy66" => 3,
    "Psy67" => 3,
    "Psy80" => 2,
    "Psy83" => 3,
    "Psy94" => 2,
    "Psy96" => 1,
    "Psy97" => 3,
    "Psy98" => 2,
    "Psy102" => 3,
    "PsyB728A" => 2,
    "PsyFF5" => 2,
    "PsyCit7" => 2,
    "Pth2598" => 1,
    "Pta11528" => 3,
    "PtoDC3000" => 1,
    "PtoT1" => 1,
    "PtoK40" => 1,
    "PtoMax13" => 1,
    "Pto1108" => 1,
    "Ptt50252" => 2
  );
  my $phylo_comp;
  my $strains = shift;
  my @fields = split(/,/, $strains);
  unless ( $phylogroups{$fields[1]} ) { warn $fields[1], " not in phylgroups\n"; return 0; }
  unless ( $phylogroups{$fields[0]} ) { warn $fields[0], " not in phylgroups\n"; return 0; }
  if ($phylogroups{$fields[1]} == $phylogroups{$fields[0]} ) {
    $phylo_comp = "Intragroup";
  }
  elsif ($phylogroups{$fields[1]} > $phylogroups{$fields[0]} ) {
    $phylo_comp = 'G' . $phylogroups{$fields[0]} . '-G' . $phylogroups{$fields[1]};
  }
  else {
    $phylo_comp = 'G' . $phylogroups{$fields[1]} . '-G' . $phylogroups{$fields[0]};
  }
  return($phylo_comp);
}

#Compare sequenced names and determine if comparison is within a subfamily or between
#This info is written to the infile for R, but is not currently plotted. This can be easily
#added to the R script
sub GetSubFamily {
  my $fam1 = shift;
  my $fam2 = shift;
  my $fam_comp;
  my $subfam1;
  my $subfam2;
  if ( $fam1 =~ /hop[A-Z][A-Z]?(\d)/ ) { $subfam1 = $1; }
  elsif ( $fam1 =~ /avr[a-zA-Z]{1,3}(\d)/ ) { $subfam1 = $1; }
  else { die "$fam1 not parsed correctly\n"; }
  if ( $fam2 =~ /hop[A-Z][A-Z]?(\d)/ ) { $subfam2 = $1; }
  elsif ( $fam2 =~ /avr[a-zA-Z]{1,3}(\d)/ ) { $subfam2 = $1; }
  else { die "$fam2 not parsed correctly\n"; }
  if ( $subfam1 == $subfam2 ) { $fam_comp = "IntraSubFam"; }
  elsif ($subfam2 > $subfam1 ) {
    $fam_comp = 'SubFam' . $subfam1 . '-SubFam' . $subfam2;
  }
  else {
    $fam_comp = 'SubFam' . $subfam2 . '-SubFam' . $subfam1;
  }
  return($fam_comp);
}

  
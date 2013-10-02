#!/usr/bin/perl -w
=head1 NAME

ParseT3SEdb.pl
1 October 2013

=head1 SYNOPSIS

ParseT3SEdb.pl -d <database> [ -o <out_file> | -t <t3se_seq_dir> ] [ -a ]


=head1 DESCRIPTION

Parses tab delimited version of T3SE database from http://pseudomonas-syringae.org/T3SS-Hops-noformat.xls
and writes fasta-formated sequences to file(s)

***Linebreaks must be converted to Unix/LF)***

Writes DNA sequences by default, but can be used for amino acid sequences if -a flag is used

Sequences can be written to a single file or to individual files for each T3SE family

In the latter case, sequences are written to files in the specified directory (or in the 
directory containing the database file if none is specified)

Sequences are APPENDED to any existing files
=head2 NOTES

The formatting of the database is a bit funky, so this may miss some sequences that are
not in the correct field.

This can be used to obtain an updated version of the query set for finding T3SEs, but
we need to exclude PlaMAFF302278_hopAE1 and add include the non-secreted homologs of
T3SEs in NonT3SE_homologs.fa

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use autodie qw(open close);


# Hard-coded hash of sequences that are messed up in the database and should not be included
# in the query set used to identify novel T3SE homologs
my %excluded_seqs = ( 'PlaMAFF302278_hopAE1' => 1 );

my $aa = 0;
my $t3se_seq_dir;
my $db_name;
my $outfilename;
GetOptions( 'amino_acid' => \$aa,
	        'database=s' => \$db_name,
	        't3se_seq_dir=s' => \$t3se_seq_dir,
	        'outfilename=s' => \$outfilename
	      )
or die "failure to communicate\n";

unless ( $db_name and -e $db_name ) {
  die "please select valid database file\n";
}

if ( $t3se_seq_dir and $outfilename ) {
  die "please select either a file for all sequences OR a folder to write individual files for each family\n";
}

#extract file base name for use in default filenames
(my $name, my $path, my $ext) = fileparse($db_name, qr/\.[^.]*/);

#Create default directory to save T3SE sequences if not provided
unless ( $t3se_seq_dir ) {
  $t3se_seq_dir = $path;
}
unless ( $t3se_seq_dir =~ /\/$/ ) {
  $t3se_seq_dir .= "/";
}

open(my $infile, "<", $db_name);

my $header;
my @t3ses;
while (<$infile>){
  chomp;
  unless ( $header ) {
    $header = $_;
    next;
  }
  $_ =~ s/[" ()]//g;
  $_ =~ /^(\d+)/;
  unless ( $1 and $1 == scalar(@t3ses) + 1 ) {
    $t3ses[-1] .= $_;
    next;
  }
  push(@t3ses, $_);
}

foreach(@t3ses) {
  my @fields = split("\t", $_);
  my $outfile;
  if ($outfilename) {
    $outfile = $outfilename;
  }
  else {
    my $fam = $fields[3];
    $fam =~ s/[\d'-]//g;
    $outfile = $t3se_seq_dir . $fam . ".fa";
  }
  my $header = $fields[4] . $fields[5] . "_" . $fields[3];
  if ( $aa ) {
    unless ($fields[16] =~ /^[ACGT]+$/i or $fields[16] =~ /\W/ or $fields[16] =~ /^$/ ) {
      if ( $excluded_seqs{$header} ) { next; }   #skip seqs in excluded set
      open(my $fh, ">>", $outfile);
      print $fh ">$header\n", uc($fields[16]), "\n";
      close($fh);
    }
  }
  else {
    if ($fields[15] =~ /^[ACGT]+$/i) {
      open(my $fh, ">>", $outfile);
      print $fh ">$header\n", uc($fields[15]), "\n";
      close($fh);
    }
  }
}
  

#!/usr/bin/perl -w

=head1 NAME

MLST.pl version 1, 22 Nov 2010

=head1 SYNOPSIS

MLST.pl filename/directory_name

=head1 DESCRIPTION

Takes a genome sequence or Blast db and blasts all MLST nucleotide seqs against it.
Hit strings are concatinated and printed to a file

Details:


=head2 NOTES

This is currently written to extract the hit string from the blast output. This means that
it often misses the very beginning and the end of the locus.

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-utoronto-dot-caE<gt>

=cut
####################################################################################################

use strict;
use warnings;
use Bio::SearchIO;
use Getopt::Long;
use File::Basename;
use autodie qw(open close);


my $mlst_file;
my $blastdatabase;
my $outfilename;
my $strain_name;

GetOptions( 'mlst_file=s' => \$mlst_file,
	        'infile=s' => \$blastdatabase,
	        'outfile=s' => \$outfilename,
	        'name=s' => \$strain_name
	      )
or die "failure to communicate\n";

#Check input files and provide default values if necessary
unless ( $mlst_file and -f $mlst_file ) {
  die "Please select valid file for reference MLST sequences\n";
}
if ( $blastdatabase and -f $blastdatabase . '.nsq' and -f $blastdatabase . '.nin' and -f $blastdatabase . '.nhr'  ) { #root of blast database supplied
  next; 
}
elsif ( $blastdatabase and -f $blastdatabase ) {         #sequence file provided. make blast database
  print `makeblastdb -in $blastdatabase -dbtype nucl`;
}
else {
  die "Please select valid genome sequence file or blast database\n";
}

(my $name, my $path, my $extension) = fileparse($blastdatabase, qr/\.[^.]*/);

unless ( $outfilename ) {
  $outfilename = $path . $name . "_mlst.fa";
}

unless ( $strain_name ) {
  $strain_name = $blastdatabase;
   $strain_name =~ s/.fa*//;
}  

#blast reference sequences agains blast database, then read in results and concatinate hit strings
print `blastn -query $mlst_file -db $blastdatabase -out temp.bl -num_threads 4 -evalue 1e-20`;

my $concat_string;
my $blastIO = new Bio::SearchIO(-format => 'blast',
                           -file   => 'temp.bl');
while ( my $blast_result = $blastIO->next_result ) {
  my $hitnum = 0;
  my %sig_hits;
  while ( my $hit = $blast_result->next_hit ) {
    while ( my $hsp = $hit->next_hsp ) {
      if ( $hsp->evalue < 1e-20 ) {
        $sig_hits{$hsp->start('query')} = $hsp;   #make hash of hsps, with start positions as keys
        $hitnum ++;
      }
    }
  }
  unless ( $hitnum ) { print "$_ ", $_->display_id(), ": no hits!\n";}
  foreach my $key (sort { $a <=> $b } (keys(%sig_hits))) {      #sort hsps by start position
    my $hit_string = $sig_hits{$key}->hit_string;
    $concat_string .= $hit_string;
  }
}

#Remove blast result file
print `rm temp.bl`;

#Save sequence string to file
open(my $out_fh, ">$outfilename");
print $out_fh ">$strain_name\n";
print $out_fh "$concat_string\n";
close($out_fh);
exit;

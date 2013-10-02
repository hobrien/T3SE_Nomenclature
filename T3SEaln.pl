#!/usr/bin/perl -w

=head1 NAME

T3SEaln.pl
9 Nov 2012

=head1 SYNOPSIS

T3SEaln.pl -i <sequence_file> [ -o <output_directory> ]


=head1 DESCRIPTION

Runs Phylogentic pipeline on specified sequence file:

-runs translatorX to align Amino Acid translations of sequences and back-translate to DNA
-adds any sequences with frameshifts to DNA alignment (must be in file called <sequence_file>_FS.fa)
-Converts alignment to phylip (and MEGA) formats
-runs PhyML on phylip alignment and copies tree and stats output files to separate directories

Files are stored in specified directory, or in same directory as input file if none specified

=head2 NOTES

AA translations stored in <outputdir>/AA/Fasta
Aligned AA sequences stored in <outputdir>/AA/Alignment/Fasta
Aligned DNA sequences stored in <outputdir>/DNA/Alignment/Fasta
Phylip formated alignment stored in <outputdir>/DNA/Alignment/Phylip
MEGA formated alignment stored in <outputdir>/DNA/Alignment/Mega
Tree file (Newick format) stored in <outputdir>/DNA/Tree
Tree stats stored in <outputdir>/DNA/Tree_stats

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################


use warnings;
use strict;

my $outputdir;
my $infilename;
use Getopt::Long;
use File::Basename;

GetOptions( 'infilename=s' => \$infilename,
	        'outfilename=s' => \$outputdir
	      )
or die "failure to communicate\n";
 
unless ( $infilename and -e $infilename ) {
  die "please select valid sequence file\n";
}

#extract file base name for use in default filenames
(my $name, my $path, my $ext) = fileparse($infilename, qr/\.[^.]*/);
my $fs = $name . "_FS";
my $aa = $name . "_AA";
my $stats = $name . "_stats.txt";

#Create default output directory name if not provided
unless ( $outputdir ) { $outputdir = $path; }
$outputdir =~ s/\/$//;

#Create folders for output files (if not present)
unless ( -e $outputdir . "/DNA" ) { print `mkdir $outputdir/DNA`; }
unless ( -e $outputdir . "/DNA/Alignment" ) { print `mkdir $outputdir/DNA/Alignment`; }
unless ( -e $outputdir . "/DNA/Alignment/Fasta" ) { print `mkdir $outputdir/DNA/Alignment/Fasta`; }
unless ( -e $outputdir . "/DNA/Alignment/Phylip" ) { print `mkdir $outputdir/DNA/Alignment/Phylip`; }
unless ( -e $outputdir . "/DNA/Alignment/Mega" ) { print `mkdir $outputdir/DNA/Alignment/Mega`; }
unless ( -e $outputdir . "/DNA/Tree" ) { print `mkdir $outputdir/DNA/Tree`; }
unless ( -e $outputdir . "/DNA/Tree_stats" ) { print `mkdir $outputdir/DNA/Tree_stats`; }
unless ( -e $outputdir . "/AA" ) { print `mkdir $outputdir/AA`; }
unless ( -e $outputdir . "/AA/Fasta" ) { print `mkdir $outputdir/AA/Fasta`; }
unless ( -e $outputdir . "/AA/Alignment" ) { print `mkdir $outputdir/AA/Alignment`; }
unless ( -e $outputdir . "/AA/Alignment/Fasta" ) { print `mkdir $outputdir/AA/Alignment/Fasta`; }

#Run translatorX:
print "translatorx_vLocal.pl -i $infilename -p F\n";
print `perl ./translatorx_vLocal.pl -i $infilename -p F`;
print "\n\n########################################################################################\n\n";
print `mv translatorx_res.aaseqs.fasta $outputdir/AA/Fasta/$aa.fa`;
print `mv translatorx_res.aa_ali.fasta $outputdir/AA/Alignment/Fasta/$aa.fa`;
#Add FS sequences:
if ( -e "$outputdir/DNA/Fasta/$fs.fa" ) {
  print "mafft --add $outputdir/DNA/Fasta/$fs.fa > $outputdir/DNA/Alignment/Fasta/$name.fa\n\n";
  print `mafft --add $outputdir/DNA/Fasta/$fs.fa translatorx_res.nt_ali.fasta > $outputdir/DNA/Alignment/Fasta/$name.fa`;
}
  else {
  print `mv translatorx_res.nt_ali.fasta $outputdir/DNA/Alignment/Fasta/$name.fa`;
}
print `rm translatorx_res* translatorx_logfile.txt`;
print "\n\n########################################################################################\n\n";

#Convert to phylip and to Mega:
print `perl ./ConvertSeq.pl -i $outputdir/DNA/Alignment/Fasta/$name.fa -f phyext -o $outputdir/DNA/Alignment/Phylip/$name.phy -r`;
print `perl ./ConvertSeq.pl -i $outputdir/DNA/Alignment/Fasta/$name.fa -f mega -o $outputdir/DNA/Alignment/Mega/$name.meg -r`;

#Make Tree
print "Making Tree...\n";
print `phyml -i $outputdir/DNA/Alignment/Phylip/$name.phy`;
print `mv $outputdir/DNA/Alignment/Phylip/$name.phy_phyml_tree.txt $outputdir/DNA/Tree/$name.nwk`;
print `mv $outputdir/DNA/Alignment/Phylip/$name.phy_phyml_stats.txt $outputdir/DNA/Tree_stats/$stats`;

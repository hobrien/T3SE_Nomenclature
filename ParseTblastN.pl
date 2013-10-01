#!/usr/bin/perl -w

=head1 NAME

ParseTblastN.pl
9 Nov 2012

=head1 SYNOPSIS

ParseTblastN.pl -i <blast_table> -s <sequence_file> [ -n <strain_name> -o <out_file> -t <t3se_seq_dir> ]


=head1 DESCRIPTION

Parses tblastn output to identify ORFs and pseudogenes homologous to query sequences

Sequences are written to fasta files named after the gene family of the query in the 
specified directory (or in the directory containing the original genome sequences)

Sequence files are APPENDED, so sequences from multiple genomes can be written to the
same files

A summary of the coordinates of each homolog and it's psuedogene status is written to the
specified outfile (or to "<seqfilename>_out.txt" )

Strain name is added to each sequence id in outfiles. If it is not provided, it is taken from the sequence file name

=head2 NOTES

Blast commands: 
makeblastdb -in <sequence_file> -dbtype nucl
tblastn -query T3SE_allAA.fa -db <sequence_file> -evalue 1e-20 -outfmt 6 -num_threads 6 -out <blast_table>

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################





use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Index::Fasta;
use List::Util qw(min max);
use autodie qw(open close);
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";


#Several parameters are going to affect how this works. I am going to hard code defaults for now
#but at some point I am going to allow modification of theses with comman line flags. 

#Proportion of overlap (relatve to shorter sequence) that will cause two results from DIFFERENT
#effector families to be merged, with the name given to the family with the higher score (if queries
#are in the same family thay are always merged as long as the spacing of the query matches)
my $min_hit_overlap = .9;

#size (bp) of gap between non-overlapping hits from queries for the SAME family that will be included
#in the effector sequence.
my $max_hit_gap = 100;

#size of gap/overlap (in bp) between the in the query and the in the hit. If this is large, it
#means that a gene is tandemly duplicated. If it is large and negative, in means that a lot of
#non-homologous sequence is inserted in the target.
my $max_query_overlap = .1;
my $max_query_gap = 100;

#Prefer ATG start codons over alternative start codons if they are within this range 
#(in base pairs) of the optimal start codon
my $atg_pref = 15;

# Maximum distance (in bp) after the start of the blast hit that start codons will be
# considered valid (start codon positions that would truncate the homologous N-terminus
# by more than this amount are ignored.
my $upstream_start = 21;

#Collect filenames and check if valid
my $infilename;
my $outfilename;
my $seqfilename;
my $t3se_seq_dir;
my $strain_name;
GetOptions( 'seqfilename=s' => \$seqfilename,
	        'infilename=s' => \$infilename,
	        'outfilename=s' => \$outfilename,
	        't3se_seq_dir=s' => \$t3se_seq_dir,
	        'name=s' => \$t3se_seq_dir
	      )
or die "failure to communicate\n";
 
unless ( $infilename and -e $infilename ) {
  die "please select valid blast table file\n";
}
unless ( $seqfilename and -e $seqfilename ) {
  die "please select valid sequence file\n";
}

#extract file base name for use in default filenames
(my $name, my $path, my $ext) = fileparse($seqfilename, qr/\.[^.]*/);

#Create default outfile name if not provided
unless ( $outfilename ) {
  $outfilename = $path . $name . "_out.txt";
}

#Create default directory to save T3SE sequences if not provided
unless ( $t3se_seq_dir ) {
  $t3se_seq_dir = $path;
}
unless ( $t3se_seq_dir =~ /\/$/ ) {
  $t3se_seq_dir .= "/";
}

#Create default strain name
unless ( $strain_name ) {
  $strain_name = $name;
}

#Create sequence index file
my $inx_name = $path . $name . ".inx";
my $inx = Bio::Index::Fasta->new(-filename => $inx_name, -write_flag => 1);
$inx->make_index($seqfilename);

#Parse blast_file and create hash of non-overlapping blast hits
my $blastIO = new Bio::SearchIO(-format => 'blasttable',
                           -file   => $infilename);
print STDERR "Parsing blast hits...\n";
my %t3se_homo = %{BinHits($blastIO, $inx)};

#Open outfile and write tabular info about blast hits
open(OUT, ">$outfilename");

foreach my $t3se ( keys %t3se_homo ) {
  my $hsp = $t3se_homo{$t3se};
  
  #Fetch the contig containing the blast hit and reverse-complement if necessary
  my $seq_name = $hsp->hit->location->seq_id;
  my $seq = $inx->fetch($seq_name) or die "could not find sequence $seq_name\n";
  unless ( $hsp->hit->strand == 1 ) {
    $seq = $seq->revcom;
    my $start = $seq->length - $hsp->hit->end + 1;
    my $end = $seq->length - $hsp->hit->start + 1;
    $hsp->hit->start($start);
    $hsp->hit->end($end);
  }
  my $id = $strain_name . "_" . $t3se;
  
  #Determine the optimal open reading frame
  my $start = FindStart($seq->seq, $hsp, $id);
  if ( $hsp->hit->start - $start > 30 ) { print STDERR "$id: first ", $hsp->hit->start - $start, " bp are non-homologous. Possible chimeric N-terminus\n"; }
  $hsp->hit->start($start);
  my $stop = FindEnd($seq->seq, $hsp, $id);
  if ( $stop - $hsp->hit->end  > 30 ) { print STDERR "$id: last ", $stop - $hsp->hit->end, " bp are non-homologous. Possible chimeric C-terminus\n"; }
  $hsp->hit->end($stop);
  
  #Extract ORF subsequence and reduce strings of Ns to 5 or less, check for frameshifts and premature stops and write to file
  my $outstring = $seq->subseq($start, $stop);
  $outstring =~ s/(NNN)+/NNN/g;
  $t3se =~ s/-\d+//;                 #Remove -1, -2 etc suffix from family name
  my $pseudo = CheckORF($outstring);
  if ( $pseudo ) {                   #Format sequence names to reflect pseudogene status
    if ( $pseudo =~ /N-terminal/ ) { $id .= "_N"; }
    elsif ( $pseudo =~ /C-terminal/ ) { $id .= "_C"; }
    elsif ( $pseudo =~ /internal fragment/ ) { $id .= "_internal_frag"; }
    if ( $pseudo =~ /Frameshift/ ) { $outfilename =~ s/\.fa/_FS.fa/; $id .= "_FS"; }
    elsif ( $pseudo =~ /stop codons/ ) { $id .= "'"; }
  }

  my $outseq_filename = $t3se_seq_dir . $t3se . '.fa';
  my $outseq_file = Bio::SeqIO->new('-file' => ">>$outseq_filename",
                             '-format' => 'fasta') or die "could not open seq file $outseq_filename\n";
  my $outseq = Bio::Seq->new(-seq => $outstring,
                             -alphabet => 'dna',
                             -id  => $id);
  $outseq_file->write_seq($outseq);
  if ( $hsp->hit->strand == 1 ) { print OUT join("\t", ($id, $seq->id, $start, $stop, 'forward', $pseudo)), "\n"; }
  else { print OUT join("\t", ($id, $seq->id, $seq->length - $stop + 1, $seq->length - $start + 1, 'reverse', $pseudo)), "\n"; }
}
close(OUT);
exit;

#Returns a hash of non-overlapping blast hits (with overlapping / adjacent hits combined 
#according to the "min_hit_overlap", "max_hit_gap" and "max_query_overlap" parameters
sub BinHits {
  my $blastIO = shift;
  my $inx = shift;
  my %hits = %{SortHits($blastIO)};
  my %t3se_homo;
  foreach ( keys %hits ) {
    my @contig = @{$hits{$_}};
    while (my $hsp = pop(@contig) ) {   #work from END of sorted array of hits
      my $print = 1;
      my $x = scalar(@contig) - 1;
      while ($x >= 0) {
        if ( HitOverlap($contig[$x], $hsp, $inx) ) {  #if hsps overlap
          if ( $hsp->bits < $contig[$x]->bits ) {   #if new hsp has higher score than current one (discard current result)
            $contig[$x] = AdjustRange($contig[$x], $hsp);
            $print = 0;    #do not print result
            $x = -1;  #skip to end of hsp array
          } 
          else {                   #if current hsp has higher score than new one (discard new one from array)
            $hsp = AdjustRange($hsp, $contig[$x]);
            splice(@contig, $x, 1);
            $x --;
          }
        } 
        else { $x --; }
      } 
      if ( $print ) {
        if ( $hsp->location->seq_id =~ /(mltB)|(LysM)|(PSPTO)|(hrp)/i ) { next; } #Query seq list includes some non-T3SEs with homology to T3SEs. Genes with these as the top hit should be skipped
        my $name = GetFamily($hsp);
        if ( $t3se_homo{$name} ) {  #Adding a second family member (need to use -1, -2, etc notation)
          $t3se_homo{$name . '-1'} = $t3se_homo{$name}; #Change the name of the first member from hopXX to hopXX-1
          delete $t3se_homo{$name};
        } 
        if ( $t3se_homo{$name . '-1'} ) { #Cycle through numbers until unused one found
          my $x = 2;
          while ( $t3se_homo{$name . '-' . $x} ) { $x ++; }
          $t3se_homo{$name . '-' . $x} = $hsp;
        } 
        else { $t3se_homo{$name} = $hsp; } #First family member. -1 not needed
      } 
    }
  } 
  return \%t3se_homo;
}

# returns true if hits are on the same strand and if they overlap by more than the 
# "min_hit_overlap" threshold or are from the same query family and are within the 
# "max_hit_gap" threshold of each other (not counting Ns introduced during scaffolding)
# unless queries overlap by more than "max_query_overlap" more than the hits overlap 
# (indicating some sort of tandem duplication). If the gap between queries is more than
# "max_query_gap" larger than the gap between hits, a warning about a internal deletion is printed
sub HitOverlap {
  my $hsp1 = shift;
  my $hsp2 = shift;
  my $inx = shift;
  my $combine = 0;
  if ( $hsp1->hit->start > $hsp2->hit->start ) { die "second hit starts before the first. hits are not sorted correctly\n"; }
  unless ( $hsp1->hit->strand == $hsp2->hit->strand ) { return $combine; } #different strands. do not combine
  my $hit_overlap = $hsp1->hit->end - $hsp2->hit->start + 1;
  if ( $hit_overlap / min($hsp1->hit->length, $hsp2->hit->length) > $min_hit_overlap or $hsp1->hit->end > $hsp2->hit->end ) { $combine = 1; } #overlap greater than "min_hit_overlap"
  elsif (GetFamily($hsp1) eq GetFamily($hsp2) ) {  #queries are from the same T3SE family. Combine if gap is less than "max_hit_gap"
    if ( $hit_overlap >= 0 ) { $combine = 1; }
    else {                                         #remove scaffold gap Ns from gap between hits before determining size
      my $id = $hsp1->hit->location->seq_id;
      my $seq = $inx->fetch($id) or die "could not find sequence $id\n";
      my $seq_string = $seq->subseq($hsp1->hit->end, $hsp2->hit->start);
      $seq_string =~ s/N+//g;
      $hit_overlap = 1 - length($seq_string);
      if ( $hit_overlap > 1 - $max_hit_gap ) { $combine = 1; }
    }
    my $query_overlap;
    if ( $hsp1->query->start < $hsp2->query->start ) {
      $query_overlap = ($hsp1->end - $hsp2->start) * 3 + 1;
    }
    else {
      $query_overlap = ($hsp2->end - $hsp1->start) * 3 + 1;
    }
    if ( ($query_overlap - $hit_overlap) / max($hsp1->hit-length, $hsp2->hit->length) > $max_query_overlap ) {
      print STDERR $hsp1->location->seq_id, ": query sequences overlap by $query_overlap bp. Hit sequences overlap by $hit_overlap bp. This may be a case of tandem duplication.\nHits will not be combined\n";
      $combine = 0;
    }
    elsif ( $query_overlap - $hit_overlap < 0 - $max_query_gap ) {
      print STDERR $hsp1->location->seq_id, ": gap between query sequences is ", 0 - $query_overlap + $hit_overlap, " bp more than hit sequence gap.  This likely indicates an internal deletion\n";
    }
  }
  return $combine;
}

# this groups hsps by contig and sorts them by start position
# it also adds query_name to hsp ( as $hsp->location->seq_id )
sub SortHits {
  my $blastIO = shift;
  my %hits;
  while ( my $result = $blastIO->next_result ) {
    while ( my $hit = $result->next_hit ) {
      while ( my $hsp = $hit->next_hsp ) {
        $hsp->hit->location->seq_id($hit->name);
        $hsp->location->seq_id($result->query_name);
        if ( $hits{$hit->name} ) {
          my @contig = @{$hits{$hit->name}};
          my $x = 0;
          while ( $x < @contig ) {
            if ( $hsp->hit->start < $contig[$x]->hit->start ) { last; }
            else { $x ++; }
          }
          splice(@contig, $x, 0, $hsp);
          $hits{$hit->name} = \@contig;
        }
        else { $hits{$hit->name} = [ $hsp ]; }
      }
    }
  }
  return \%hits;
}

# Combines two blast hits into a single hit with the minimum start and maximum end of the individual hits
sub AdjustRange {
  my $hsp1 = shift;
  my $hsp2 = shift;
  (my $start, my $stop, my $strand) = $hsp1->hit->location->union($hsp2->hit->location);
  $hsp1->hit->location->start($start);
  $hsp1->hit->location->end($stop);
  ($start, $stop, $strand) = $hsp1->query->location->union($hsp2->query->location);
  $hsp1->query->location->start($start);
  $hsp1->query->location->end($stop);
  return $hsp1;
}

# Parse the hit name to determine which T3SE family it belongs to
sub GetFamily {
  my $hsp = shift;
  my $location = $hsp->location->seq_id;
  unless ( $location =~ /((hop[A-Z]{1,3})|(avr[A-Z]{1,3})|(mltB)|(LysM)|(PSPTO)|(hrp[A-Z]))/i) { 
    die "$location not parsed correctly\n" ;
  }
  return $1;
}

# Select start codon closest to corresponding start codon in query sequence (ATG codons are
# preferred if they are within "atg_pref" codons of an alternative start codon.
# Search starts "upstream_start" bp from the end of the end of the blast search
sub FindStart {
  my $contig = shift;
  my $hsp = shift;
  my $id = shift;
  my $start;
  my $x = $hsp->hit->start;
  my $atg;
  my $target_start = $x - ($hsp->query->start * 3) + 3;
  $x += $upstream_start;                            #start looking for a start codon "upstream_start" pb before the start of the blast hit
  while ($x > 1 ) {
    my $codon = substr($contig, $x-1, 3);
    if ( $codon =~ /(TAA)|(TAG)|(TGA)|(NNN)/ ) {  #stop codon
      unless ( $start ) {
        $start = $hsp->hit->start;
        print STDERR "$id: in-frame stop codon before start codon ( C-terminal truncation )\n";
      }
      last;
    }
    if ( $codon =~ /ATG/ ) {                     #preferred start
      unless ( $start ) { $start = $x; }
      if ( $atg ) {
        unless (abs($x - $target_start) > abs($start - $target_start ) ) { $start = $x; }
      }
      else {
        unless (abs($start - $target_start) - abs($x - $target_start ) > $atg_pref  ) { $start = $x; }
      }
      $atg = 1;
    }
    elsif ( $codon =~ /(AT[ACGT])|([CGT]TG)/ ) { #alternative start
      unless ( $start ) { $start = $x; }
      if ( $atg ) {
        if (abs($start - $target_start) - abs($x - $target_start ) > $atg_pref  ) { $start = $x; }
      }
      else {
        unless (abs($x - $target_start) > abs($start - $target_start ) ) { $start = $x; }
      }
      $atg = 0;
    }
    $x -=3;
  }
  unless ( $start ) {
    print STDERR "$id: no start codon, setting to start of blast hit\n";
    $start = $hsp->hit->start;
  }
  return $start;
}

# This is a straightforward function to find the first in-frame stop codon (if one is present)
sub FindEnd {
  my $contig = shift;
  my $hsp = shift;
  my $id = shift;
  my $stop;
  my $x = $hsp->hit->end;
  $x -=21;
  while ($x < length($contig) ) {
    if ( substr($contig, $x - 3, 3) =~ /(TAA)|(TAG)|(TGA)/ ) { $stop = $x; last; }
    if ( substr($contig, $x - 3, 3) =~ /NNN/ ) { last; }
    $x +=3;
  }
  unless ( $stop ) {
    print STDERR "$id: no stop codon, setting to end of blast hit\n";
    $stop = $hsp->hit->end;
  }
  return $stop;
}

#This will test if the sequence has start and stop condons, if it is in frame and if there are internal stops
sub CheckORF {
  my $seq = shift;
  my $pseudo = '';
  if ($seq =~ /[^ACGTN]/i ) { die "non-canonical bases in seq\n"; }
  unless ( $seq =~ s/(TGA)$|(TAA)$|(TAG)$// ) { $pseudo = 'N-terminal fragment'; } #remove stop codon
  unless ( $seq =~ /^(AT[ACGT])|^([CGT]TG)/ ) {                                    #check for start codon
    if ( $pseudo ) { $pseudo = 'internal fragment'; }
    else { $pseudo = 'C-terminal fragment'; }
  }
  if ( length($seq) % 3 and $seq !~ /NNN/ ) { $pseudo .= ", Frameshift"; }         #Check for frameshifts in sequences without Ns
  #Check for in-frame stop codons in first section before contig break
  $seq =~ s/([ACGT]+)N*//i;
  if ( $1 =~ /^([ACGT]{3})+((TGA)|(TAA)|(TAG))/ ) { 
    $pseudo .= ", Internal stop codon";
  } 
  #Check each internal contig of scaffolds (sequence between runs of Ns) for stop codons
  my @segments;
  while ( $seq =~ s/([ACGT]+)N+//i ) { push(@segments, $1); } 
  push(@segments, $seq);
  foreach my $subseq ( @segments ) {
    my $length;
    my $frame1stop;
    my $frame2stop;
    my $frame3stop;
    $subseq =~ s/(TGA)|(TAA)|(TAG)/XXX/g; #Convert stop codons to XXX (makes it easier to search for them)
    while  ( $subseq =~ s/([^X]+XXX)// ) { #For each stop codon in any frame, 
      $length += length($1);
      if ( $length % 3 == 0 ) { $frame1stop ++; }
      elsif ( ($length + 2 ) % 3 == 0 ) { $frame2stop ++; }
      else { $frame3stop ++ ; }
    }
    if ( $frame1stop and $frame2stop and $frame3stop ) { $pseudo .= ", Internal stop codon"; }
  }
  #Check last section for stop codon in frame with the END of the sequence
  if ( $seq and $seq =~ /((TGA)|(TAA)|(TAG))([ACGT]{3}+)$/ ) { $pseudo .= ", Internal stop codon"; }
  if ( $pseudo ) { $pseudo =~ s/^, //; } #remove leading comma from start of list of pseudogene notes
  return $pseudo;
}




#!/usr/bin/perl -w
use warnings;
use strict;

my $name = shift or die "T3SE name not provided\n";
my $fs = $name . "_FS";
my $aa = $name . "_AA";

my $folder = '/Users/HeathOBrien/BTSync/Pseudomonas/Publications/Nomenclature/';
#Run translatorX:
print "translatorx_vLocal.pl -i $folder/T3SE/DNA/Fasta/$name.fa -p F\n";
print `translatorx_vLocal.pl -i $folder/T3SE/DNA/Fasta/$name.fa -p F`;
print "\n\n########################################################################################\n\n";
print "mv translatorx_res.aaseqs.fasta $folder/T3SE/AA/Fasta/$aa.fa\n";
print `mv translatorx_res.aaseqs.fasta $folder/T3SE/AA/Fasta/$aa.fa`;
print "mv translatorx_res.aa_ali.fasta $folder/T3SE/AA/Alignment/Fasta/$aa.fa\n";
print `mv translatorx_res.aa_ali.fasta $folder/T3SE/AA/Alignment/Fasta/$aa.fa`;
print "ConvertSeq.pl -i $folder/T3SE/AA/Alignment/Fasta/$aa.fa -f phyext -o $folder/T3SE/AA/Alignment/Phylip/$aa.phy -r\n";
print `ConvertSeq.pl -i $folder/T3SE/AA/Alignment/Fasta/$aa.fa -f phyext -o $folder/T3SE/AA/Alignment/Phylip/$aa.phy -r`;
print `ConvertSeq.pl -i $folder/T3SE/AA/Alignment/Fasta/$aa.fa -f aln -o $folder/T3SE/AA/Alignment/ClustalW/$aa.aln -r`;

#Add FS sequences:
if ( -e "$folder/T3SE/DNA/Fasta/$fs.fa" ) {
  print "mafft --add $folder/T3SE/DNA/Fasta/$fs.fa > $folder/T3SE/DNA/Alignment/Fasta/$name.fa\n\n";
  print `mafft --add $folder/T3SE/DNA/Fasta/$fs.fa translatorx_res.nt_ali.fasta > $folder/T3SE/DNA/Alignment/Fasta/$name.fa`;
}
  else {
  print `mv translatorx_res.nt_ali.fasta $folder/T3SE/DNA/Alignment/Fasta/$name.fa`;
}
print "\n\n########################################################################################\n\n";

#Convert to phylip:
print `ConvertSeq.pl -i $folder/T3SE/DNA/Alignment/Fasta/$name.fa -f phyext -o $folder/T3SE/DNA/Alignment/Phylip/$name.phy -r`;
print "ConvertSeq.pl -i $folder/T3SE/DNA/Alignment/Fasta/$name.fa -f aln -o $folder/T3SE/DNA/Alignment/ClustalW/$name.aln -r\n";
print `ConvertSeq.pl -i /$folder/T3SE/DNA/Alignment/Fasta/$name.fa -f aln -o $folder/T3SE/DNA/Alignment/ClustalW/$name.aln -r`;
print "ConvertSeq.pl -i$folder/T3SE/DNA/Alignment/Fasta/$name.fa -f mega -o $folder/T3SE/DNA/Alignment/Mega/$name.meg -r\n";
print `ConvertSeq.pl -i $folder/T3SE/DNA/Alignment/Fasta/$name.fa -f mega -o $folder/T3SE/DNA/Alignment/Mega/$name.meg -r`;

#Make Tree
print "Making Tree...\n";
my $stats = $name . "_stats.txt";
print `phyml -i $folder/T3SE/DNA/Alignment/Phylip/$name.phy`;
print `mv $folder/T3SE/DNA/Alignment/Phylip/$name.phy_phyml_tree.txt $folder/T3SE/DNA/Tree/$name.nwk`;
print `mv $folder/T3SE/DNA/Alignment/Phylip/$name.phy_phyml_stats.txt $folder/T3SE/DNA/Tree_stats/$stats`;

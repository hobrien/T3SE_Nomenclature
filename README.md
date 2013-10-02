T3SE_Nomenclature
=================

Scripts and datafiles necessary to run the Type 3 Secreted Effector Pipeline the I've developed.

Prerequisites:

 blast+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 
 BioPerl (http://www.bioperl.org/wiki/Main_Page)
 
 translatorx (http://pc16141.mncn.csic.es/cgi-bin/translatorx_vLocal.pl; slightly modified version included here)
 
 MAFFT (http://mafft.cbrc.jp/alignment/software/)
 
 PhyML (http://www.atgc-montpellier.fr/phyml/binaries.php)
 
 MEGA (http://www.megasoftware.net/)
 
 R (http://www.r-project.org/)
 
 ggplot2 (http://ggplot2.org/; R library)
 
Scripts:

 ParseTblastN.pl (takes blast output and genome sequences and extracts corresponding ORFs / psuedogenes)
 
 T3SEaln.pl (takes fasta file of new sequences from a gene family, adds it to fasta file of reference sequences and builds tree)
 
 translatorx_Local.pl (see Prerequisites above)

Datafiles:

 Psy108.fa (draft genome sequences from a T3SE rich strain for testing purposes)
 
 T3SE_allAA.fa (reference amino acid sequences to use as query for tBlastn)
 

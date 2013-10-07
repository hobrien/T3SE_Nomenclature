T3SE_Nomenclature: 
=================
*Scripts and datafiles necessary to run the Type 3 Secreted Effector Pipeline that I've developed.*
Prerequisites:
--------------

 * blast+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 * BioPerl (http://www.bioperl.org/wiki/Main_Page)
 * translatorx (http://pc16141.mncn.csic.es/cgi-bin/translatorx_vLocal.pl; slightly modified version included here)
 * MAFFT (http://mafft.cbrc.jp/alignment/software/)
 * PhyML (http://www.atgc-montpellier.fr/phyml/binaries.php)
 * MEGA (http://www.megasoftware.net/)
 * R (http://www.r-project.org/)
 * ggplot2 (http://ggplot2.org/; R library)
 
 
 
Scripts:
--------

* ConvertSeq.pl (script to convert among sequence and alignment file formats)
* MakeDistPlot.pl (combines distance data from T3SE alignment and from MLST alignment and makes a scatterplot)
* MLST.pl (blast MLST ref sequences against genome and extract hit sequences)
* ParseT3SEdb.pl (parses T3SS-Hops-noformat.txt and extracts DNA or AA sequences)
* ParseTblastN.pl (takes blast output and genome sequences and extracts corresponding ORFs / psuedogenes)
* T3SEaln.pl (takes fasta file of new sequences from a gene family, adds it to fasta file of reference sequences and builds tree)
* translatorx_Local.pl (see Prerequisites above)
 
 

Datafiles:
----------

* NonT3SE_homologs.fa (Homologs of T3SEs that need to be added to query set)
* Psy108.fa (draft genome sequences from a T3SE rich strain for testing purposes)
* T3SE_allAA.fa (reference amino acid sequences to use as query for tBlastn)
* T3SS-Hops-noformat.txt (T3SE db downloaded from http://pseudomonas-syringae.org/T3SS-Hops-noformat.xls, 
 converted to tab delimited with Unix/LF linebreaks)
* PsyMLST.fa (Alignment of MLST genes for SOME of the strains in the database (it would be more useful if more were added))
* PsyMLSTref.fa (reference MLST sequences to blast against genome when using MLST.pl)



Workflow:
---------
        perl ParseT3SEdb.pl -d T3SS-Hops-noformat.txt -o T3SEaa.fa -a

        cat NonT3SE_homologs.fa >> T3SEaa.fa

        perl ParseT3SEdb.pl -d T3SS-Hops-noformat.txt

        makeblastdb -in Psy108 -dbtype nucl

        tblastn -query T3SEaa.fa -db Psy108 -evalue 1e-20 -outfmt 6 -num_threads 6 -out <blast_table>

        perl ParseTblastN.pl -i Psy108_t3se.bl -s Psy108.fa

        perl MLST.pl -i Psy108.fa -m PSYMLSTref.fa

        mafft --add Psy108_mlst.fa PsyMLST.fa >PsyMLST_aln.fa

        perl ConvertSeq.pl -i PsyMLST_aln.fa -f mega 
 
For each T3SE family:

        T3SEaln.pl -i T3SE

        Calculate distances in Mega for effector and for MLST genes (see below)

        perl MakeDistPlot.pl -m PsyMLSTdist.txt -t DNA/Distances/T3SE.dist.txt



Running MEGA:
-------------

Currently, this step hasn't been automated. It is probably possible to calculate synonymous
distances in R, but I only know how to do total genetic distances. There's also something
called "MEGA Computational Core", which is scriptable, but it only works on Windows.

As I'm going through this, I'm also realising that there is no way to automatically
correct the reading frame when sequences with franeshifts cause the 3' end of the matrix to
be out of frame, so if we are going to stick with synonymous distances, we should either
add that in or skip the part of T3SEaln.pl that adds frameshifted sequences to the alignment

Current Steps:
* Open Mega file in MEGA. Select "File"->"Open A File/Session..."
* Select file then select "OK"
* Select "Nucleotide Sequences" and "OK" 
* Select "Yes" when asked if sequences are protein coding
* Select "Genetic Code" "Standard" or "Bacterial Plastid" and select "OK" (the only differences is the permitted start codons)
* Return to main window and select Distance->Compute Pairwise Distances. Select "Yes"
* Change "Substitution Type" to "Syn-Nonsynonymous"
* Change "Gaps/Missing Data Treatment" to "Pairwise Deletion"
* Change "Substitutions to Include" to "s: Synonymous only"
* Other parameters can be set to default:
    - "Variance Estimation Method": "None" 
    - "Model/Method": "Nei-Gojobori method (Jukes-Cantor)"
* Select "Compute".
* Select "File"->"Export/Print Distances" then change "Export Type" from "Matrix" to "Column".
* "Output Format" should be "CSV: Comma-separated file"
* Select "Print/Save Matrix" then save the resulting window




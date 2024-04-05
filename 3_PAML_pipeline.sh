##This code was run on a Cornell BioHPC server line-by-line in the terminal
##Here are my steps for PI;  other genes were the same but I changed the filenames

##Making the codon-aware alignment
export PATH=/programs/mafft/bin:$PATH
mafft --maxiterate 5000 --auto --adjustdirectionaccurately --thread 8 --op 3 --leavegappyregion PI_OG0004848_cds.fa > PI_OG0004848_cds.mafft.fa
java -jar /local/workdir/arj66/macse_software/macse_v2.07.jar -prog alignSequences -out_NT PI_OG0004848_cds.macse.mafft.fa -seq PI_OG0004848_cds.mafft.fa
##macse introduced exclamation points ie frameshift sequencing errors-- changed these to N's
sed -i 's/!/N/g' PI_OG0004848_cds.macse.mafft.fa
/programs/trimal-1.4/source/trimal -in PI_OG0004848_cds.macse.mafft.fa -out PI_OG0004848_cds.trimal.macse.mafft.fa -gappyout
java -jar /local/workdir/arj66/macse_software/macse_v2.07.jar -prog exportAlignment -align PI_OG0004848_cds.trimal.macse.mafft.fa -codonForFinalStop NNN -codonForInternalStop NNN -keep_gap_only_sites_ON -out_NT PI_OG0004848_cds.new.ready.trimal.macse.mafft.fa

##RAxML for trees
nohup /programs/raxml-ng_v1.2.0/raxml-ng --all --msa PI_OG0004848_cds.new.ready.trimal.macse.mafft.fa --msa-format FASTA --model GTR+G --prefix PI_OG0004848_cds.new.gene_tree --threads 12 --bs-trees 1000 &

##PAML to get dN/dS for each branch to plot
/programs/paml-4.10.6/bin/codeml 3_1_codeml_model1_NSSites0.ctl

##PAML to get which sites are under selection
/programs/paml-4.10.6/bin/codeml 3_2_codeml_model2_NSSites2.ctl

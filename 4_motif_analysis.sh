##This code was run on a Cornell BioHPC server line-by-line in the terminal

#got the LFY targets from Supplementary Data 1 from Jin et al 2021
#converting the bed to fasta of the target sequences
bedtools getfasta -fi TAIR10_chr_all.fas -bed LFY_ChIPseq.bed -fo LFY_ChIPseq.fa

#doing meme to identify the canonical LFY binding motifs in the Jin data
export PATH=/programs/meme-5.5.2/bin:/programs/meme-5.5.2/bin2:$PATH
nohup meme-chip -meme-p 6 -maxw 25 LFY_ChIPseq.fa &

#Got the mLUBS and dLUBS by email from Rieu et al 2023 lab group
/programs/meme_4.11.2_1/bin/jaspar2meme dLUBS_PFM_file -pfm > dLUBS.meme
export PATH=/programs/meme-5.5.2/bin:/programs/meme-5.5.2/bin2:$PATH
meme2images dLUBS.meme dLUBS_test_out
/programs/meme_4.11.2_1/bin/jaspar2meme mLUBS_PFM_file -pfm > mLUBS.meme
meme2images mLUBS.meme mLUBS_test_out

#Decided to investigate AP1 and AG which should be activated by LFY
#^(AP1 is questionable as there are multiple homologs and function may not be conserved)
#Also AP3 and PI which should be activated by LFY-UFO
#Also PHYA which should somewhat be a "control"

##looking for motifs 3kb upstream in Arabidopsis
export PATH=/programs/seqkit-0.15.0:$PATH
seqkit grep -f genes_of_interest.txt Araport11_upstream_3000_20220914 -o genes_of_interest_3kb_upstream.fa
#manually added the gene names

export PATH=/programs/meme-5.5.2/bin:/programs/meme-5.5.2/bin2:$PATH
nohup fimo RWCCAAYGGWCAAWR.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out Arabidopsis_canonical_LFY

nohup fimo mLUBS.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out Arabidopsis_mLUBS

nohup fimo dLUBS.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out Arabidopsis_dLUBS


##Looking for motifs 3kb upstream in Euphorbia peplus
grep "Ep_chr2_g04448\|Ep_chr6_g17595\|Ep_chr4_g09342\|Ep_chr8_g24863\|Ep_chr8_g24861\|Ep_chr8_g24862\|Ep_chr8_g25689\|Ep_chr8_g25690\|Ep_chr4_g08890" Euphorbia_peplus.gff > genes_of_interest.gff
grep "\bgene\b" genes_of_interest.gff > genes_of_interest_gene_lines.gff
#manually did the bed file
#selecting 3kb upstream based on the genes_of_interest_gene_lines.gff
##for ones on the plus strand, "upstream" is 3000 bases lower than the lowest #
##for the ones on the minus strand, "upstream" is 3000 bases higher than the highest #
bedtools getfasta -nameOnly -fi Euphorbia_peplus.fa -bed E_peplus_3kb_upstream.bed -fo genes_of_interest_3kb_upstream.fa

export PATH=/programs/meme-5.5.2/bin:/programs/meme-5.5.2/bin2:$PATH
nohup fimo RWCCAAYGGWCAAWR.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out E_peplus_canonical_LFY

nohup fimo mLUBS.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out E_peplus_mLUBS

nohup fimo dLUBS.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out E_peplus_dLUBS


##Looking for motifs 3kb upstream in Hevea brasiliensis (in Euphorbiaceae family)

grep "KAF2298949.1\|KAF2299252.1\|KAF2313700.1\|KAF2291016.1\|KAF2312883.1\|KAF2314627.1\|KAF2310046.1\|KAF2312359.1\|KAF2290821.1\|KAF2291335.1\|KAF2296683.1\|KAF2324094.1\|KAF2307360.1" genomic.gff > genes_of_interest.gff
#manually did the bed file
#selecting 3kb upstream based on the genes_of_interest_gene_lines.gff
##for ones on the plus strand, "upstream" is 3000 bases lower than the lowest #
##for the ones on the minus strand, "upstream" is 3000 bases higher than the highest #
bedtools getfasta -nameOnly -fi GCA_010458925.1_ASM1045892v1_genomic.fna -bed Hevea_3kb_upstream.bed -fo genes_of_interest_3kb_upstream.fa

export PATH=/programs/meme-5.5.2/bin:/programs/meme-5.5.2/bin2:$PATH
nohup fimo RWCCAAYGGWCAAWR.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out Hevea_canonical_LFY

nohup fimo mLUBS.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out Hevea_mLUBS

nohup fimo dLUBS.meme genes_of_interest_3kb_upstream.fa &
mv fimo_out Hevea_dLUBS

for file in seqs/pdb/*.pdb; do
	bname=${file##*/}
	bname=${bname%.fastq}
	grep CA seqs/pdb/$bname.pdb | cut -c18-20 > seqs/aa/$bname.seqaa 
done

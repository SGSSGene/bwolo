all:
	make -C fasta2Fmi
	make -C bwolo

clean:
	make -C fasta2Fmi clean
	make -C bwolo clean
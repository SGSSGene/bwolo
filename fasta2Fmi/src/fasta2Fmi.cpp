/*
 * Copyright 2014 Bonsai Bioinformatics Research Group
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <getopt.h>

#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include "BwoloTypeConfig.h"

using namespace std;
using namespace seqan;

void printHelp() {
	cout << "fasta2Fmi [option] -f [fasta] -i [out]" << endl;
	cout << "fasta2Fmi will index a sequence in a FM-index and save it in files"
			<< endl;
	cout
			<< "\t-f --fasta [required][required_argument:text] : Fasta file of the sequences."
			<< endl;
	cout
			<< "\t-i --index [required][required_argument:text] : Output file for the index."
			<< endl;
	cout << "\t-v --verbose [optional] : Verbose mode." << endl;
	cout << "\t-h --help print this message" << endl;
}

int main(int argc, char *argv[]) {
	char* indexFile = NULL;
	char* fastaFile = NULL;
	bool isVerbose = false;
	const char *shortOptions = "f:i:vh";
	static struct option longOptions[] = {
			{ "fasta", required_argument, NULL,	'f' },
			{ "index", required_argument, NULL, 'i' },
			{ "verbose", no_argument, NULL, 'v' },
			{ "help", no_argument, NULL, 'h' },
			{ NULL, 0, NULL, 0 } };
	int optionIndex = 0;
	int nextOption = -1;
	nextOption = getopt_long(argc, argv, shortOptions, longOptions,
			&optionIndex);
	while (nextOption != -1) {
		switch (nextOption) {
		case 'f':
			fastaFile = optarg;
			break;
		case 'i':
			indexFile = optarg;
			break;
		case 'v':
			isVerbose = true;
			break;
		case 'h':
			printHelp();
			exit(EXIT_SUCCESS);
			break;
		case '?':
			printHelp();
			exit(EXIT_FAILURE);
			break;
		}
		nextOption = getopt_long(argc, argv, shortOptions, longOptions,
				&optionIndex);
	}
	if (optind < argc) {
		if (isVerbose) cout << ("this arguments are not recognized and will be ignored : ");
		while (optind < argc)
			if (isVerbose) cout << '\t' << argv[optind++] << endl;
		printHelp();
	}
	if(indexFile == NULL || fastaFile == NULL){
		cout << "Fasta input file or Index output file are missing"<<endl;
		printHelp();
		exit(EXIT_FAILURE);
	}

	if(isVerbose)printDebugLevel(cout);

	if (isVerbose) cout << "opening fasta file : " << fastaFile << " : ";
	SequenceStream seqStream(fastaFile);
	if (!isGood(seqStream)) {
		if (isVerbose) cout << "Error : could not open " << fastaFile << endl;
		exit(EXIT_FAILURE);
	} else {
		if (isVerbose) cout << "OK" << endl;
	}

	if (isVerbose) cout << "reading the all sequences : ";

	CharString id;
	BwoloTText seq;
	BwoloTText tmp_seq;
	size_t sequences_ct = 0;
	while (readRecord(id, tmp_seq, seqStream) == 0) {
		seq += tmp_seq;
		sequences_ct += 1;
	}
	if (isVerbose) cout << "OK - " << sequences_ct << " sequences; total length: " << length(seq) << endl;

	if (isVerbose) cout << "Indexing sequence : " << id << " : ";
	BwoloTIndex myIndex(seq);
	//les index étant créés à la demande, il faut forcer la création des différentes "fibres"
	if(!indexCreate(myIndex)) {
		if (isVerbose) cout << "NOK" << endl;
		cout << "Error : could not indexing the sequence" << endl;
		close(seqStream);
		exit(EXIT_FAILURE);
	}
	if (isVerbose) cout << "OK" << endl;

	if (isVerbose) cout << "saving into file : "<<indexFile;
	if(!save(myIndex, indexFile)){
		if (isVerbose) cout << "Error : There was an error writing : " << indexFile << endl;
		close(seqStream);
		exit(EXIT_FAILURE);
	}else{
		if (isVerbose) cout << "OK" << endl;
	}
	close(seqStream);
	exit(EXIT_SUCCESS);
}

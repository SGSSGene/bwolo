/*
 * Copyright 2014 Bonsai Bioinformatics Research Group
 *
 * This file is part of Bwolo.
 *
 * Bwolo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bwolo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Bwolo.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <getopt.h>

#include "../include/interfaceMainFonction.h"
/*
 * test
 */
//#include "../tests/test_testAll.h"


using namespace std;
using namespace seqan;


void printHelp(){
	cout << "bwolo [option] -p [pattern] -i [index]" << endl;
	cout << "bwolo [option] -f [fasta] -i [index]" << endl;
	cout << "bwolo will find the seed from a pattern in the indexed sequence"
			<< endl;
	cout
	<< "\t-p --pattern [optional*][required_argument:text] : pattern to search."
	<< endl;
	cout << "\t-f --fasta [optional*][required_argument:text] : fasta file containing the patterns to search."
	<< endl;
	cout
	<< "\t-i --index [required][required_argument:text] : Output file for the index."
	<< endl;
	cout
	<< "\t-s --score [optional][requiered_argument:int] : minimal alignment score. Default : -3"
	<< endl;
	cout << "\t-r --rev-compl [optional] : add inverse complement" << endl;
	cout << "\t-v --verbose [optional] : Verbose mode." << endl;
	cout << "\t-h --help print this message" << endl;
	cout << "* : required one of them" << endl;
}


int main(int argc, char *argv[])
{
	/*    testAll();*/

	char* indexFile = NULL;
	char* patternCstr = NULL;
	char* fastaFile = NULL;
	int minScore=-3;
	bool revCompl = false;
	bool isVerbose = false;
	const char *shortOptions = "p:f:s:i:rvh";
	static struct option longOptions[] = {
			{ "pattern", required_argument, NULL,	'p' },
			{ "fasta", required_argument, NULL,	'f' },
			{ "score", required_argument, NULL,	's' },
			{ "index", required_argument, NULL, 'i' },
			{ "rev-compl", no_argument, NULL, 'r'},
			{ "verbose", no_argument, NULL, 'v' },
			{ "help", no_argument, NULL, 'h' },
			{ NULL, 0, NULL, 0 } };
	int optionIndex = 0;
	int nextOption = -1;
	nextOption = getopt_long(argc, argv, shortOptions, longOptions,
			&optionIndex);
	while (nextOption != -1) {
		switch (nextOption) {
		case 'p':
			patternCstr = optarg;
			break;
		case 'f':
			fastaFile = optarg;
			break;
		case 's':
			minScore = (int) strtol(optarg, NULL, 0);;
			break;
		case 'i':
			indexFile = optarg;
			break;
		case 'r':
			revCompl = true;
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
	if(indexFile == NULL || (patternCstr == NULL && fastaFile == NULL ) ){
		cout << "Index input file or pattern are missing"<<endl;
		printHelp();
		exit(EXIT_FAILURE);
	}

	unsigned short errNumber=(unsigned short)-minScore;
	if(fastaFile != NULL){
		printPotentialOccurences_PatternFasta(fastaFile, indexFile, errNumber, revCompl, isVerbose);
	}
	if(patternCstr != NULL){
		printPotentialOccurences_PatternSequences(patternCstr, indexFile, errNumber, revCompl, isVerbose);
	}
	exit(EXIT_SUCCESS);
}

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

#include "../include/interfaceMainFonction.h"

#include <iostream>
#include <set>

#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index/index_fm.h>
#include <seqan/file.h>

#include <seqan/index.h>
#include <seqan/seeds.h>
/////test/////////
#include <seqan/align.h>

#include "../include/ExtendRight.h"
///fin test///////

#include "../include/PatternSegment.h"
#include "../include/PatternFiltrer.h"
#include "../include/PatternSegment.h"

#include "../include/AlignmentInFmi.h"
using namespace std;
using namespace seqan;

template <typename TText> class FinalOccurence{
public:
	typename Position<TText>::Type beginPosition;
	typename Position<TText>::Type endPosition;
	unsigned short score;
	//constructeur
	FinalOccurence(typename Position<TText>::Type beginPosition, typename Position<TText>::Type endPosition, unsigned short score):
		beginPosition(beginPosition),
		endPosition(endPosition),
		score(score){}
};

/*
 * les occurences sont triées d'abord par leur position de début puis si egal, par leur position de fin.
 */
template <typename TText> class LessFinalOccurenceBeginPositionFirst
{
public:
	bool operator()(FinalOccurence<TText> const& x, FinalOccurence<TText> const& y){
		//XXX : les parenthèses sont necessaires pour compiler : x.beginPosition < y.beginPosition ne fonctionne pas.
		return ( (x.beginPosition) < (y.beginPosition) ) || (x.beginPosition == y.beginPosition && (x.endPosition) < (y.endPosition) );
	}
};

/*
 * les occurences sont triées d'abord par leur position de fin puis si egal, par leur position de début.
 */
template <typename TText> class LessFinalOccurenceEndPositionFirst
{
public:
	bool operator()(FinalOccurence<TText> const& x, FinalOccurence<TText> const& y){
		//XXX : les parenthèses sont necessaires pour compiler : x.beginPosition < y.beginPosition ne fonctionne pas.
		return (x.endPosition) < (y.endPosition) || (x.endPosition == y.endPosition && (x.beginPosition) < (y.endPosition) );
	}
};

//////////////////////TEST//////////////////////////


void extendLeft(AlignmentsMap & finalAligments, BwoloTIterator const&  it, Iterator<BwoloTPattern, Rooted>::Type const & patternIt, unsigned int seedPosition, unsigned int seedSize , unsigned short errNumber, bool ){
	ExtendingSeedState seed(patternIt, it, seedSize ,seedPosition,seedSize);//la seed commence avec une penalité de "seedSize" (normal)
	backtracking(finalAligments, seed, errNumber-seedPosition);//le nombre max d'erreur tiens compte du fait qu'une graine plus a droite ne doit pas exister
}


void extendLeft(AlignmentsMap & AligmentsByPattern, PotentialCandidates const& results,BwoloTPattern const & pattern, unsigned short errNumber, bool isVerbose){
	//préparation des iterateur dans le pattern
	BwoloTPattern localPattern = pattern;
	String<Iterator<BwoloTPattern, Rooted>::Type> patternIterators = getBeginSegmentPatternIterator(localPattern, errNumber);
	Iterator<BwoloTPattern, Rooted>::Type endPatternIterator = end(localPattern);
	appendValue(patternIterators, endPatternIterator);//pour le cas ou la graine est sur le dernier segment.
	String<String<candidate::Candidate> > beginPositionDependantResults;
	String<candidate::Candidate> sizeDependantResults;
	for(unsigned int seedPosition = 0 ; seedPosition < length(results); ++seedPosition){
		beginPositionDependantResults = results[seedPosition];
		for(unsigned int seedSize = 0 ; seedSize < length(beginPositionDependantResults); ++seedSize){
			sizeDependantResults = beginPositionDependantResults[seedSize];
			for(unsigned int seedCpt = 0 ; seedCpt < length(sizeDependantResults); ++seedCpt){
				candidate::Candidate candidateLocal = sizeDependantResults[seedCpt];
				Iterator<BwoloTPattern, Rooted>::Type patternIt=patternIterators[seedPosition+seedSize+2];
				extendLeft(AligmentsByPattern, candidateLocal.it,patternIt,seedPosition,seedSize , errNumber,isVerbose );
			}
		}
	}
}

void cleanResults(std::set<FinalOccurence<BwoloTText>, LessFinalOccurenceBeginPositionFirst<BwoloTText> > const & cleanResults,
		std::set<FinalOccurence<BwoloTText>, LessFinalOccurenceBeginPositionFirst<BwoloTText> > const & nonCleanResults){
	//je prend la premiere occurence o1
	//je regarde à k position plus loin. Si je trouve une occurence "erreur compatible" o2 :
	//	- moins bonne, je garde o1, je recommence avec l'occurence au moins k+1 eloigné
	//	- egal, je garde o1 et o2
	//	- meilleur : je prend o2, je verifie si elle est la meilleur dans sa zone "erreur compatible". je recommence avec
}



void extendRight(AlignmentsMap & AligmentsByPattern, BwoloTPattern const & pattern, BwoloTIndex const & index, unsigned short errNumber, bool /*isVerbose*/){
	BwoloTText const & text = indexText(index);
	typedef std::set<FinalOccurence<BwoloTText>, LessFinalOccurenceBeginPositionFirst<BwoloTText> > TSet;
	TSet mySet;

	for(auto it = AligmentsByPattern.begin() ; it !=AligmentsByPattern.end() ; ++it ){
		ExtendingSeedState seed = it->second;
		Infix<Fibre<BwoloTIndex, EsaSA>::Type const>::Type occurences = getOccurrences(seed.it);
		BwoloTText OccRepresentative = representative(seed.it);
		Size<BwoloTIndex>::Type OccLength = repLength(seed.it);

		for (unsigned int occPos = 0; occPos < length(occurences); occPos++){
			Position<BwoloTIndex>::Type occDBPosition = occurences[occPos];
			Position<BwoloTIndex>::Type OccRightPosition = occDBPosition+ OccLength;
			Position<BwoloTPattern>::Type OccPatternBeginPosition = endPosition(pattern) - getPatternSegmentBeginPos(pattern, seed.beginPatternSegment, errNumber); //le segment 0 est à droite dans le pattern "a l'endroit"
			unsigned short maxScore = errNumber - seed.score;
			BestAlignEndPos< Position<BwoloTText>::Type , unsigned short> endDBPos;
			if(extendRight(endDBPos, text, pattern, OccRightPosition, OccPatternBeginPosition,  maxScore) ){
				Position<BwoloTIndex>::Type occEndPosition = posNext(endDBPos.dbEndPos);
				unsigned short scoreFinal = endDBPos.score + seed.score;
//				BwoloTText occurence = infix(text, occDBPosition, occEndPosition);
//				cout<< occDBPosition <<":"<<occEndPosition<<":-"<< scoreFinal << ":"<< occurence << ":"<< pattern << endl;
				//TODO : ajouter le resultat dans un set
				mySet.insert(FinalOccurence<BwoloTText>(occDBPosition, occEndPosition,scoreFinal) );
			}
		}
	}
	//TODO filtration du set ?
	//TODO mettre l'affichage du set dans une fonction séparer.
	typedef TSet::iterator setIt;
	setIt mySetIt=mySet.begin();
	setIt endSetIt = mySet.end();
	while(mySetIt != endSetIt){
		FinalOccurence<BwoloTText> occurence = *mySetIt;
		BwoloTText occurenceSequence = infix(text, occurence.beginPosition, occurence.endPosition);
		cout<< occurence.beginPosition <<":"<< occurence.endPosition <<":-"<< occurence.score << ":"<< occurenceSequence << ":"<< pattern << endl;
		mySetIt++;
	}

}

////////////////////FIN TEST////////////////////////

static mySequence reverseComplement(mySequence seq) {
	mySequence seq_rev;
	seq_rev.id = seq.id;
	appendValue(seq_rev.id, '_');
	for (size_t i{0}; i < length(seq.seq); ++i) {
		char v = seq.seq[length(seq)-1-i];
		if (v == 'A') v = 'T';
		else if (v == 'a') v = 't';
		else if (v == 'C') v = 'G';
		else if (v == 'c') v = 'g';
		else if (v == 'G') v = 'C';
		else if (v == 'g') v = 'c';
		else if (v == 'T') v = 'A';
		else if (v == 't') v = 'a';

		appendValue(seq_rev.seq, v);
	}
	return seq_rev;
}

void getPatternsFromfasta(char* fastaFile, StringSet<mySequence> & patterns, bool addRevCompl, bool isVerbose){
	if (isVerbose) cout << "opening fasta file : ";
	SequenceStream seqStream(fastaFile);
	if (!isGood(seqStream)) {
		if (isVerbose) cout << "FAILED" << std::endl;
		if (isVerbose) cerr << "Error while opening the fasta file " << fastaFile << std::endl;
		exit(EXIT_FAILURE); // ouverture du fichier raté
	}
	if (isVerbose) cout << "OK"<<std::endl;
	if (isVerbose) cout << "reading sequences... " << std::endl;
	CharString id;
	BwoloTPattern seq;
	while (!atEnd(seqStream)){
		if (readRecord(id, seq, seqStream) != 0) {
			//Error : There was an error reading
			if (isVerbose) cout << "FAILED" << std::endl;
			if (isVerbose) cerr << "Error while reading sequence in fasta file " << fastaFile << std::endl;
			close(seqStream);
			exit(EXIT_FAILURE);
		} else {
			mySequence pattern;
			pattern.id = id;
			pattern.seq = seq;
			appendValue(patterns, pattern);
			if (isVerbose) cout << pattern.id << ':'<< pattern.seq << std::endl;

			if (addRevCompl) {
				mySequence pattern_rev = reverseComplement(pattern);
				appendValue(patterns, pattern_rev);
				if (isVerbose) cout << pattern_rev.id << ':'<< pattern_rev.seq << std::endl;
			}
		}
	}
	close(seqStream);
}

void getPatternsFromCStr(char* patternCstr, StringSet<mySequence> & patterns, bool addRevCompl, bool isVerbose){
	BwoloTPattern seq = patternCstr;
	CharString id = "input_Seq";
	mySequence pattern;
	pattern.id = id;
	pattern.seq = seq;
	appendValue(patterns, pattern);
	if (isVerbose) cout << pattern.id << ':'<< pattern.seq << std::endl;

	if (addRevCompl) {
		mySequence pattern_rev = reverseComplement(pattern);
		appendValue(patterns, pattern_rev);
		if (isVerbose) cout << pattern_rev.id << ':'<< pattern_rev.seq << std::endl;
	}
}

void printPotentialOccurences(unsigned int id, mySequence & pattern, PotentialCandidates & occurences, bool){
	String<String<candidate::Candidate> > beginPositionDependantResults;
	String<candidate::Candidate> sizeDependantResults;
	for(unsigned int i = 0 ; i < length(occurences); i++){
		beginPositionDependantResults = occurences[i];
		for(unsigned int j = 0 ; j < length(beginPositionDependantResults); j++){
			sizeDependantResults = beginPositionDependantResults[j];
			for(unsigned int k = 0 ; k < length(sizeDependantResults); k++){
				candidate::Candidate candidateLocal = sizeDependantResults[k];
				Infix<Fibre<BwoloTIndex, EsaSA>::Type const>::Type occurences =
						getOccurrences(candidateLocal.it);
				for (unsigned int occPos = 0; occPos < length(occurences); occPos++) {
					std::cout << id << ':' << pattern.id << ':' << pattern.seq //banalité sur l'id du pattern
							<< ':' << i << ':' << j //la position de début et la taille de la graine
							<< ':' << representative(candidateLocal.it) << ':' << occurences[occPos]<<std::endl;; // la valeur de l'occurence et sa position
				}
			}
		}
	}
}

void printPotentialOccurencesMultiple(StringSet<mySequence> & patterns, String<PotentialCandidates> & occurencesByPattern, bool isVerbose){
	size_t nbPattern = length(patterns);
	for (unsigned int i = 0; i < nbPattern; i++){
		printPotentialOccurences(i, patterns[i], occurencesByPattern[i], isVerbose);
	}
}

void shearchAndPrintPotentialOccurences_PatternString(StringSet<mySequence> & patterns, char* indexFile, unsigned short errNumber, bool isVerbose){

	if (isVerbose) cout<<"opening the index from file : " << indexFile;
	BwoloTIndex myIndex;
	if(!open(myIndex, indexFile)){
		if (isVerbose) cout<<"FAILED" << endl;
		cerr<< "Error while opening index from file " << indexFile <<std::endl;
		exit(EXIT_FAILURE);
	}
	if (isVerbose) cout<<"OK" << endl;

	if (isVerbose) cout<<"getting the root index iterator : ";
	BwoloTIterator indexRootIt(myIndex);
	if (isVerbose) cout<<"OK" << endl;

	if (isVerbose) cout<<"Filtering patterns..." << std::endl;
	size_t nbPattern = length(patterns);
	//String<PotentialCandidates> occurencesByPattern;

	clock_t start = clock();
	for (unsigned int i = 0; i < nbPattern; i++){
		BwoloTPattern pattern;
		pattern = patterns[i].seq;
		if (isVerbose) cout<<"pattern : "<< patterns[i].id << " : "<< pattern;
		reverse(pattern); //ne jamais oublier d'inverser le pattern pour un FM-index.
		StringSet<BwoloTPattern>  segmentedPattern;
		segmentedPattern = segmentPattern(pattern, errNumber);
		BwoloTIterator localIndexRoot =  indexRootIt;
		PotentialCandidates results;
		if (isVerbose) cout<<"... filtration...";
		PatternFiltrer::filtrate(results, localIndexRoot, segmentedPattern);
		if (isVerbose) cout<<"OK"<< std::endl;
		if (isVerbose) cout<<"... left extention...";
		AlignmentsMap AligmentsByPattern;
		extendLeft(AligmentsByPattern, results, pattern, errNumber, isVerbose);
		if (isVerbose) cout<<"OK"<<std::endl;
		if (isVerbose) cout<<"... right extention...";
		extendRight(AligmentsByPattern, patterns[i].seq, myIndex, errNumber, isVerbose);//attention, ici c'est le pattern dans le "bon sens"
		if (isVerbose) cout<<"OK"<<std::endl;
		//appendValue(occurencesByPattern, results);
	}
	clock_t diff = clock() - start;
	double msec = diff * 1. / CLOCKS_PER_SEC;
	std::cerr << "Search time: " << msec << "\n";
//	printPotentialOccurencesMultiple(patterns, occurencesByPattern, isVerbose);
}

void printPotentialOccurences_PatternSequences(char* patternSequence, char* indexFile, unsigned short errNumber, bool addRevCompl, bool isVerbose){
	StringSet<mySequence> patterns;
	if (isVerbose) cout << "get pattern from parameter" << std::endl;
	getPatternsFromCStr(patternSequence, patterns, addRevCompl, isVerbose);
	shearchAndPrintPotentialOccurences_PatternString(patterns, indexFile, errNumber, isVerbose);
}

void printPotentialOccurences_PatternFasta(char* patternFile, char* indexFile, unsigned short errNumber, bool addRevCompl, bool isVerbose){
	StringSet<mySequence> patterns;
	if (isVerbose) cout << "get pattern list from fasta file" << std::endl;
	getPatternsFromfasta(patternFile, patterns, addRevCompl, isVerbose);
	shearchAndPrintPotentialOccurences_PatternString(patterns, indexFile, errNumber,  isVerbose);
}






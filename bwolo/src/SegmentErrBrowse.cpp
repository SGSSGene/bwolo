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

#include "../include/SegmentErrBrowse.h"

namespace SegmentErrBrowse{
bool exactSegmentIterate(
		candidate::Candidate & candidate,
		BwoloTPattern const & patternSegment){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::exactSegmentIterate :"<< std::endl;
	std::cout << "\tcandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(candidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(candidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegment : " << patternSegment<< std::endl;
#endif
	bool rep = goDown(candidate.it, patternSegment);
	if(rep){
		candidate.firstNcOpSeg = ALL_NCOP_ALLOW;
	}
	return rep;
}

void filtrateExactSegment(
		String<candidate::Candidate> & candidateList,
		candidate::Candidate const actualPotentialCandidate,
		BwoloTPattern const & patternSegment
){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::filtrateExactSegment :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegment : " << patternSegment<< std::endl;
#endif
	candidate::Candidate localCandidate = actualPotentialCandidate;
	if(exactSegmentIterate(localCandidate,patternSegment)){
		appendValue(candidateList, localCandidate);
	}
}

bool deletionErrIterate(
		BwoloTIterator & candidateIt,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::deletionErrIterate :"<< std::endl;
	std::cout<<"\tcandidateIt : " << representative(candidateIt) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
#endif
	Iterator<BwoloTPattern, Rooted>::Type localPatternSegmentIt =
			patternSegmentIt;
	//j'avance de 1 dans le pattern = d??l??tion;
	goNext(localPatternSegmentIt);
	//on finit de descendre dans l'index
	while (!atEnd(localPatternSegmentIt)) {
		BwoloTIndexAlphabet ItValuePatternSeg = getValue(localPatternSegmentIt);
		if (!goDown(candidateIt, ItValuePatternSeg)) {
			return false;
		}
		goNext(localPatternSegmentIt);
	}
	return true;
}

bool insertionErrIterate(
		BwoloTIterator & candidateIt,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		BwoloTPatternAlphabet alphabetChar){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::insertionErrIterate :"<< std::endl;
	std::cout<<"\tcandidateIt : " << representative(candidateIt) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
	std::cout<<"\talphabetChar : " << alphabetChar << std::endl;
#endif

	Iterator<BwoloTPattern, Rooted>::Type localPatternSegmentIt =
			patternSegmentIt;
	//j'avance dans l'index avec le charact??re alphabetChar sans avancer dans le pattern si c'est possible
	if (!goDown(candidateIt, alphabetChar)) {
		return false;
	}
	for (; !atEnd(localPatternSegmentIt); goNext(localPatternSegmentIt)) {
		BwoloTPatternAlphabet ItValuePatternSeg = getValue(localPatternSegmentIt);
		if (!goDown(candidateIt, ItValuePatternSeg)) {
			return false;
		}
	}
	return true;
}

bool substitutionErrIterate(
		BwoloTIterator & candidateIt,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		BwoloTPatternAlphabet alphabetChar){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::substitutionErrIterate :"<< std::endl;
	std::cout<<"\tcandidateIt : " << representative(candidateIt) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
	std::cout<<"\talphabetChar : " << alphabetChar << std::endl;
#endif
	Iterator<BwoloTPattern, Rooted>::Type localPatternSegmentIt = patternSegmentIt;
	//je v??rifie que je suis bien dans un cas de substitution : que je ne remplace pas par le m??me caract??re;
	BwoloTPatternAlphabet ItValuePatternSeg = getValue(localPatternSegmentIt);
	if (ItValuePatternSeg != alphabetChar) {
		//j'avance dans l'index avec le charact??re alphabetChar si c'est possible et j'avance de 1 dans le pattern
		if (!goDown(candidateIt, alphabetChar)) {
			return false;
		}
		goNext(localPatternSegmentIt);
		//je fini d'avancer jusqu'?? la fin du pattern tant que c'est possible
		for (; !atEnd(localPatternSegmentIt); goNext(localPatternSegmentIt)) {
			BwoloTPatternAlphabet ItValuePatternSeg = getValue(localPatternSegmentIt);
			if (!goDown(candidateIt, ItValuePatternSeg)) {
				return false;
			}
		}
		return true;
	}else{
		//n'est pas une substitution : m??me charact??re
		return false;
	}
}

void browseErr(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool isSegLastNucleotide
){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::browseErr :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
	std::cout<<"\tisSegLastNucleotide : " << (isSegLastNucleotide?"true":"false") << std::endl;
#endif
	//deletion
	candidate::Candidate deletionCandidate = actualPotentialCandidate;
	if(deletionErrIterate(deletionCandidate.it, patternSegmentIt)){
		deletionCandidate.firstNcOpSeg = isSegLastNucleotide ? NO_INSERT : ALL_NCOP_ALLOW;
		appendValue(potentialCandidates,deletionCandidate);
	}
	browseErrWoDeletion(potentialCandidates, actualPotentialCandidate, patternSegmentIt, isSegLastNucleotide);
}

void browseErrWoDeletion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool
		){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::browseErrWoDeletion :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
#endif
	BwoloTPatternAlphabet alphabetChar;
	for (alphabetChar = MinValue<BwoloTPatternAlphabet>::VALUE;
			alphabetChar < +ValueSize<BwoloTPatternAlphabet>::VALUE;
			++alphabetChar) { //?? chaque caract??re, on introduit une erreur
		//insertion
		candidate::Candidate insertionCandidate = actualPotentialCandidate;
		if(insertionErrIterate(insertionCandidate.it, patternSegmentIt, alphabetChar)){
			//l'insertion ??tant r??aliser juste avant le dernier, pas de restriction pour la prochaine section
			insertionCandidate.firstNcOpSeg = ALL_NCOP_ALLOW;
			appendValue(potentialCandidates,insertionCandidate);
		}
		//substitution
		candidate::Candidate substitutionCandidate = actualPotentialCandidate;
		if(substitutionErrIterate(substitutionCandidate.it, patternSegmentIt, alphabetChar)){
			substitutionCandidate.firstNcOpSeg = ALL_NCOP_ALLOW;
			appendValue(potentialCandidates,substitutionCandidate);
		}
	}
}

void browseErrWoInsertion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool isSegLastNucleotide
		){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::browseErrWoInsertion :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
	std::cout<<"\tisSegLastNucleotide : " << (isSegLastNucleotide?"true":"false") << std::endl;
#endif
	//deletion
	candidate::Candidate deletionCandidate = actualPotentialCandidate;
	if(deletionErrIterate(deletionCandidate.it, patternSegmentIt)){
		deletionCandidate.firstNcOpSeg = isSegLastNucleotide ? NO_INSERT : ALL_NCOP_ALLOW;
		appendValue(potentialCandidates,deletionCandidate);
	}
	BwoloTPatternAlphabet alphabetChar;
	for (alphabetChar = MinValue<BwoloTPatternAlphabet>::VALUE;
			alphabetChar < +ValueSize<BwoloTPatternAlphabet>::VALUE;
			++alphabetChar) { //?? chaque caract??re, on introduit une erreur
		//substitution
		candidate::Candidate substitutionCandidate = actualPotentialCandidate;
		if(substitutionErrIterate(substitutionCandidate.it, patternSegmentIt, alphabetChar)){
			substitutionCandidate.firstNcOpSeg = ALL_NCOP_ALLOW;
			appendValue(potentialCandidates,substitutionCandidate);
		}
	}
}

void browseLastInsertion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate
		){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::browseLastInsertion :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
#endif
	BwoloTPatternAlphabet alphabetChar;
	for (alphabetChar = MinValue<BwoloTPatternAlphabet>::VALUE;
			alphabetChar < +ValueSize<BwoloTPatternAlphabet>::VALUE;
			++alphabetChar) { //?? chaque caract??re, on introduit une erreur
		//insertion
		candidate::Candidate insertionCandidate = actualPotentialCandidate;
		if(goDown(insertionCandidate.it, alphabetChar)){
			insertionCandidate.firstNcOpSeg = NO_DELET;
			appendValue(potentialCandidates,insertionCandidate);
		}
	}
}


void browseErrFirstInsertion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool
		){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::browseErrFirstInsertion :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
#endif
	BwoloTPatternAlphabet alphabetChar;
	for (alphabetChar = MinValue<BwoloTPatternAlphabet>::VALUE;
			alphabetChar < +ValueSize<BwoloTPatternAlphabet>::VALUE;
			++alphabetChar) { //?? chaque caract??re, on introduit une erreur
		//insertion
		candidate::Candidate insertionCandidate = actualPotentialCandidate;
		if(insertionErrIterate(insertionCandidate.it, patternSegmentIt, alphabetChar)){
			//l'insertion ??tant r??aliser juste avant le dernier, pas de restriction pour la prochaine section
			insertionCandidate.firstNcOpSeg = ALL_NCOP_ALLOW;
			appendValue(potentialCandidates,insertionCandidate);
		}
	}
}

void browseErrFirstDeletion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool isSegLastNucleotide
		){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::browseErrFirstDeletion :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
	std::cout<<"\tisSegLastNucleotide : " << (isSegLastNucleotide?"true":"false") << std::endl;
#endif
	candidate::Candidate deletionCandidate = actualPotentialCandidate;
	if(deletionErrIterate(deletionCandidate.it, patternSegmentIt)){
		deletionCandidate.firstNcOpSeg = isSegLastNucleotide ? NO_INSERT : ALL_NCOP_ALLOW;
		appendValue(potentialCandidates,deletionCandidate);
	}
}

void browseErrFirstSubstitution(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool
		){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::browseErrFirstSubstitution :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegmentIt : " << value(patternSegmentIt)<< std::endl;
#endif
	BwoloTPatternAlphabet alphabetChar;
	for (alphabetChar = MinValue<BwoloTPatternAlphabet>::VALUE;
			alphabetChar < +ValueSize<BwoloTPatternAlphabet>::VALUE;
			++alphabetChar) { //?? chaque caract??re, on introduit une erreur
		//substitution
		candidate::Candidate substitutionCandidate = actualPotentialCandidate;
		if(substitutionErrIterate(substitutionCandidate.it, patternSegmentIt, alphabetChar)){
			substitutionCandidate.firstNcOpSeg = ALL_NCOP_ALLOW;
			appendValue(potentialCandidates,substitutionCandidate);
		}
	}
}

void filtrateFirstNcSegmentWithError(
		String<candidate::Candidate> & candidateList,
		candidate::Candidate const & actualPotentialCandidate,
		NcOp::NcOp operationAllowed,
		BwoloTPattern const & patternSegment
){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::filtrateFirstNcSegmentWithError :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\toperationAllowed : " << int(operationAllowed) << std::endl;
	std::cout<<"\tpatternSegment : " << patternSegment << std::endl;
#endif
	candidate::Candidate localCandidate = actualPotentialCandidate;
	//BwoloTPatternAlphabet alphabetChar;
	BwoloTPattern localPatternSegment = patternSegment;
	Iterator<BwoloTPattern, Rooted>::Type nextPatternSegmentIt = begin(localPatternSegment);
	Iterator<BwoloTPattern, Rooted>::Type patternSegmentIt = nextPatternSegmentIt;
	if(!atEnd(nextPatternSegmentIt)){
		goNext(nextPatternSegmentIt);
	}
	bool lastNc = atEnd(nextPatternSegmentIt);
	//premier Nucl??otide
	if(IS_INSERT_ALLOW(operationAllowed)){
		//browseInsertion
		browseErrFirstInsertion(candidateList, localCandidate, patternSegmentIt, lastNc);
	}
	if(IS_DELET_ALLOW(operationAllowed)){
		//browseDeletion
		browseErrFirstDeletion(candidateList, localCandidate, patternSegmentIt, lastNc);
	}
	if(IS_SUBSTIT_ALLOW(operationAllowed)){
		//browseSubstitution
		browseErrFirstSubstitution(candidateList, localCandidate, patternSegmentIt, lastNc);
	}
}

void filtrateSegmentWithError(
		String<candidate::Candidate> & candidateList,
		candidate::Candidate const & actualPotentialCandidate,
		BwoloTPattern const & patternSegment
){
#ifdef DEBUG_VERBOSE
	std::cout <<"SegmentErrBrowse::filtrateSegmentWithError :"<< std::endl;
	std::cout << "\tactualPotentialCandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(actualPotentialCandidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(actualPotentialCandidate.firstNcOpSeg) << std::endl;
	std::cout<<"\tpatternSegment : " << patternSegment << std::endl;
#endif
	candidate::Candidate localCandidate = actualPotentialCandidate;
	BwoloTPatternAlphabet alphabetChar;
	BwoloTPattern localPatternSegment = patternSegment;
	Iterator<BwoloTPattern, Rooted>::Type nextPatternSegmentIt = begin(localPatternSegment);
	Iterator<BwoloTPattern, Rooted>::Type patternSegmentIt = nextPatternSegmentIt;
	if(!atEnd(nextPatternSegmentIt)){
		goNext(nextPatternSegmentIt);
	}
	bool lastNc = atEnd(nextPatternSegmentIt);
	//premier Nucl??otide
	switch(actualPotentialCandidate.firstNcOpSeg){
	case NO_DELET :
		//browseWOdeletion
		browseErrWoDeletion(candidateList, localCandidate, patternSegmentIt, lastNc);
		break;
	case NO_INSERT :
		//browseWOinsertion
		browseErrWoInsertion(candidateList, localCandidate, patternSegmentIt, lastNc);
		break;
	default :
		//browseErr
		browseErr(candidateList, localCandidate, patternSegmentIt, lastNc);
		break;
	}
	//passage du premier nucl??otide
	alphabetChar = getValue(patternSegmentIt);
	//TODO : v??rifier si une insertion n'est pas possible
	if(!goDown(localCandidate.it, alphabetChar)){ // si je ne peux pas descendre de mani??re exact == stop ici;
		return;
	}
	//nucl??otides centraux
	patternSegmentIt = nextPatternSegmentIt;
	if(!lastNc){
		goNext(nextPatternSegmentIt);
	}
	while(!atEnd(nextPatternSegmentIt)){
		alphabetChar = getValue(patternSegmentIt);
		//browseErr
		browseErr(candidateList, localCandidate, patternSegmentIt, false);
		//passage au nucl??otide suivant
		if(!goDown(localCandidate.it, alphabetChar)){ // si je ne peux pas descendre de mani??re exact == stop ici;
			return;
		}
		patternSegmentIt = nextPatternSegmentIt;
		goNext(nextPatternSegmentIt);
	}
	//dernier Nucl??otide
	//browseErr
	browseErr(candidateList, localCandidate, patternSegmentIt, true);
	//lastInsertion
	alphabetChar = getValue(patternSegmentIt);
	//passage du dernier nucl??otide
	if(!goDown(localCandidate.it, alphabetChar)){ // si je ne peux pas descendre de mani??re exact == stop ici;
		return;
	}
	//derniere erreur possible : apr??s le dernier nucl??otide
	browseLastInsertion(candidateList, localCandidate);
}

}




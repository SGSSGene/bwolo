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



#include "../include/PatternFiltrer.h"

#include "../include/CandidateLog.h"
#include "../include/SegmentErrBrowse.h"

namespace PatternFiltrer{
void filtrateFromPatternSegmentPosition(
			String<String<candidate::Candidate> > & positionDependantResults,
			BwoloTIterator const & indexRootIt,
			Iterator<StringSet<BwoloTPattern> >::Type const & patternSegmentIt){
#ifdef DEBUG_VERBOSE
	std::cout <<"PatternFiltrer::filtrateFromPatternSegmentPosition :"<< std::endl;
	std::cout << "\tindexRootIt : "<< representative(indexRootIt) << std::endl;
	std::cout<< "\tpatternSegmentIt : " << value(patternSegmentIt) <<std::endl;
#endif
	//Mise en place des variables locales
	candidateLog::CandidateLogger logger;
	String<candidate::Candidate> potentialCandidates;
	String<candidate::Candidate> nextSegmentPotentialCandidates;
	Iterator<StringSet<BwoloTPattern> >::Type localPatternSegmentIt;
	Iterator<StringSet<BwoloTPattern> >::Type nextpatternSegmentIt = patternSegmentIt;
	//prépartion du hasNext
	localPatternSegmentIt = nextpatternSegmentIt;
	if(!atEnd(nextpatternSegmentIt)){
		goNext(nextpatternSegmentIt);
	}
	BwoloTPattern patternSegmentValue = value(localPatternSegmentIt);
	candidate::Candidate firstCandidate;

	//préparartion de l'iteration : premier candidate
	firstCandidate.it = indexRootIt;
	firstCandidate.firstNcOpSeg = ALL_NCOP_ALLOW;
	//iterration du premier segment exact et ajout dans la liste des prochains candidats potentiels à traiter
	SegmentErrBrowse::filtrateExactSegment(nextSegmentPotentialCandidates, firstCandidate, patternSegmentValue);
	localPatternSegmentIt = nextpatternSegmentIt;
	if(!atEnd(nextpatternSegmentIt)){
		goNext(nextpatternSegmentIt);
	}
	//iteration sur chaque segment du pattern
	while(!atEnd(localPatternSegmentIt)){

		//préparation de l'iteration
		clear(potentialCandidates); //récupère la nouvelle liste de candidats potentiels
		potentialCandidates = nextSegmentPotentialCandidates;
		clear(nextSegmentPotentialCandidates);
		Iterator<String<candidate::Candidate> , Rooted>::Type candidatesIt = begin(potentialCandidates);
		candidateLog::clear(logger);//réinitialise le logger des candidats potentiels déja traité.
		String<candidate::Candidate> sizeDependantResults;//liste des candidats filtrés pour cette longueur de graine
		patternSegmentValue = value(localPatternSegmentIt);
		//itération de la liste des canditats
		while(!atEnd(candidatesIt)){
			candidate::Candidate myCandidate = value(candidatesIt);
			NcOp::NcOp remainsOps = NO_NCOP_ALLOW;
			if(!candidateLog::contains(remainsOps, logger, myCandidate)){ // on ne traite un candidate seulement si un homonyme avec des conditions similaires n'a pas déjà été traité.
				candidateLog::add(logger, myCandidate);
				if(!atEnd(nextpatternSegmentIt)){ //filtration du segment avec 1 erreur si il y a encore au moins un segment après celui qui est traité --> futurs candidats potentiels
					SegmentErrBrowse::filtrateSegmentWithError(nextSegmentPotentialCandidates, myCandidate, patternSegmentValue);
				}else if(remainsOps != NO_NCOP_ALLOW){//si il y a encore des opérations autorisées
					SegmentErrBrowse::filtrateFirstNcSegmentWithError(nextSegmentPotentialCandidates, myCandidate, remainsOps, patternSegmentValue);
					candidateLog::update(logger, myCandidate, remainsOps);
				}
				//cloture de la graines par un segment exact --> candidats filtrés pour cette longueur de graine
				SegmentErrBrowse::filtrateExactSegment(sizeDependantResults, myCandidate, patternSegmentValue);
			}
			goNext(candidatesIt);
		}
		appendValue(positionDependantResults, sizeDependantResults);
		localPatternSegmentIt = nextpatternSegmentIt;
		if(!atEnd(nextpatternSegmentIt)){
			goNext(nextpatternSegmentIt);
		}
	}
}

void filtrate(
		String<String<String<candidate::Candidate> > > & results,
		BwoloTIterator const & indexRootIt,
		StringSet<BwoloTPattern> const & pattern){
#ifdef DEBUG_VERBOSE
	std::cout <<"PatternFiltrer::filtrate :"<< std::endl;
	std::cout << "\tindexRootIt : " <<  "^" << representative(indexRootIt)<<std::endl;
	std::cout <<"\tpattern : " << "?" <<std::endl;
#endif
	StringSet<BwoloTPattern> localPattern = pattern;
	Iterator<StringSet<BwoloTPattern>, Rooted  >::Type patternSegmentIt;
	Iterator<StringSet<BwoloTPattern>, Rooted >::Type nextpatternSegmentIt = begin(localPattern);
	//prépartion du hasNext
	patternSegmentIt = nextpatternSegmentIt;
	if(!atEnd(nextpatternSegmentIt)){
		goNext(nextpatternSegmentIt);
	}
#ifdef DEBUG_WHILE_VERBOSE
		unsigned whileCpt1 = 0;
#endif
	while(!atEnd(nextpatternSegmentIt)){
#ifdef DEBUG_WHILE_VERBOSE
		std::cout<< "whileCpt1 = " << whileCpt1<<std::endl;
		whileCpt1++;
#endif
		String<String<candidate::Candidate> > positionDependantResults;
		filtrateFromPatternSegmentPosition( positionDependantResults, indexRootIt, patternSegmentIt);
		appendValue(results, positionDependantResults);
		patternSegmentIt = nextpatternSegmentIt;
		goNext(nextpatternSegmentIt);
	}
}
}

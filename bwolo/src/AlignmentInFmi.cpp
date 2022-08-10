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

#include "../include/AlignmentInFmi.h"
using namespace seqan;

AligmentPosition::AligmentPosition( AligmentPosition const& k):
		patternIt(k.patternIt),
		beginPatternSegment(k.beginPatternSegment),
		it(k.it){}

AligmentPosition::AligmentPosition(Iterator<BwoloTPattern, Rooted>::Type const& patternIt, unsigned short beginPatternSegment, BwoloTIterator const& it):
		patternIt(patternIt),
		beginPatternSegment(beginPatternSegment),
		it(it){}

void AligmentPosition::init(Iterator<BwoloTPattern, Rooted>::Type const& patternIt, unsigned short beginPatternSegment, BwoloTIterator const& it)
{
		this->patternIt = Iterator<BwoloTPattern, Rooted>::Type(patternIt);
		this->beginPatternSegment = beginPatternSegment;
		this->it = BwoloTIterator(it);
	}

ExtendingSeedState::ExtendingSeedState(ExtendingSeedState const& k) :
						AligmentPosition(k.patternIt,k.beginPatternSegment,k.it),
						score(k.score),
						initialSeedCoreSize(k.initialSeedCoreSize){}

ExtendingSeedState::ExtendingSeedState(Iterator<BwoloTPattern, Rooted>::Type const& patternIt,
			BwoloTIterator const&  it,
			unsigned short score,
			unsigned short beginPatternSegment,
			unsigned short initialSeedCoreSize) :
				AligmentPosition(patternIt,beginPatternSegment,it),
				score(score),
				initialSeedCoreSize(initialSeedCoreSize){}

void ExtendingSeedState::init(unsigned short score,
	unsigned short initialSeedCoreSize){
		this->score = score;
		this->initialSeedCoreSize = initialSeedCoreSize;
	}

AlignmentMapEntry newAlignmentMapEntry(ExtendingSeedState seed){
	return AlignmentMapEntry(seed, seed);
}

//insert l'AligmentMap la nouvelle entrée ou met a jour cette dernière avec un meilleur score si c'est le cas.
// return si l'entrée a été ajouté ou mise a jour.
bool AlignmentMapInsert(AlignmentsMap & map, AlignmentMapEntry entry){
	AlignmentsMap::iterator repIt = map.find(entry.first);
	if(repIt != map.end()){
		if (repIt->second.score > entry.second.score){
			repIt->second.score = entry.second.score;
			repIt->second.initialSeedCoreSize = entry.second.initialSeedCoreSize;
			return true;
		}else{
			return false;
		}
	}else{
		map.insert(entry);
		return true;
	}
}

void backtrakingNextStep(AlignmentsMap & finalAligments,
		AlignmentsMap & tPlusOneAlignments,
		AlignmentsMap &  tPlusTwoAlignments,
		ExtendingSeedState const & seedState,
		unsigned short maxScore){
	if(atEnd(seedState.patternIt)){
		AlignmentMapInsert(finalAligments, newAlignmentMapEntry(seedState));
		return;
	}
	else{
		Iterator<BwoloTPattern, Rooted>::Type nextPatternIt = seedState.patternIt;
		BwoloTPatternAlphabet patternChar = getValue(nextPatternIt);
		bool matchOnly = true;
		goNext(nextPatternIt);
		if(seedState.score < maxScore){//déplacement avec erreurs possible
			for (BwoloTPatternAlphabet alphabetChar = MinValue<BwoloTPatternAlphabet>::VALUE;
								alphabetChar < +ValueSize<BwoloTPatternAlphabet>::VALUE;
								++alphabetChar){
				BwoloTIterator nextIndexIt = seedState.it;
				if(goDown(nextIndexIt, alphabetChar)){
					if(patternChar != alphabetChar){
						matchOnly = false;
						//parcours par insertion;
						ExtendingSeedState insertionSeed(seedState);
						insertionSeed.score++;
						insertionSeed.it = nextIndexIt;
						AlignmentMapInsert(tPlusOneAlignments, newAlignmentMapEntry(insertionSeed));
						//parcours par substitution
						ExtendingSeedState substitutionSeed(seedState) ;
						substitutionSeed.score++;
						substitutionSeed.it = nextIndexIt;
						substitutionSeed.patternIt = nextPatternIt;
						AlignmentMapInsert(tPlusTwoAlignments, newAlignmentMapEntry(substitutionSeed));
					}else{
						//parcours par match
						ExtendingSeedState matchSeed(seedState) ;
						matchSeed.it = nextIndexIt;
						matchSeed.patternIt = nextPatternIt;
						AlignmentMapInsert(tPlusTwoAlignments, newAlignmentMapEntry(matchSeed));
					}
				}
			}
			if(!matchOnly){
				//parcours par deletion
				ExtendingSeedState deletionSeed(seedState) ;
				deletionSeed.score++;
				deletionSeed.patternIt = nextPatternIt;
				AlignmentMapInsert(tPlusOneAlignments, newAlignmentMapEntry(deletionSeed));
			}
		}else{
			//parcours par match uniquement
			BwoloTIterator localIt = seedState.it;
			if(goDown(localIt, patternChar)){
				ExtendingSeedState matchSeed = seedState ;
				matchSeed.it = localIt;
				matchSeed.patternIt = nextPatternIt;
				AlignmentMapInsert (tPlusTwoAlignments, newAlignmentMapEntry(matchSeed));
			}
		}
	}
}

void backtracking(AlignmentsMap & finalAligments,
		ExtendingSeedState const&  seed,
		unsigned short maxScore){
	AlignmentsMap alignments;
	AlignmentsMap tPlusOneAlignments;//ici seront insérées les insertions et les déletions
	AlignmentsMap tPlusTwoAlignments;// ici seront insérés les matches et les substitutions
	ExtendingSeedState localSeed(seed);
	AlignmentMapInsert(alignments, newAlignmentMapEntry(ExtendingSeedState(localSeed)));
	while (! (alignments.empty() && tPlusOneAlignments.empty())){
		for(auto alignIt = alignments.begin(); alignIt != alignments.end(); ++alignIt){
			ExtendingSeedState seedState = alignIt->second;
			backtrakingNextStep(finalAligments, tPlusOneAlignments, tPlusTwoAlignments, seedState, maxScore);
		}
		alignments = tPlusOneAlignments;
		tPlusOneAlignments = tPlusTwoAlignments;
		tPlusTwoAlignments.clear();
	}
}
void printResult(AlignmentsMap & finalAligments){
	for(auto it = finalAligments.begin() ; it !=finalAligments.end() ; ++it ){
		ExtendingSeedState seed = it->second;
		Infix<Fibre<BwoloTIndex, EsaSA>::Type const>::Type occurences = getOccurrences(seed.it);
		BwoloTText OccRepresentative = representative(seed.it);
		unsigned int OccLength = repLength(seed.it);
		for (unsigned int occPos = 0; occPos < length(occurences); occPos++) {
			std::cout << occurences[occPos]<< ':'  << occurences[occPos] + OccLength << ':' << (-seed.score) << ':' << OccRepresentative << std::endl; // la valeur de l'occurence et sa position
		}
	}
}



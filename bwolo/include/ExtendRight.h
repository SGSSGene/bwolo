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

#ifndef EXTENDRIGHT_H_
#define EXTENDRIGHT_H_

#include <algorithm>

#include <seqan/sequence.h>
using namespace seqan;

template <typename TPosition, typename TScore = unsigned short>
struct BestAlignEndPos{
	TPosition dbEndPos;
	TScore score;
};

template <typename T> inline T MaxV(){
	return std::numeric_limits<T>::max();
}

template <typename TText, typename TPosition, typename TAlphabet, typename TScoreMatrice>
	void fillMatriceLine(
		BestAlignEndPos<TPosition, typename Value<TScoreMatrice>::Type > & bestAlignEndPos,
		TScoreMatrice & matriceLine,
		TScoreMatrice const& previousLine,
		TAlphabet PatternValue,
		TText const& dataBase,
		TPosition dbPos,
		TPosition dbEnd){
	typedef typename Value<TScoreMatrice>::Type TScore;
	typedef typename Position<TScoreMatrice>::Type TMatricePos;

//	const TScore NAR = std::numeric_limits<TScore>::max();
	bestAlignEndPos.dbEndPos = dbPos;
	bestAlignEndPos.score = MaxV<TScore>();
	TMatricePos matriceEndPos = endPosition(matriceLine);
	TScore score = (matriceEndPos-1)/2;
	dbPos -= score;//astuce pour synchronyser le début de la matrice avec dbPos; Sinon dbPos correspond au centre de la matrice.
	TScore insertionScore = MaxV<TScore>();
	TScore deletionScore = MaxV<TScore>();
	TScore matchSubstition = MaxV<TScore>();
	TScore bestValue = MaxV<TScore>();

	for(TMatricePos cpt = 0 ; cpt<matriceEndPos; cpt++){
		//replissage d'une case. Attention avec les cases "NAR".
		if(cpt != 0 ){
			insertionScore= (matriceLine[cpt-1]== MaxV<TScore>()) ? MaxV<TScore>() : (matriceLine[cpt-1])+1;
		}else{
			insertionScore = MaxV<TScore>();
		}
		if((cpt+1)  < matriceEndPos){
			deletionScore = previousLine[cpt+1] == MaxV<TScore>() ? MaxV<TScore>() : (previousLine[cpt+1])+1;
		}else{
			deletionScore = MaxV<TScore>();
		}
		if(dbPos+cpt < dbEnd && dbPos+cpt>=0){ //cas ou j'essaie d'alligner avec la fin du texte et cas (faible risque) ou avec le début du texte (impossible dans le cas de bwolo)
			matchSubstition = dataBase[dbPos+cpt] == PatternValue ?
					previousLine[cpt] :
					(previousLine[cpt]==MaxV<TScore>() ? MaxV<TScore>() : previousLine[cpt]+1);
		}else{
			matchSubstition = MaxV<TScore>();
		}
		bestValue = std::min(insertionScore, std::min(deletionScore, matchSubstition));
		matriceLine[cpt]=bestValue;
		if(bestAlignEndPos.score > bestValue){
			bestAlignEndPos.score = bestValue;
			bestAlignEndPos.dbEndPos = dbPos+cpt;
		}
	}
}


template <typename TText,typename TPattern , typename TScore>
	bool extendRight(
		BestAlignEndPos<typename Position<TText >::Type, TScore> & bestAlignEndPos,
		TText const & dataBase,
		TPattern const& pattern,
		typename Position<TText >::Type dbBeginPos,
		typename Position<TPattern >::Type patternBeginPos,
		TScore maxScore){
	typedef typename Position<TText >::Type TPosition1;
	typedef typename Position<TPattern >::Type TPosition2;


	String<TScore> matriceLine; //Score pour un segment du text necessaire pour l'alignement à une position dans le pattern
	resize(matriceLine, maxScore*2+1, MaxV<TScore>());
	String<TScore> previousLine = matriceLine;//position précédente dans le pattern.
	for(int cpt = 0 ; cpt <= maxScore ; cpt++){
		previousLine[cpt+maxScore]=cpt; // initialisation previousLine = [X:...:X:X:X:0:1:2:3:...:maxScore] (previousLine[maxScore] = 0)
	}
	TPosition1 dbPos = dbBeginPos;
	TPosition1 dbEnd = endPosition(dataBase);
	TPosition2 patternPos = patternBeginPos;
	TPosition2 patternEnd = endPosition(pattern);
	//a partir de maintenant : MatriceAlignement[dbPos][patternPos] = matriceLine[maxScore];
	//ne pas oublier d'initialiser correctement bestAlignEndPos si patternPos < patternEnd.
	bestAlignEndPos.dbEndPos = dbPos-1;
	bestAlignEndPos.score = 0;
	while(patternPos < patternEnd){
		typename Value<TPattern>::Type PatternValue = pattern[patternPos];
		fillMatriceLine(bestAlignEndPos, matriceLine, previousLine, PatternValue, dataBase, dbPos, dbEnd);
		if(bestAlignEndPos.score > maxScore){
			break;
		}
		patternPos++;
		dbPos++;
		previousLine = matriceLine;
	}
	return bestAlignEndPos.score <= maxScore;
}


#endif /* EXTENDRIGHT_H_ */

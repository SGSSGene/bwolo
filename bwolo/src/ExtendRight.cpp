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

//#include "../include/ExtendRight.h"
//
//#include <algorithm>
//
//#include <iostream>
//
//template<typename TValue> BestAlignEndPos fillMatriceLine(String<unsigned short> & matriceLine,
//		String<unsigned short> const& previousLine,
//		TValue PatternValue,
//		BwoloTText const& dataBase,
//		Size<BwoloTPattern>::Type dbPos,
//		Size<BwoloTPattern>::Type dbEnd){
//	const unsigned short NAR = std::numeric_limits<unsigned short>::max();
//	BestAlignEndPos bestPos;
//	bestPos.dbEndPos = dbPos;
//	bestPos.score = NAR;
//	unsigned int matriceLength = length(matriceLine);
//	unsigned int score = (matriceLength-1)/2;
//	dbPos -= score;//astuce pour synchronyser le début de la matrice avec dbPos; Sinon dbPos correspond au centre de la matrice.
//	unsigned short insertionScore = NAR;
//	unsigned short deletionScore = NAR;
//	unsigned short matchSubstition = NAR;
//	unsigned short bestValue = NAR;
//
//	for(unsigned int cpt = 0 ; cpt<matriceLength; cpt++){
//		//replissage d'une case. Attention avec les cases "NAR".
//		if(cpt != 0 ){
//			insertionScore= (matriceLine[cpt-1]== NAR) ? NAR : (matriceLine[cpt-1])+1;
//		}else{
//			insertionScore = NAR;
//		}
//		if((cpt+1)  < matriceLength){
//			deletionScore = previousLine[cpt+1] == NAR ? NAR : (previousLine[cpt+1])+1;
//		}else{
//			deletionScore = NAR;
//		}
//		if(dbPos+cpt < dbEnd && dbPos+cpt>=0){ //cas ou j'essaie d'alligner avec la fin du texte et cas (faible risque) ou avec le début du texte (impossible dans le cas de bwolo)
//			matchSubstition = dataBase[dbPos+cpt] == PatternValue ?
//					previousLine[cpt] :
//					(previousLine[cpt]==NAR ? NAR : previousLine[cpt]+1);
//		}else{
//			matchSubstition = NAR;
//		}
//		bestValue = std::min(insertionScore, std::min(deletionScore, matchSubstition));
//		matriceLine[cpt]=bestValue;
//		if(bestPos.score > bestValue){
//			bestPos.score = bestValue;
//			bestPos.dbEndPos = dbPos+cpt;
//		}
//	}
//	return bestPos;
//}
//
//
//BestAlignEndPos extendRight(BwoloTText const & dataBase,
//		BwoloTPattern const& pattern,
//		Size<BwoloTText>::Type dbBeginPos,
//		Size<BwoloTPattern>::Type patternBeginPos,
//		unsigned short maxScore){
//
////	std::cout<< "pattern : "<<suffix(pattern, patternBeginPos) << std::endl;
////	std::cout<< "text : " << infix(dataBase, dbBeginPos, dbBeginPos+20) << std::endl;
//
//	typedef typename Value<BwoloTPatternIterator>::Type TValue;
//	const short NAR = std::numeric_limits<unsigned short>::max();
//
//	String<unsigned short> matriceLine;
//	resize(matriceLine, maxScore*2+1, NAR);
//	String<unsigned short> previousLine = matriceLine;
//	for(int cpt = 0 ; cpt <= maxScore ; cpt++){
//		previousLine[cpt+maxScore]=cpt;
//	}
//	///////////////////////test//////////////////
////	for(unsigned cptLength = 0 ; cptLength < length(previousLine); cptLength++){
////		std::cout << previousLine[cptLength] << ":";
////	}
////	std::cout << std::endl;
///////////////////////////////////////////////////
//
//	Size<BwoloTText>::Type dbPos = dbBeginPos;
//	Size<BwoloTText>::Type dbEnd = length(dataBase);
//	Size<BwoloTPattern>::Type patternPos = patternBeginPos;
//	Size<BwoloTPattern>::Type patternEnd = length(pattern);
//	BestAlignEndPos bestAlign;
//	while(patternPos < patternEnd){
//		TValue PatternValue = pattern[patternPos];
//		bestAlign = fillMatriceLine(matriceLine, previousLine, PatternValue, dataBase, dbPos, dbEnd);
//////////////////////// test ////////////////////////
////		for(unsigned cptLength = 0 ; cptLength < length(matriceLine); cptLength++){
////			std::cout << matriceLine[cptLength] << ":";
////		}
////		std::cout << std::endl;
//////////////////////////////////////////////////////
//		if(bestAlign.score > maxScore){
//			break;
//		}
//		patternPos++;
//		dbPos++;
//
//		previousLine = matriceLine;
//	}
//	return bestAlign;
//}
//
//
//
//

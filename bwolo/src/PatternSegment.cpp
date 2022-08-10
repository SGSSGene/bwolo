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

#include "../include/PatternSegment.h"

Segment<BwoloTPattern>
getPatternSegment(BwoloTPattern & pattern, unsigned short segmentPosition, unsigned short errNumber)
{
#ifdef DEBUG_VERBOSE
	std::cout <<"getPatternSegment :"<< std::endl;
	std::cout << "\tgetPatternSegment : " << pattern << std::endl;
	std::cout << "\tsegmentPosition : " << segmentPosition << std::endl;
#endif
	unsigned int beginPos = getPatternSegmentBeginPos(pattern, segmentPosition, errNumber);
	unsigned int endPos = getPatternSegmentEndPos(pattern, segmentPosition, errNumber);
	Segment<BwoloTPattern> seg = infix(pattern, beginPos, endPos);
	return seg;
}

StringSet<BwoloTPattern> segmentPattern(BwoloTPattern & pattern, unsigned short errNumber){
#ifdef DEBUG_VERBOSE
	std::cout <<"segmentPattern :"<< std::endl;
	std::cout << "\tpattern : " << pattern << std::endl;
#endif

	StringSet<BwoloTPattern> patternSegmentate;
	BwoloTPattern patternSegment;
	unsigned short patternSegmentNumber = errNumber+2;
	for(unsigned short i(0); i< patternSegmentNumber; i++){
		patternSegment = getPatternSegment(pattern, i, errNumber);
		appendValue(patternSegmentate,patternSegment );
	}

	return patternSegmentate;
}


String<Iterator<BwoloTPattern, Rooted>::Type> getBeginSegmentPatternIterator(BwoloTPattern & pattern, unsigned short errNumber){
	String<Iterator<BwoloTPattern, Rooted>::Type> result;
	unsigned short patternSegmentNumber = errNumber+2;
	for(unsigned short segmentPosition(0); segmentPosition< patternSegmentNumber; segmentPosition++){
		unsigned int beginPos = getPatternSegmentBeginPos(pattern, segmentPosition, errNumber);
		Iterator<BwoloTPattern, Rooted>::Type myIter = iter(pattern, beginPos, Rooted());
		appendValue(result, myIter);
	}
	return result;
}

String<Iterator<BwoloTPattern, Rooted>::Type> getEndSegmentPatternIterator(BwoloTPattern & pattern, unsigned short errNumber){
	String<Iterator<BwoloTPattern, Rooted>::Type> result;
	unsigned short patternSegmentNumber = errNumber+2;
	for(unsigned short segmentPosition(0); segmentPosition< patternSegmentNumber; segmentPosition++){
		unsigned int endPos = getPatternSegmentEndPos(pattern, segmentPosition, errNumber);
		Iterator<BwoloTPattern, Rooted>::Type myIter = iter(pattern, endPos, Rooted());
		appendValue(result, myIter);
	}
	return result;
}






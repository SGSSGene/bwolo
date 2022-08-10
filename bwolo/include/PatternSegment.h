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

#ifndef PATTERNSEGMENT_H_
#define PATTERNSEGMENT_H_
#ifdef DEBUG_VERBOSE
#include <iostream>
#endif

#include <seqan/sequence.h>
#include "BwoloTypeConfig.h"

using namespace seqan;



inline unsigned int getPatternSegmentBeginPos(BwoloTPattern const& pattern, unsigned short segmentPosition, unsigned short errNumber){
	return (length(pattern)*segmentPosition)/(errNumber+2);
}

inline unsigned int getPatternSegmentEndPos(BwoloTPattern const& pattern, unsigned short segmentPosition, unsigned short errNumber){
	return (length(pattern)*(segmentPosition+1))/(errNumber+2);
}

/**
 * \fn Segment<BwoloTPattern> getPatternSegment(BwoloTPattern & pattern, unsigned int segmentPosition)
 * \brief retourne le segment à la position donnée du pattern. Pour l'instant, la découpe se fait en 5 parties egales
 * \param[in] pattern le pattern à découpé
 * \param[in] segmentPosition la position du segment à revoyer, entre 0 et 4 compris
 * \return le segment à la position donnée
 */
Segment<BwoloTPattern>
getPatternSegment(BwoloTPattern & pattern, unsigned short segmentPosition, unsigned short errNumber);

/**
 * \fn StringSet<BwoloTPattern> segmentPattern(BwoloTPattern & pattern)
 * \brief découpe le pattern donnée en 5 segments egaux et retourne l'ensemble sous la forme d'un StringSet
 * \param [in] pattern pattern à découper
 * \return le pattern découpé sous la forme d'un StringSet ou chaque sous chaine correspond à une découpe.
 */
StringSet<BwoloTPattern> segmentPattern(BwoloTPattern & pattern, unsigned short errNumber);


String<Iterator<BwoloTPattern, Rooted>::Type> getBeginSegmentPatternIterator(BwoloTPattern & pattern, unsigned short errNumber);
String<Iterator<BwoloTPattern, Rooted>::Type> getEndSegmentPatternIterator(BwoloTPattern & pattern, unsigned short errNumber);

#endif /* PATTERNSEGMENT_H_ */

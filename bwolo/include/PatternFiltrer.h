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

#ifndef PATTERNFILTRER_H_
#define PATTERNFILTRER_H_

#if (defined DEBUG_VERBOSE || defined DEBUG_WHILE_VERBOSE)
#include <iostream>
#endif

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "BwoloTypeConfig.h"
#include "candidate.h"
#include "CandidateLog.h"

using namespace seqan;

namespace PatternFiltrer{
/**
 * \fn void filtrateFromPatternSegmentPosition( String<String<candidate::Candidate> > & positionDependantResults, BwoloTIterator const & indexRootIt, Iterator<StringSet<BwoloTPattern> >::Type const & patternSegmentIt)
 * \brief filtrate all candidate from a segment of the pattern. Search with all size of seed from this segment.
 * \param[out] positionDependantResults all the return candidates. They are store in a String of String of Candidate structure. First dimension is the size, begin from 2 segment large to the max size allow by the segment position. The second dimension is the order in which they are found.
 * \param[in] indexRootIt iterator pointing the root of the FM-index where begin the search.
 * \param[in] patternSegmentIt iterator pointing the first segment in the partern of the different seeds use for filtering. The different seeds will begin by this segment and will be extend by 2, 3.. segment until the last segment of the pattern.
 */
	void filtrateFromPatternSegmentPosition(
			String<String<candidate::Candidate> > & positionDependantResults,
			BwoloTIterator const & indexRootIt,
			Iterator<StringSet<BwoloTPattern> >::Type const & patternSegmentIt);
/**
 * \fn void filtrate(String<String<String<candidate::Candidate> > > & results, BwoloTIterator const & indexRootIt, StringSet<BwoloTPattern> const & pattern)
 * \brief filtrate candidates following the tested algorithm
 * The filtrate algorithm. The pattern is cut in many segments (5). Then different aproximate seed is search.
 * This seed is create as : two exact segment at each extremity and a segments with exactly 1 error in the seed. The size of the seed varies from 2 exact segments to all segments. All seed are search.
 * A import work is done to prevent unnecessary search : insertion next to deletion between 2 contiguous nucleotides (in contiguous segments) or if two seeds started at the same segment but by different way of error are the same.
 * The results are return in a structure that store candidate firstly by their first segment position in the pattern then by their size (in segment) and finally their order of discovery (quite random).
 * \param[out] results return filtered candidates
 * \param[in] indexRootIt iterator pointing the root of the FM-index where begin the search.
 * \param[in] pattern pattern cut in segment.
 */
	void filtrate(
			String<String<String<candidate::Candidate> > > & results,
			BwoloTIterator const & indexRootIt,
			StringSet<BwoloTPattern> const & pattern);
}

#endif /* PATTERNFILTRER_H_ */

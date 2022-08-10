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

#ifndef SEGMENTERRBROWSE_H_
#define SEGMENTERRBROWSE_H_

#ifdef DEBUG_VERBOSE
#include <iostream>
#endif

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "BwoloTypeConfig.h"
#include "CandidateLog.h"
#include "candidate.h"

using namespace seqan;

namespace SegmentErrBrowse{

/**
 * \fn bool exactSegmentIterate(candidate::Candidate & candidate, BwoloTPattern const & patternSegment)
 * \brief search exactly the pattern segment from the candidate position. It return if the pattern segment exists. It iterate throw the candidate iterator, the candidate is modified.
 * \param[in, out] candidate candidate to extend with the pattern segment
 * \param[in] patternSegment pattern segment to search exactly
 */
bool exactSegmentIterate(
		candidate::Candidate & candidate,
		BwoloTPattern const & patternSegment);
/**
 * \fn void filtrateExactSegment( String<candidate::Candidate> & candidateList, candidate::Candidate actualPotentialCandidate, BwoloTPattern const & patternSegment)
 * \brief search exactly the pattern segment from the actual candidate in the index. The new candidate, if this one exists, is added in the String (use as vector).
 * \param[out] candidateList where adding the new candidate if exists.
 * \param[in] actualPotentialCandidate candidate to extend with the pattern segment
 * \param[in] patternSegment pattern segment to search exactly
 */
void filtrateExactSegment(
		String<candidate::Candidate> & candidateList,
		candidate::Candidate const actualPotentialCandidate,
		BwoloTPattern const & patternSegment
);

/**
 * \fn void filtrateSegmentWithError( String<candidate::Candidate> & candidateList, candidate::Candidate const & actualPotentialCandidate, BwoloTPattern const & patternSegment)
 * \brief search the pattern segment from the candidate position and add new candidate in the String used as an vector. It does the allowed operation (error) on the first nucleotide and all the operation on the other nucleotide of the pattern.
 * \param[out] candidateList new candidate found after ading the segment
 * \param[in] actualPotentialCandidate actual candidate used to find next candidate.
 * \patternSegment patternSegment to search
 */
void filtrateSegmentWithError(
		String<candidate::Candidate> & candidateList,
		candidate::Candidate const & actualPotentialCandidate,
		BwoloTPattern const & patternSegment
);

/**
 * \fn void filtrateFirstNcSegmentWithError( String<candidate::Candidate> & candidateList, candidate::Candidate const & actualPotentialCandidate, NcOp::NcOp operationAllowed, BwoloTPattern const & patternSegment)
 * \brief search the pattern segment from the candidate position and add new candidate in the String used as an vector. It only does the allowed operation (error) on the first nucleotide. It is used for adding the operation not ever done on an existing candidate if this one have ever being seen in the logger
 * \param[out] candidateList new candidate found after ading the segment
 * \param[in] actualPotentialCandidate actual candidate used to find next candidate.
 * \operationAllowed operation to do in the first nucleotide to find new candidate
 * \patternSegment patternSegment to search
 */
void filtrateFirstNcSegmentWithError(
		String<candidate::Candidate> & candidateList,
		candidate::Candidate const & actualPotentialCandidate,
		NcOp::NcOp operationAllowed,
		BwoloTPattern const & patternSegment
);

/**
 * \fn void browseErr(String<candidate::Candidate> & potentialCandidates, candidate::Candidate const & actualPotentialCandidate, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt, bool isSegLastNucleotide )
 * \brief search the pattern segment suffix by beginning by an error (insertion, deletion, substitution) in the actual branch of the index where the potential candidate is. It insert all the new potential candidate in the given String (used as a vector).
 * \param[out] potentialCandidates String used as a vector where the new potential candidates are adding.
 * \param[in] actualPotentialCandidate potential candidate from where the new pattern segment will be add with an error.
 * \param[in] patternSegmentIt pattern suffix to search with a substitution or a deletion on the first nucleotide or an insertion before this.
 * \param isSegLastNucleotide if the suffix is only 1 nucleotide length (last nucleotide of the segment). It is needed to give correct further allow operation in the first nucleotide of the next segment in the new potential candidate.
 */
void browseErr(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool isSegLastNucleotide
);

/**
 * \fn void browseErrWoDeletion(String<candidate::Candidate> & potentialCandidates, candidate::Candidate const & actualPotentialCandidate, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt, bool isSegLastNucleotide )
 * \brief search the pattern segment suffix by beginning by a substitution or an insertion in the actual branch of the index where the potential candidate is. It insert all the new potential candidate in the given String (used as a vector).
 * \param[out] potentialCandidates String used as a vector where the new potential candidates are adding.
 * \param[in] actualPotentialCandidate potential candidate from where the new pattern segment will be add with an error.
 * \param[in] patternSegmentIt pattern suffix to search with a substitution on the first nucleotide or an insertion before this.
 * \param isSegLastNucleotide if the suffix is only 1 nucleotide length (last nucleotide of the segment). It is needed to give correct further allow operation in the first nucleotide of the next segment in the new potential candidate.
 */
void browseErrWoDeletion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool isSegLastNucleotide
		);

/**
 * \fn void browseErrWoInsertion(String<candidate::Candidate> & potentialCandidates, candidate::Candidate const & actualPotentialCandidate, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt, bool isSegLastNucleotide )
 * \brief search the pattern segment suffix by beginning by a deletion or a substitution in the actual branch of the index where the potential candidate is. It insert all the new potential candidate in the given String (used as a vector).
 * \param[out] potentialCandidates String used as a vector where the new potential candidates are adding.
 * \param[in] actualPotentialCandidate potential candidate from where the new pattern segment will be add with an error.
 * \param[in] patternSegmentIt pattern suffix to search with a deletion or a substitution on the first nucleotide.
 * \param isSegLastNucleotide if the suffix is only 1 nucleotide length (last nucleotide of the segment). It is needed to give correct further allow operation in the first nucleotide of the next segment in the new potential candidate.
 */
void browseErrWoInsertion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool isSegLastNucleotide
		);
/**
 * \fn void browseLastInsertion( String<candidate::Candidate> & potentialCandidates, candidate::Candidate const & actualPotentialCandidate )
 * \brief search new potential candidate by adding a nucleotide in the actual branch of the index where the potential candidate is. It is used to do the last insertion cause other functions insert before going down.
 * \param[out] potentialCandidates String used as a vector where the new potential candidates are adding.
 * \param[in] actualPotentialCandidate potential candidate from where it will be adding an insertion.
 */
void browseLastInsertion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate
		);
/**
 * \fn void browseErrFirstInsertion( String<candidate::Candidate> & potentialCandidates, candidate::Candidate const & actualPotentialCandidate, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt, bool );
 * \brief search pattern segment suffix by beginning by a insertion in the actual branch of the index where the potential candidate is. It insert all the new potential candidate in the given String (used as a vector).
 * \param[out] potentialCandidates String used as a vector where the new potential candidates are adding.
 * \param[in] actualPotentialCandidate potential candidate from where the new pattern segment will be add with a insertion before the first nucleotide.
 * \param[in] patternSegmentIt pattern suffix to search with a insertion before the first nucleotide and add for having new potential candidate.
 */
void browseErrFirstInsertion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool
		);
/**
 * \fn void browseErrFirstDeletion( String<candidate::Candidate> & potentialCandidates,	candidate::Candidate const & actualPotentialCandidate, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt, bool isSegLastNucleotide );
 * \brief search pattern segment suffix by beginning by a deletion (of the first nucleotide in the pattern) in the actual branch of the index where the potential candidate is. It insert all the new potential candidate in the given String (used as a vector).
 * \param[out] potentialCandidates String used as a vector where the new potential candidates are adding.
 * \param[in] actualPotentialCandidate potential candidate from where the new pattern segment will be add with a deletion of the first nucleotide.
 * \param[in] patternSegmentIt pattern suffix to search with a deletion of the first nucleotide and add for having new potential candidate.
 * \param isSegLastNucleotide if the suffix is only 1 nucleotide length (last nucleotide of the segment). It is needed to give correct further allow operation in the first nucleotide of the next segment in the new potential candidate.
 */
void browseErrFirstDeletion(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool isSegLastNucleotide
		);

/**
 * \fn void browseErrFirstSubstitution( String<candidate::Candidate> & potentialCandidates, candidate::Candidate const & actualPotentialCandidate, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt, bool );
 * \brief search pattern segment suffix by beginning by a substitution in the actual branch of the index where the potential candidate is. It insert all the new potential candidate in the given String (used as a vector).
 * \param[out] potentialCandidates String used as a vector where the new potential candidates are adding.
 * \param[in] actualPotentialCandidate potential candidate from where the new pattern segment will be add with a substitution in the first nucleotide.
 * \param[in] patternSegmentIt pattern suffix to search with a substitution in the first nucleotide and add for having new potential candidate.
 */
void browseErrFirstSubstitution(
		String<candidate::Candidate> & potentialCandidates,
		candidate::Candidate const & actualPotentialCandidate,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		bool
		);

/**
 * \fn bool deletionErrIterate( BwoloTIterator & candidateIt, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt)
 * \brief search the pattern segment suffix throw the given index branch by begin by a deletion of the first nucleotide of the segment suffix. candidateIt is modified after this.
 * \return return true if it is possible to go down by begin by the deletion of the first letter of the segment suffix. False it it is impossible.
 * \param[in, out] candidateIt index iterator pointing the branch where begin the search
 * \param[in] patternSegmentIt segment suffix to search by begin by the deletion of the first letter
 */
bool deletionErrIterate(
		BwoloTIterator & candidateIt,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt);

/**
 * \fn bool insertionErrIterate( BwoloTIterator & candidateIt, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt, BwoloTAlphabet alphabetChar)
 * \brief search the pattern segment suffix throw the given index branch by begin by a insertion of the first nucleotide by the given nucleotide. candidateIt is modified after this.
 * \return return true if it is possible to go down by begin by the insertion of the first letter of the segment suffix of the pattern by the letter given. False it it is impossible.
 * \param[in, out] candidateIt index iterator pointing the branch where begin the search
 * \param[in] patternSegmentIt segment suffix to search after the insertion of the given letter
 * \param  alphabetChar letter to insert before the segment suffix in the search
 */
bool insertionErrIterate(
		BwoloTIterator & candidateIt,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		BwoloTPatternAlphabet alphabetChar);
/**
 * \fn bool substitutionErrIterate( BwoloTIterator & candidateIt, Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt, BwoloTAlphabet alphabetChar)
 * \brief search the pattern segment suffix throw the given index branch by begin by a substitution of the first nucleotide by the given nucleotide. Only substitution is done, if the first letter of the segment suffix and the given letter are same, nothing is done, the function return false. candidateIt is modified after this.
 * \return return true if it is possible to go down by begin by the substitution of the first letter of the segment suffix of the pattern by the letter given. False it it is impossible or if the given letter is the same than the first letter of the segment suffix.
 * \param[in, out] candidateIt index iterator pointing the branch where begin the search
 * \param[in] patternSegmentIt segment suffix to search after the substitution of the first letter by the given letter
 * \param  alphabetChar letter to substitute instead of the first letter of the segment suffix in the search
 */
bool substitutionErrIterate(
		BwoloTIterator & candidateIt,
		Iterator<BwoloTPattern, Rooted>::Type const & patternSegmentIt,
		BwoloTPatternAlphabet alphabetChar);
}




#endif /* SEGMENTERRBROWSE_H_ */

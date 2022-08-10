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

#ifndef CANDIDATELOG_H_
#define CANDIDATELOG_H_

#ifdef DEBUG_VERBOSE
#include <iostream>
#endif
#include <unordered_map>

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "BwoloTypeConfig.h"
#include "candidate.h"
#include "FirstNucleotideOperation.h"
#include "CandidateFingerPrint.h"
using namespace seqan;

namespace candidateLog {

/**
 * Hashmap (C++11) use in the log system of candidates already seen
 */
typedef std::unordered_map<CandidateFingerPrint, NcOp::NcOp> hashMapLog;

/*
 * \struct  CandidateLogger
 * \brief log system of candidates already seen
 */
typedef struct CadidateLogger_s {
	hashMapLog hashMap; /* !< hash map use to store candidates already seen*/
} CandidateLogger;

/**
 * \fn bool contains(NcOp::NcOp & remainingOperation ,CandidateLogger const & logger, candidate::Candidate const & candidate)
 * \brief return true if the candidate is already store. Return also the remaining operation on the first nucleotide that have not even done
 * \param[out] remainingOperation remaining operation that have not done for the first nucleotide.
 * \param[in] logger CandidateLogger where the previous candidates is stored.
 * \param[in] candidate we want to know if it has already been processed and if there are still some operations to do.
 */
bool contains(NcOp::NcOp & remainingOperation ,CandidateLogger const & logger, candidate::Candidate const & candidate);

/**
 * \fn bool add(CandidateLogger & logger, candidate::Candidate const & candidate);
 * \brief insert a new candidate in the log system only if the candidate don't exist in it.
 * \Return indicating whether the element was successfully inserted or not.
 * \param[in,out] logger log system where insert the new candidate.
 * \param[in] candidate new candidate to add in the logger.
 */
bool add(CandidateLogger & logger, candidate::Candidate const & candidate);

/**
 * \fn bool update(CandidateLogger & logger, candidate::Candidate const & candidate, NcOp::NcOp opToAdd)
 * \brief update the operation done on the candidate in the logger. If the candidate is not found in the logger, nothing is done.
 * \return if the candidate was found in the logger an correctly update.
 * \param[in, out] logger log system where the candidate need to be update.
 * \param[in] candidate candidate to update the operation done on the first nucleotide.
 * \param[in] opToAdd operation done on the first nucleotide to add for the candidate.
 */
bool update(CandidateLogger & logger, candidate::Candidate const & candidate, NcOp::NcOp opToAdd);

/**
 * \fn void clear(CandidateLogger & logger)
 * \brief Reset the logger. All candidate log are removed and destroyed.
 *
 */
void clear(CandidateLogger & logger);
}

#endif /* CANDIDATELOG_H_ */

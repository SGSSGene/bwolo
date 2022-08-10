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

#ifndef CANDIDATE_H_
#define CANDIDATE_H_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "BwoloTypeConfig.h"
#include "FirstNucleotideOperation.h"

using namespace seqan;

namespace candidate{
/**
 * \enum NucleotidePosition
 * \brief NucleotidePosition est la position de l'erreur dans un segment
 */
typedef enum{
	noneErrPos, /*!< no error in the segment */
	firstNucleotide,/*!<  the error is in the first nucleotide of the segment  */
	centralNucleotide,/*!< the error is not at the beginning or end of the segment   */
	lastNucleotide/*!<  the error is in the last nucleotide of the segment  */
} NucleotidePosition;

/*typedef enum{
	nonePreviousErr,
	insertionErr,
	deletionErr,
	substitutionErr
} PreviousNucleotideErr; */

/**
 * \struct Candidate
 * \brief Candidate is a structure represent a potential candidate in construction as the filtration algorithm define it
 */
typedef struct Candidate_s{
/*	NucleotidePosition previousNucleotideErrPos;
	PreviousNucleotideErr previousNucleotideErrType;*/
	NcOp::NcOp firstNcOpSeg; /*!< firstNcOpSeg is a flag that define which operation is authorized in the first nucleotide of the new segment*/
	BwoloTIterator it; /*!< the iterator from the root of the index represent the actual candidate*/
}Candidate;

}



#endif /* CANDIDATE_H_ */

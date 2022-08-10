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

#ifndef INTERFACEMAINFONCTION_H_
#define INTERFACEMAINFONCTION_H_

#include <seqan/sequence.h>

#include "BwoloTypeConfig.h"
#include "candidate.h"
#include "MySequence.h"
#include "PotentialCandidates.h"


void getPatternsFromfasta(char* fastaFile, StringSet<mySequence> & patterns, bool isVerbose);

void getPatternsFromCStr(char* patternCstr, StringSet<mySequence> & patterns, bool isVerbose);

void printPotentialOccurences(unsigned int id, mySequence & pattern, PotentialCandidates & occurences, bool isVerbose);

void printPotentialOccurencesMultiple(StringSet<mySequence> & patterns, String<PotentialCandidates> & occurencesByPattern, bool isVerbose);

void shearchAndPrintPotentialOccurences_PatternString(StringSet<mySequence> & patterns, char* indexFile, unsigned short errNumber, bool isVerbose);

void printPotentialOccurences_PatternSequences(char* patternSequence, char* indexFile, unsigned short errNumber, bool isVerbose);

void printPotentialOccurences_PatternFasta(char* patternFile, char* indexFile, unsigned short errNumber, bool isVerbose);

#endif /* INTERFACEMAINFONCTION_H_ */

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

#ifndef BWOLOTYPECONFIG_H_
#define BWOLOTYPECONFIG_H_

#include <functional>

#include <seqan/sequence.h>
#include <seqan/index.h>
using namespace seqan;
/*
 * Type d'Alphabet qui sera utilisé pour le texte
 */
typedef Dna BwoloTIndexAlphabet;

/*
 * Type d'Alphabet qui sera utilisé pour le Pattern
 */
typedef Dna BwoloTPatternAlphabet;

/*
 * Stratégie pour la chaine de caractères du Texte
 */
typedef Alloc<> BwoloTLongTxtSpec;

/*
 * Stratégie pour la chaine de caractères du pattern
 */
typedef Alloc<> BwoloTPatternSpec;

/*
 * Stratégie d'index à utiliser
 */
typedef FMIndex<> BwoloTIdxSpec;

/*
 * type iterateur à utiliser
 */
typedef TopDown<> BwoloTItSpec;




//Ne pas modifier

/*
 * Type de chaine pour le texte
 */
typedef String<BwoloTIndexAlphabet, BwoloTLongTxtSpec> BwoloTText;

/*
 * Type de chaine pour le Pattern
 */
typedef String<BwoloTPatternAlphabet, BwoloTPatternSpec> BwoloTPattern;

/*
 * Type d'Index utilisé
 */
typedef Index<BwoloTText, BwoloTIdxSpec> BwoloTIndex;

/*
 * Type d'Iterator
 */
typedef typename Iterator< BwoloTIndex, BwoloTItSpec >::Type BwoloTIterator;

typedef typename Size<BwoloTIterator>::Type BwoloTIteratorSize;

typedef typename Iterator< BwoloTPattern, Rooted >::Type BwoloTPatternIterator;

#endif /* BWOLOTYPECONFIG_H_ */

///*TODO : changer ces variables globales en paramêtres */
//extern unsigned short ERR_NUMBER;
//
//extern unsigned short PATTERN_SEGMENT_NUMBER;

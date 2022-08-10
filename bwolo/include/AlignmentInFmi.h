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

#ifndef ALIGNMENTINFMI_H_
#define ALIGNMENTINFMI_H_

#include "BwoloTypeConfig.h"
#include <unordered_map>

class AligmentPosition{
public:
	Iterator<BwoloTPattern, Rooted>::Type patternIt;
	unsigned short beginPatternSegment;
	BwoloTIterator it;
	//constructeur de copie
	AligmentPosition( AligmentPosition const& k);
	AligmentPosition(Iterator<BwoloTPattern, Rooted>::Type const& patternIt, unsigned short beginPatternSegment, BwoloTIterator const& it);

private :
	void init(Iterator<BwoloTPattern, Rooted>::Type const& patternIt, unsigned short beginPatternSegment, BwoloTIterator const& it);
};

namespace std {
//fonction hash BwoloTIterator
template<typename TIndex, typename TSpec> struct hash<Iter<TIndex, VSTree<TSpec> > > {
	std::size_t operator()(const Iter<TIndex, VSTree<TSpec> > & k) const {
		typedef typename Size<Iter<TIndex, VSTree<TSpec> > >::Type TSize;
		return std::hash<TSize>()(range(k).i1)
				^ (std::hash<TSize>()(range(k).i2)
						^ (std::hash<TSize>()(repLength(k)) << 1) << 1);
	}
};

//fonction hash TPattern Iterator
template<typename TContainer, typename TIterator, typename TSpec> struct hash<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > > {
	std::size_t operator()(const Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & k) const {
		typedef typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type TPosition;
		return std::hash<TPosition >()(position(k)); // = position(k) ?
	}
};

//fonction equal_to BwoloTIterator
// vrai si x et y sont dans le même index
template< typename TIndex, typename TSpec > struct equal_to<Iter< TIndex, VSTree<TopDown<TSpec> > > > {
	bool operator()(const Iter< TIndex, VSTree<TopDown<TSpec> > >& x,
			const Iter< TIndex, VSTree<TopDown<TSpec> > >& y) const {
		return (range(x).i1 == range(y).i1 && range(x).i2 == range(y).i2 && repLength(x) == repLength(y));
	}
};


//fonction de Hash AligmentPosition necessaire pour la hashmap
template<> struct hash<AligmentPosition> {
	std::size_t operator()(const AligmentPosition& k) const {
		return std::hash<BwoloTIterator>()(k.it)
				^ (std::hash<Iterator<BwoloTPattern, Rooted>::Type>()(k.patternIt)
					^ (std::hash<unsigned short>()(k.beginPatternSegment) << 1) <<1);
	}
};
//fonction equal_to AligmentPosition necessaire pour la hashmap
template<> struct equal_to<AligmentPosition> {
	bool operator()(const AligmentPosition& x,
			const AligmentPosition& y) const {
		return (equal_to<BwoloTIterator>()(x.it, y.it)
				&& x.beginPatternSegment == y.beginPatternSegment
				&& (x.patternIt == y.patternIt) );
	}
};
}

class ExtendingSeedState : public AligmentPosition{
public:
	unsigned short score;
	unsigned short initialSeedCoreSize;
	//constructeur de copie
	ExtendingSeedState(ExtendingSeedState const& k);
	ExtendingSeedState(Iterator<BwoloTPattern, Rooted>::Type const& patternIt,
			BwoloTIterator const&  it,
			unsigned short score,
			unsigned short beginPatternSegment,
			unsigned short initialSeedCoreSize);
private :
	void init(unsigned short score,
	unsigned short initialSeedCoreSize);
};

typedef std::unordered_map<AligmentPosition, ExtendingSeedState> AlignmentsMap;
typedef std::pair<AligmentPosition,ExtendingSeedState> AlignmentMapEntry;

AlignmentMapEntry newAlignmentMapEntry(ExtendingSeedState seed);

//insert l'AligmentMap la nouvelle entrée ou met a jour cette dernière avec un meilleur score si c'est le cas.
// return si l'entrée a été ajouté ou mise a jour.
bool AlignmentMapInsert(AlignmentsMap & map, AlignmentMapEntry entry);


void backtrakingNextStep(AlignmentsMap & finalAligments,
		AlignmentsMap & tPlusOneAlignments,
		AlignmentsMap &  tPlusTwoAlignments,
		ExtendingSeedState const & seedState,
		unsigned short maxScore);

void backtracking(AlignmentsMap & finalAligments,
		ExtendingSeedState const&  seed,
		unsigned short maxScore);

void printResult(AlignmentsMap & finalAligments);



#endif /* ALIGNMENTINFMI_H_ */

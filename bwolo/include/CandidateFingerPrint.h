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

#ifndef CANDIDATEFINGERPRINT_H_
#define CANDIDATEFINGERPRINT_H_

/**
 * \struct CandidateFingerPrint
 * \brief minimalist information to characterize a unique potential candidate in a suffix array.
 */
typedef struct CandidateFingerPrint_s {
	BwoloTIteratorSize rang1; /* !< inferior rank in the suffix array */
	BwoloTIteratorSize rang2;/* !< superior rank in the suffix array */
	BwoloTIteratorSize length;/* !< lenth of the candidate */
} CandidateFingerPrint;

namespace std {
//fonction de Hash necessaire pour la hashmap
template<> struct hash<CandidateFingerPrint> {
	std::size_t operator()(const CandidateFingerPrint& k) const {
		return std::hash<BwoloTIteratorSize>()(k.rang1)
				^ (std::hash<BwoloTIteratorSize>()(k.rang2)
						^ (std::hash<BwoloTIteratorSize>()(k.length) << 1) << 1);
	}
};
//fonction Equal necessaire pour la hashmap
template<> struct equal_to<CandidateFingerPrint> {
	bool operator()(const CandidateFingerPrint& lhs,
			const CandidateFingerPrint& rhs) const {
		return (lhs.rang1 == rhs.rang1 && lhs.rang2 == rhs.rang2
				&& lhs.length == rhs.length);
	}

};
}
#endif /* CANDIDATEFINGERPRINT_H_ */

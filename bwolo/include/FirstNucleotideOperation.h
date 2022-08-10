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

#ifndef FIRSTNUCLEOTIDEOPERATION_H_
#define FIRSTNUCLEOTIDEOPERATION_H_

namespace NcOp{
/**
 * \def INSERT_ALLOW 0000 0001 Flag mask for allowing insertion in the first nucleotide of a segment
 */
#define INSERT_ALLOW 0x1    // 0000 0001
/**
 * \def DELET_ALLOW 0000 0010 Flag mask for allowing deletion in the first nucleotide of a segment
 */
#define DELET_ALLOW 0x2     // 0000 0010
/**
 * \def SUBSTIT_ALLOW 0000 0100 Flag mask for allowing substitution in the first nucleotide of a segment
 */
#define SUBSTIT_ALLOW 0x4   // 0000 0100

//type prédéfini
/**
 * \def NO_NCOP_ALLOW 0000 0000 Preformatted flag word for allowing no operation in the first nucleotide of a segment
 */
#define NO_NCOP_ALLOW 0x0  // 0000 0111
/**
 * \def ALL_NCOP_ALLOW 0000 0111 Preformatted flag word for allowing all operation in the first nucleotide of a segment
 */
#define ALL_NCOP_ALLOW 0x7  // 0000 0111
/**
 * \def NO_INSERT 0000 0110 Preformatted flag word for allowing all operation except insertion in the first nucleotide of a segment
 */
#define NO_INSERT 0x6       // 0000 0110
/**
 * \def NO_DELET 0000 0101 Preformatted flag word for allowing all operation except deletion in the first nucleotide of a segment
 */
#define NO_DELET 0x5        // 0000 0101
/**
 * \def NO_SUBSTIT 0000 0011 Preformatted flag word for allowing all operation except substitution in the first nucleotide of a segment
 */
#define NO_SUBSTIT 0x3      // 0000 0011

//macro
/**
 * \def IS_INSERT_ALLOW(a) return if insertion is allow
 * \param a NcOp flag word
 */
#define IS_INSERT_ALLOW(a) ((a & INSERT_ALLOW) !=0 )
/**
 * \def IS_DELET_ALLOW(a) return if deletion is allow
 * \param a NcOp flag word
 */
#define IS_DELET_ALLOW(a) ((a & DELET_ALLOW) !=0 )
/**
 * \def IS_DELET_ALLOW(a) return if substitution is allow
 * \param a NcOp flag word
 */
#define IS_SUBSTIT_ALLOW(a) ((a & SUBSTIT_ALLOW) != 0)

/**
 * \def NcOp flag word use to indicate allowed operations in the first nucleotide of a segment
 */
typedef Byte NcOp;


}


#endif /* FIRSTNUCLEOTIDEOPERATION_H_ */

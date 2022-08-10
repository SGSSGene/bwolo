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

#include "../include/CandidateLog.h"

namespace candidateLog{

bool contains(NcOp::NcOp & remainingOperation ,CandidateLogger const & logger, candidate::Candidate const & candidate){
#ifdef DEBUG_VERBOSE
	std::cout <<"candidateLog::contains :"<< std::endl;
	std::cout << "\tremainingOperation : " << int(remainingOperation) << std::endl;
	std::cout << "\tcandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(candidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(candidate.firstNcOpSeg) << std::endl;
#endif
	CandidateFingerPrint key;
	remainingOperation = NO_NCOP_ALLOW;
	key.length = repLength(candidate.it);
	key.rang1 = range(candidate.it).i1;
	key.rang2 = range(candidate.it).i2;
	hashMapLog::const_iterator rep = logger.hashMap.find(key);
	if(rep == logger.hashMap.end()){
		remainingOperation = candidate.firstNcOpSeg;
		return false;
	}
	NcOp::NcOp OpDone = rep->second;
	remainingOperation = (~OpDone) & candidate.firstNcOpSeg; //operations restantes = celles qui n'ont pas encore été faite et qui sont a faire pour ce candidat
	return true;
}

bool add(CandidateLogger & logger, candidate::Candidate const & candidate){
#ifdef DEBUG_VERBOSE
	std::cout <<"candidateLog::add :"<< std::endl;
	std::cout << "\tcandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(candidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(candidate.firstNcOpSeg) << std::endl;
#endif
	CandidateFingerPrint key;
	key.length = repLength(candidate.it);
	key.rang1 = range(candidate.it).i1;
	key.rang2 = range(candidate.it).i2;
	std::pair<CandidateFingerPrint,NcOp::NcOp> newEntry(key, candidate.firstNcOpSeg);
	std::pair<hashMapLog::iterator,bool> rep= logger.hashMap.insert(newEntry); //je le rajoute
	return rep.second;
}

bool update(CandidateLogger & logger, candidate::Candidate const & candidate, NcOp::NcOp opToAdd){
#ifdef DEBUG_VERBOSE
	std::cout <<"candidateLog::update :"<< std::endl;
	std::cout << "\tcandidate : "<<std::endl;
	std::cout<<"\t\tit : " << representative(candidate.it) << std::endl;
	std::cout<<"\t\tfirstNcOpSeg : " << int(candidate.firstNcOpSeg) << std::endl;
	std::cout<<"\topToAdd"<<int(opToAdd)<<std::endl;
#endif
	CandidateFingerPrint key;
	key.length = repLength(candidate.it);
	key.rang1 = range(candidate.it).i1;
	key.rang2 = range(candidate.it).i2;
	hashMapLog::iterator rep = logger.hashMap.find(key);
	if(rep != logger.hashMap.end()){
		rep->second |= opToAdd; //je mets à jour les opérations effectués.
	}else{
		return false; //le candidat n'existe pas dans la table.
	}
	return true;
}

void clear(CandidateLogger & logger){
#ifdef DEBUG_VERBOSE
	std::cout <<"candidateLog::clear :"<< std::endl;
#endif
	logger.hashMap.clear();
}

}

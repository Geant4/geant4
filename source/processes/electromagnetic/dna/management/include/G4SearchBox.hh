// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, so do not hesitate to send us your feedback!
//
// In order for Geant4-DNA to be maintained and still open-source, article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we ask that you please cite the following papers reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifdef DNADEV
#ifndef G4SEARCHBOX_HH_
#define G4SEARCHBOX_HH_

#include "G4IT.hh"
#include <map>

class G4IT_Base;

class G4SearchBoxResult
{
public:
	typedef G4SearchBoxResult* Ptr;
	struct ResultState
	{
		G4IT_Base* fITBase;
		double fDistance;
	};
};

template<typename SearchT>
class G4SearchBox
{
public:
	G4SearchBox();
	virtual ~G4SearchBox();

	void Rebuild();

	template <typename PushedT> void Insert(const PushedT* toBeInserted)
	{
		fSearchMap[*toBeInserted]->Insert(toBeInserted);
	}

	template <typename LookInType, typename PositionT> G4SearchBoxResult::Ptr Nearest(const PositionT& pos)
	{
		return fSearchMap[LookInType()]->Nearest(pos);
	}

	template <typename LookInType, typename PositionT> G4SearchBoxResult::Ptr NearestInRange(const PositionT& pos, double range)
	{
		return fSearchMap[LookInType()]->NearestInRange(pos, range);
	}

protected:
	// TODO
	typedef std::map<G4IT_Base, SearchT*> SearchMap;
	typedef SearchT SearchType;

	SearchMap fSearchMap;
};

#endif /* G4SEARCHBOX_HH_ */
#endif

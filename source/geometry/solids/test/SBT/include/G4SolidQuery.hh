//
// G4SolidQuery.hh
//
// Declaration of the pure virtual class G4SolidQuery: a class
// that gives an inherited object a standard method in which
// to be queried for a (pointer to a) G4VSolid
//

#ifndef G4SolidQuery_hh
#define G4SolidQuery_hh

#include "G4VSolid.hh"

class G4SolidQuery {
	public:
	G4SolidQuery() {;}
	virtual ~G4SolidQuery() {;}
	
	virtual G4VSolid *GetSolid() const = 0;
};

#endif

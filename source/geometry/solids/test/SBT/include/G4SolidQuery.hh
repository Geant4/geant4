//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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

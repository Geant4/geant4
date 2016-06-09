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

#ifndef G4VFieldPropagation_h
#define  G4VFieldPropagation_h 1

#include "G4KineticTrackVector.hh"
#include "G4V3DNucleus.hh"

class G4VFieldPropagation
{
public:
  G4VFieldPropagation();
  virtual ~G4VFieldPropagation();

private:
  G4VFieldPropagation(const  G4VFieldPropagation &right);
  const G4VFieldPropagation & operator=(const G4VFieldPropagation & right);
  G4int operator==(const G4VFieldPropagation & right) const;
  G4int operator!=(const G4VFieldPropagation & right) const;

public:
  virtual void Init(G4V3DNucleus * theNucleus) = 0;
  virtual void Transport(G4KineticTrackVector &theActive,
			 const G4KineticTrackVector &theSpectators,
			 G4double theTimeStep) = 0;
  virtual G4ThreeVector GetMomentumTransfer() const =0;
};

#endif


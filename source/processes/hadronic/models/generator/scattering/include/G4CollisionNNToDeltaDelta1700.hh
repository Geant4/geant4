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

#ifndef G4CollisionNNToDeltaDelta1700_h
#define G4CollisionNNToDeltaDelta1700_h

#include "globals.hh"
#include "G4GeneralNNCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include "g4std/vector"

class G4CollisionNNToDeltaDelta1700 : public G4GeneralNNCollision
{

public:

  G4CollisionNNToDeltaDelta1700();

  virtual ~G4CollisionNNToDeltaDelta1700() {}


  virtual G4String GetName() const { return "NN -> Delta Delta(1700) Collision"; }
  virtual const G4std::vector<G4String>& GetListOfColliders(G4int whichOne) const
  {
    G4Exception("Tried to call G4CollisionNNToDeltaDelta1700::GetListOfColliders. Please find out why!");
    G4std::vector<G4String> * aList = new G4std::vector<G4String>;
    return *aList;
  } 
  

protected:
  
  virtual const G4CollisionVector* GetComponents() const { return components; } 

private:  

  G4CollisionVector* components;

};

#endif

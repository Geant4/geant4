//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, A. Feliciello, 30th June 1998
// -----------------------------------------------------------------------------

#include "G4CollisionInitialState.hh"
#include "G4SystemOfUnits.hh"
#include "G4BCAction.hh"
// Class G4CollisionInitialState

G4CollisionInitialState::G4CollisionInitialState() :
   theCollisionTime(DBL_MAX), thePrimary(0), theTarget(0), theFSGenerator(0)
{
}


G4CollisionInitialState::G4CollisionInitialState(G4double time,
   G4KineticTrack * aPrimary, G4KineticTrack * aTarget)
{
   theCollisionTime = time;
   thePrimary       = aPrimary;
   theTarget        = aTarget;
   theTs.clear();
   theFSGenerator = 0;
}

// +new interface post pion:
G4CollisionInitialState::G4CollisionInitialState(G4double time,
   G4KineticTrack * aPrimary, const G4KineticTrackVector & aTarget,
   G4BCAction * aFSGenerator)
{
   theCollisionTime = time;
   thePrimary       = aPrimary;
   theTarget = 0;
   for (size_t i=0; i<aTarget.size(); i++) theTs.push_back(aTarget[i]);
   theFSGenerator = aFSGenerator;
}
// -new interface post pion:


G4CollisionInitialState::G4CollisionInitialState(const G4CollisionInitialState & right)
{
     theCollisionTime = right.theCollisionTime;
     thePrimary       = right.thePrimary;
     theTarget        = right.theTarget;
     for (size_t i=0; i<right.theTs.size(); i++) theTs.push_back(right.theTs[i]);
     theFSGenerator   = right.theFSGenerator;
}

G4CollisionInitialState & G4CollisionInitialState::operator=(const G4CollisionInitialState& right)
{
   if (this != &right)
   {
      theCollisionTime = right.theCollisionTime;
      thePrimary       = right.thePrimary;
      theTarget        = right.theTarget;
      for (size_t i=0; i<right.theTs.size(); i++)
          theTs.push_back(right.theTs[i]);
      theFSGenerator   = right.theFSGenerator;
   }

   return *this;
}

//#include <typeinfo>

  G4KineticTrackVector * G4CollisionInitialState::
  GetFinalState()
  {
//    G4cerr << "what is the FS generator? "
//           << typeid(*theFSGenerator).name()
//	   <<G4endl;
    return theFSGenerator->GetFinalState(thePrimary, theTs);
  }

//#include <typeinfo>

void G4CollisionInitialState::Print() const
{
   G4int tgtPdg=theTarget ?
     theTarget->GetDefinition()->GetPDGEncoding() : 0;
   G4cout << "  collision " << this << " time: "
     << theCollisionTime/second << " proj: "
     << thePrimary << "/pdg=" << thePrimary->GetDefinition()->GetPDGEncoding()
     << " tgt: " << theTarget << "/pdg=" << tgtPdg
     << " Collision type: "<< typeid(*theFSGenerator).name();

}

  /*
std::ostream& G4CollisionInitialState::operator<<(std::ostream& out, const G4CollisionInitialState & collision)
{
  G4int tgtPdg=collision.GetTarget() ?
    collision.GetTarget()->GetDefinition()->GetPDGEncoding() : 0;
  out << "  collision " << collision << " time: "
    << collision.GetCollisionTime()/second << " proj: "
    << collision.GetPrimary() << "/pdg=" << collision.GetPrimary()->GetDefinition()->GetPDGEncoding()
    << " tgt: " << collision.GetTarget() << "/pdg=" << tgtPdg
    << " Collision type: "<< typeid(*collision.GetGenerator()).name();

  return out;
}
*/

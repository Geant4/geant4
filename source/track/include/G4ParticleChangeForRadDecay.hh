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
// $Id: G4ParticleChangeForRadDecay.hh 68795 2013-04-05 13:24:46Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// 
// ------------------------------------------------------------
//   Implemented for the new scheme                 25 Jan. 2000  H.Kurahige
//
// Class Description
//  This class is a concrete class for ParticleChange for RadDecay
//        
#ifndef G4ParticleChangeForRadDecay_h
#define G4ParticleChangeForRadDecay_h 1

#include "globals.hh"
#include "G4ios.hh"
class G4VTouchable;
#include "G4ParticleChangeForDecay.hh"

class G4ParticleChangeForRadDecay: public G4ParticleChangeForDecay
{ 
  public:
    // default constructor
    G4ParticleChangeForRadDecay(){}

    // destructor
    virtual ~G4ParticleChangeForRadDecay(){}

  protected:
    // hide copy constructor and assignment operaor as protected
    G4ParticleChangeForRadDecay(const G4ParticleChangeForRadDecay &) : G4ParticleChangeForDecay() {}
    G4ParticleChangeForRadDecay & operator=(const G4ParticleChangeForRadDecay &){return *this;}


  public: // with description
  void AddSecondary(G4Track* aSecondary);
    //  Add a secondary particle to theListOfSecondaries.
    // ------------------------------------------------------   


};

inline void G4ParticleChangeForRadDecay::AddSecondary(G4Track *aTrack)
{
  // add a secondary after size check
  if (theSizeOftheListOfSecondaries > theNumberOfSecondaries) {
    theListOfSecondaries->SetElement(theNumberOfSecondaries, aTrack);
    theNumberOfSecondaries++;
  } else {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VParticleChange::AddSecondary() Warning  ";
      G4cout << "theListOfSecondaries is full !! " << G4endl;
      G4cout << " The track is deleted " << G4endl;
    }
#endif
    G4Exception("G4ParticleChangeForRadDecay::AddSecondary",
                "TRACK101", JustWarning,
                "Secondary Bug is full. The track is deleted"); 
  }
}



#endif

















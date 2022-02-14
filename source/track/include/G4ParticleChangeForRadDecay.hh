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
// G4ParticleChangeForRadDecay
//
// Class description:
//
// Concrete class for ParticleChange for Radioactive Decay.

// Author: Hisaya Kurashige, 25 January 2000
// --------------------------------------------------------------------
#ifndef G4ParticleChangeForRadDecay_hh
#define G4ParticleChangeForRadDecay_hh 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleChangeForDecay.hh"

class G4VTouchable;

class G4ParticleChangeForRadDecay : public G4ParticleChangeForDecay
{
  public:

    G4ParticleChangeForRadDecay() {}
      // Default constructor

    virtual ~G4ParticleChangeForRadDecay() {}
      // Destructor

    inline void AddSecondary(G4Track* aSecondary);
      // Add a secondary particle to theListOfSecondaries

  protected:

    inline G4ParticleChangeForRadDecay(const G4ParticleChangeForRadDecay&);
    inline G4ParticleChangeForRadDecay& operator=(const G4ParticleChangeForRadDecay&);
      // Hidden copy constructor and assignment operator
};

// ----------------------
// Inline methods
// ----------------------

inline G4ParticleChangeForRadDecay::
G4ParticleChangeForRadDecay(const G4ParticleChangeForRadDecay&)
  : G4ParticleChangeForDecay() {}

inline
G4ParticleChangeForRadDecay&
G4ParticleChangeForRadDecay::operator=(const G4ParticleChangeForRadDecay&)
{
  return *this;
}

inline
void G4ParticleChangeForRadDecay::AddSecondary(G4Track* aTrack)
{
  // add a secondary after size check
  if(theSizeOftheListOfSecondaries > theNumberOfSecondaries)
  {
    theListOfSecondaries->SetElement(theNumberOfSecondaries, aTrack);
    ++theNumberOfSecondaries;
  }
  else
  {
    G4Exception("G4ParticleChangeForRadDecay::AddSecondary()", "TRACK101",
                JustWarning, "Secondaries buffer is full. Track is deleted!");
  }
}

#endif

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LFission.hh,v 1.2 1999-11-11 15:37:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Low-energy Model: Fission -- header file
// F.W. Jones, TRIUMF, 03-DEC-96
//  
// For further comments see G4LFission.cc.
//
// use -scheme for elastic scattering: HPW, 20th June 1997
// most of the code comes from the old Low-energy Fission class
//


#ifndef G4LFission_h
#define G4LFission_h 1
 
#include "g4rw/tphdict.h"
#include "globals.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4HadronicInteraction.hh"


class G4LFission : public G4HadronicInteraction
{
public:

   G4LFission();
   ~G4LFission();
 
   G4VParticleChange* ApplyYourself(const G4Track& aTrack,
                                    G4Nucleus& targetNucleus);

   static G4double Atomas(const G4double A, const G4double Z);

private:

   void init();

   G4double spneut[10];
};
#endif

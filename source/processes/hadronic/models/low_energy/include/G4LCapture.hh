// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LCapture.hh,v 1.2 1999-11-11 15:37:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Low energy model: neutron capture -- header file
// F.W. Jones, TRIUMF, 03-DEC-96
//  
// For further comments see G4LCapture.cc.
//
// use -scheme for elastic scattering: HPW, 20th June 1997
// most of the code comes from the old Low-energy Capture class
//


#ifndef G4LCapture_h
#define G4LCapture_h 1
 
#include "g4rw/tphdict.h"
#include "globals.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4LFission.hh"
#include "G4HadronicInteraction.hh"


class G4LCapture : public G4HadronicInteraction
{
public:

   G4LCapture();

   ~G4LCapture();
 
   G4VParticleChange* ApplyYourself(const G4Track& aTrack,
                                    G4Nucleus& targetNucleus);

private:

// Computes atomic mass in GeV using method from G4LFission
   inline
   G4double Atomas(const G4double A, const G4double Z)
   {
      return G4LFission::Atomas(A, Z)/GeV;
   }
};
#endif

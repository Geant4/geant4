// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuonMinusCaptureAtRest physics process --------
//                   by Vladimir Ivanchenko
//                     E-mail: Vladimir.Ivantchenko@cern.ch
//                            April 2000
// **************************************************************
//-----------------------------------------------------------------------------

#ifndef G4StopElementSelector_h
#define G4StopElementSelector_h 1
 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4ParticleDefinition.hh"
#include "g4std/iomanip"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4MuonMinus.hh"

class G4StopElementSelector
 
{ 
  private:
  // hide assignment operator as private 
      G4StopElementSelector& operator=(const G4StopElementSelector &right);
      G4StopElementSelector(const G4StopElementSelector& );
   
  public:
 
     G4StopElementSelector();
 
    ~G4StopElementSelector();

     G4Element* GetElement(const G4Material* aMaterial);
     G4double  GetMuonCaptureRate(G4double Z, G4double A);
     G4double  GetMuonDecayRate(G4double Z, G4double A);

  private:


};

#endif
 

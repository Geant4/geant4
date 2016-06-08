// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QGSMFragmentation.hh,v 1.2 1999/02/19 13:55:48 hpw Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#ifndef G4QGSMFragmentation_h
#define G4QGSMFragmentation_h 1

#include "G4VLongitudinalStringDecay.hh"

//******************************************************************************
class G4QGSMFragmentation:public G4VLongitudinalStringDecay
   {
public:
      G4QGSMFragmentation();
      ~G4QGSMFragmentation();

      const G4QGSMFragmentation & operator=(const G4QGSMFragmentation &right);
      int operator==(const G4QGSMFragmentation &right) const;
      int operator!=(const G4QGSMFragmentation &right) const;

  private:
   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py);      
   G4QGSMFragmentation(const G4QGSMFragmentation &right);

  private:
    // model parameters
    const G4double arho; 
    const G4double aphi;  
    const G4double an; 
    const G4double ala;  
    const G4double aksi; 
    const G4double alft;

  };

// Class G4QGSMFragmentation 
#endif



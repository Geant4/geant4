// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QGSMFragmentation.hh,v 1.3 1998/12/01 15:35:46 maxim Exp $
// GEANT4 tag $Name: geant4-00 $
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
      G4QGSMFragmentation(const G4QGSMFragmentation &right);
      ~G4QGSMFragmentation();

      const G4QGSMFragmentation & operator=(const G4QGSMFragmentation &right);
      int operator==(const G4QGSMFragmentation &right) const;
      int operator!=(const G4QGSMFragmentation &right) const;

  private:
   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py);      

  };
//******************************************************************************

// Class G4QGSMFragmentation 
#endif



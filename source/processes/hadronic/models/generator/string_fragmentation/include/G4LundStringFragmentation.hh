// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LundStringFragmentation.hh,v 1.1 1999/01/07 16:12:17 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $ Maxim Komogorov
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------

#ifndef G4LundStringFragmentation_h
#define G4LundStringFragmentation_h 1

#include "G4VLongitudinalStringDecay.hh"

//**************************************************************************************************************

class G4LundStringFragmentation: public G4VLongitudinalStringDecay
    {
public:
    G4LundStringFragmentation();
//    G4LundStringFragmentation(G4double sigmaPt);
    G4LundStringFragmentation(const G4LundStringFragmentation &right);
   ~G4LundStringFragmentation();

public:
    const G4LundStringFragmentation & operator=(const G4LundStringFragmentation &right);
    int operator==(const G4LundStringFragmentation &right) const;
    int operator!=(const G4LundStringFragmentation &right) const;


private:
   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py);      
};

//**************************************************************************************************************
// Class G4LundStringFragmentation 
#endif



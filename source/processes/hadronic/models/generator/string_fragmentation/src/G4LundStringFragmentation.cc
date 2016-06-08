// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LundStringFragmentation.cc,v 1.1.10.1 1999/12/07 20:51:56 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#include "G4LundStringFragmentation.hh"
#include "Randomize.hh"

// Class G4LundStringFragmentation 
//****************************************************************************************

G4LundStringFragmentation::G4LundStringFragmentation()
   {
   }

// G4LundStringFragmentation::G4LundStringFragmentation(G4double sigmaPt)
// : G4VLongitudinalStringDecay(sigmaPt)
//    {
//    }

G4LundStringFragmentation::G4LundStringFragmentation(const G4LundStringFragmentation &right)
   {
   }


G4LundStringFragmentation::~G4LundStringFragmentation()
   { 
   }

//****************************************************************************************

const G4LundStringFragmentation & G4LundStringFragmentation::operator=(const G4LundStringFragmentation &right)
   {
     G4Exception("G4LundStringFragmentation::operator= meant to not be accessable");
     return *this;
   }

int G4LundStringFragmentation::operator==(const G4LundStringFragmentation &right) const
   {
   return !memcmp(this, &right, sizeof(G4LundStringFragmentation));
   }

int G4LundStringFragmentation::operator!=(const G4LundStringFragmentation &right) const
   {
   return memcmp(this, &right, sizeof(G4LundStringFragmentation));
   }

//****************************************************************************************

G4double G4LundStringFragmentation::GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py)
    {
    const G4double  alund = 0.7/GeV/GeV; 

//    If blund get restored, you MUST adapt the calculation of zOfMaxyf.
//    const G4double  blund = 1;

    G4double z, yf;
    G4double Mass = pHadron->GetPDGMass();
    
    G4double Mt2 = Px*Px + Py*Py + Mass*Mass;
    G4double zOfMaxyf=alund*Mt2/(alund*Mt2 + 1.);
    G4double maxYf=(1-zOfMaxyf)/zOfMaxyf * exp(-alund*Mt2/zOfMaxyf);
    do
       {
       z = zmin + G4UniformRand()*(zmax-zmin);
//       yf = pow(1. - z, blund)/z*exp(-alund*Mt2/z);
	 yf = (1-z)/z * exp(-alund*Mt2/z);
       } 
    while (G4UniformRand()*maxYf > yf); 
    return z;
    }

//****************************************************************************************

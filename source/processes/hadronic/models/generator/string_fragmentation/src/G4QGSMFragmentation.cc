// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QGSMFragmentation.cc,v 1.3 1998/12/01 15:35:29 maxim Exp $
// GEANT4 tag $Name: geant4-00 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#include "G4QGSMFragmentation.hh"
#include "Randomize.hh"
#include "G4ios.hh"

// Class G4QGSMFragmentation 
//****************************************************************************************
 
G4QGSMFragmentation::G4QGSMFragmentation()
   {
   }

G4QGSMFragmentation::G4QGSMFragmentation(const G4QGSMFragmentation &right)
   {
   }

G4QGSMFragmentation::~G4QGSMFragmentation()
   {
   }

//****************************************************************************************

const G4QGSMFragmentation & G4QGSMFragmentation::operator=(const G4QGSMFragmentation &right)
   {
    G4Exception("G4QGSMFragmentation::operator= meant to not be accessable");
    return *this;
   }

int G4QGSMFragmentation::operator==(const G4QGSMFragmentation &right) const
   {
   return !memcmp(this, &right, sizeof(G4QGSMFragmentation));
   }

int G4QGSMFragmentation::operator!=(const G4QGSMFragmentation &right) const
   {
   return memcmp(this, &right, sizeof(G4QGSMFragmentation));
   }
 
//****************************************************************************************

G4double G4QGSMFragmentation::GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py)
    {
    const G4double arho = 0.5; 
    const G4double aphi = 0.;  
    const G4double an   = -0.5; 
    const G4double ala  = -0.75;  
    const G4double aksi = -1.; 
    const G4double alft = 0.5;
    
    G4double z;    
    G4double theA, d1, d2, yf;
    if (abs(PartonEncoding) < 10) 
        {
        switch(abs(PartonEncoding))
           {
        case 1:
	  theA = arho;
	  break;
	case 2:
	  theA = arho;
	  break;
	case 3:
	  theA = aphi;
	  break;
	default:
	  G4Exception("Unknown PDGencoding in G4QGSMFragmentation::G4LightConeZ");
          }
	do 	{
	    z  = zmin + G4UniformRand() * (zmax - zmin);
	    d1 =  (1. - z);
	    d2 =  (alft - theA);
	    yf = pow(d1, d2);
	    } 
	while (G4UniformRand() > yf);
	return z;
	}    
    switch(abs(PartonEncoding))
        {
    case 1103: 
        d2 =  (alft - (2.*an - arho));
        break;

    case 2101:   case 2103: 
        d2 = (alft - (2.*an - arho));
        break;

    case 3101:   case 3103: 
        d2 =  (alft - (2.*ala - arho));
        break;

    case 2203:
        d2 =  (alft - (2.*an  - arho));
        break;

    case 3201:   case 3203: 
        d2 =  (alft - (2.*ala - arho));
        break;

    case 3303: 
    default:    
        d2 =  (alft - (2.*aksi - arho));
        break;
        }
    do  {
        z = zmin + G4UniformRand() * (zmax - zmin);
        d1 =  (1. - z);
        yf = pow(d1, d2);
        } 
    while (G4UniformRand() > yf); 
    return z;
    }
    
//*********************************************************************************************

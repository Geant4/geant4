//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4LundStringFragmentation.cc,v 1.4 2005/06/04 13:47:01 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
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

G4LundStringFragmentation::G4LundStringFragmentation(const G4LundStringFragmentation &) : G4VLongitudinalStringDecay()
   {
   }


G4LundStringFragmentation::~G4LundStringFragmentation()
   { 
   }

//****************************************************************************************

const G4LundStringFragmentation & G4LundStringFragmentation::operator=(const G4LundStringFragmentation &)
   {
     throw G4HadronicException(__FILE__, __LINE__, "G4LundStringFragmentation::operator= meant to not be accessable");
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

G4double G4LundStringFragmentation::GetLightConeZ(G4double zmin, G4double zmax, G4int ,  G4ParticleDefinition* pHadron, G4double Px, G4double Py)
    {
    const G4double  alund = 0.7/GeV/GeV; 

//    If blund get restored, you MUST adapt the calculation of zOfMaxyf.
//    const G4double  blund = 1;

    G4double z, yf;
    G4double Mass = pHadron->GetPDGMass();
    
    G4double Mt2 = Px*Px + Py*Py + Mass*Mass;
    G4double zOfMaxyf=alund*Mt2/(alund*Mt2 + 1.);
    G4double maxYf=(1-zOfMaxyf)/zOfMaxyf * std::exp(-alund*Mt2/zOfMaxyf);
    do
       {
       z = zmin + G4UniformRand()*(zmax-zmin);
//       yf = std::pow(1. - z, blund)/z*std::exp(-alund*Mt2/z);
	 yf = (1-z)/z * std::exp(-alund*Mt2/z);
       } 
    while (G4UniformRand()*maxYf > yf); 
    return z;
    }

//****************************************************************************************

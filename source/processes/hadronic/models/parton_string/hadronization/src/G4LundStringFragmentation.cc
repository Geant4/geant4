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
//
// $Id: G4LundStringFragmentation.cc,v 1.5 2006/06/29 20:55:03 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
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

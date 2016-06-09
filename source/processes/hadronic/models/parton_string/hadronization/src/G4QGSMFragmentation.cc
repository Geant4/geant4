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
// $Id: G4QGSMFragmentation.cc,v 1.5 2006/06/29 20:55:05 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#include "G4QGSMFragmentation.hh"
#include "Randomize.hh"
#include "G4ios.hh"

// Class G4QGSMFragmentation 
//****************************************************************************************
 
G4QGSMFragmentation::G4QGSMFragmentation() :
arho(0.5), aphi(0.), an(-0.5), ala(-0.75), aksi(-1.), alft(0.5)
   {
   }

G4QGSMFragmentation::G4QGSMFragmentation(const G4QGSMFragmentation &) : G4VLongitudinalStringDecay(),
arho(0.5), aphi(0.), an(-0.5), ala(-0.75), aksi(-1.), alft(0.5)
   {
   }

G4QGSMFragmentation::~G4QGSMFragmentation()
   {
   }

//****************************************************************************************

const G4QGSMFragmentation & G4QGSMFragmentation::operator=(const G4QGSMFragmentation &)
   {
    throw G4HadronicException(__FILE__, __LINE__, "G4QGSMFragmentation::operator= meant to not be accessable");
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

G4double G4QGSMFragmentation::GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* , G4double , G4double )
{    
  G4double z;    
  G4double theA(0), d1, d2, yf;
  G4int absCode = std::abs( PartonEncoding );
  if (absCode < 10)
  { 
    if(absCode == 1 || absCode == 2) theA = arho;
    else if(absCode == 3) theA = aphi;
    else throw G4HadronicException(__FILE__, __LINE__, "Unknown PDGencoding in G4QGSMFragmentation::G4LightConeZ");

    do 	
    {
      z  = zmin + G4UniformRand() * (zmax - zmin);
      d1 =  (1. - z);
      d2 =  (alft - theA);
      yf = std::pow(d1, d2);
    } 
    while (G4UniformRand() > yf);
  }
  else
  {       
    if(absCode == 1103 || absCode == 2101 || 
       absCode == 2203 || absCode == 2103)
    {
      d2 =  (alft - (2.*an - arho));
    }
    else if(absCode == 3101 || absCode == 3103 ||
            absCode == 3201 || absCode == 3203)
    {
      d2 =  (alft - (2.*ala - arho));
    }
    else
    {
      d2 =  (alft - (2.*aksi - arho));
    }

    do  
    {
      z = zmin + G4UniformRand() * (zmax - zmin);
      d1 =  (1. - z);
      yf = std::pow(d1, d2);
    } 
    while (G4UniformRand() > yf); 
  }
  return z;
}
    
//*********************************************************************************************

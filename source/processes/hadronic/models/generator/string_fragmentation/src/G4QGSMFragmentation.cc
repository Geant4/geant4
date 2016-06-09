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
// $Id: G4QGSMFragmentation.cc,v 1.8 2002/12/12 19:17:56 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
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

G4QGSMFragmentation::G4QGSMFragmentation(const G4QGSMFragmentation &right) :
arho(0.5), aphi(0.), an(-0.5), ala(-0.75), aksi(-1.), alft(0.5)
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
  G4double z;    
  G4double theA(0), d1, d2, yf;
  G4int absCode = abs( PartonEncoding );
  if (absCode < 10)
  { 
    if(absCode == 1 || absCode == 2) theA = arho;
    else if(absCode == 3) theA = aphi;
    else G4Exception("Unknown PDGencoding in G4QGSMFragmentation::G4LightConeZ");

    do 	
    {
      z  = zmin + G4UniformRand() * (zmax - zmin);
      d1 =  (1. - z);
      d2 =  (alft - theA);
      yf = pow(d1, d2);
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
      yf = pow(d1, d2);
    } 
    while (G4UniformRand() > yf); 
  }
  return z;
}
    
//*********************************************************************************************

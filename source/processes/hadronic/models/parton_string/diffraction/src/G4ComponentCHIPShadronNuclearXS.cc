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
//  Calculation of the total, elastic and inelastic cross-sections
//  of hadron (proton, neutron, pi+, pi-, K+, K-, anti_proton, anti_neutron
//  interactions with nuclei based on CHIPS model
//
//   Created by V. Uzhinsky, 31.05.2011  

#include "G4ComponentCHIPShadronNuclearXS.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4VQCrossSection.hh"

#include "G4QProtonElasticCrossSection.hh"
#include "G4QProtonNuclearCrossSection.hh"

#include "G4QNeutronElasticCrossSection.hh"
#include "G4QNeutronNuclearCrossSection.hh"

#include "G4QAntiBaryonElasticCrossSection.hh" 
#include "G4QAntiBaryonNuclearCrossSection.hh"

#include "G4QPionMinusElasticCrossSection.hh" 
#include "G4QPionMinusNuclearCrossSection.hh"

#include "G4QPionPlusElasticCrossSection.hh" 
#include "G4QPionPlusNuclearCrossSection.hh"

#include "G4QKaonMinusElasticCrossSection.hh"
#include "G4QKaonMinusNuclearCrossSection.hh"

#include "G4QKaonPlusElasticCrossSection.hh" 
#include "G4QKaonPlusNuclearCrossSection.hh"
///////////////////////////////////////////////////////////////////////////////


G4ComponentCHIPShadronNuclearXS::G4ComponentCHIPShadronNuclearXS() 
: fUpperLimit( 10000 * GeV ),
  fLowerLimit( 10 * MeV )
{
  PxsManagerEl      = G4QProtonElasticCrossSection::GetPointer();
  PxsManagerInEl    = G4QProtonNuclearCrossSection::GetPointer();

  NxsManagerEl      = G4QNeutronElasticCrossSection::GetPointer();
  NxsManagerInEl    = G4QNeutronNuclearCrossSection::GetPointer();

  PBARxsManagerEl   = G4QAntiBaryonElasticCrossSection::GetPointer();
  PBARxsManagerInEl = G4QAntiBaryonNuclearCrossSection::GetPointer();

  PIPxsManagerEl    = G4QPionPlusElasticCrossSection::GetPointer();
  PIPxsManagerInEl  = G4QPionPlusNuclearCrossSection::GetPointer(); 

  PIMxsManagerEl    = G4QPionMinusElasticCrossSection::GetPointer();
  PIMxsManagerInEl  = G4QPionMinusNuclearCrossSection::GetPointer();

  KPxsManagerEl     = G4QKaonPlusElasticCrossSection::GetPointer();
  KPxsManagerInEl   = G4QKaonPlusNuclearCrossSection::GetPointer();

  KMxsManagerEl     = G4QKaonMinusElasticCrossSection::GetPointer();
  KMxsManagerInEl   = G4QKaonMinusNuclearCrossSection::GetPointer();

}

///////////////////////////////////////////////////////////////////////////////////////
G4ComponentCHIPShadronNuclearXS::~G4ComponentCHIPShadronNuclearXS()
{
}

////////////////////////////////////////////////////////////////////////////////
G4double G4ComponentCHIPShadronNuclearXS::GetTotalElementCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4double N)
{
  G4double momentum = std::sqrt(kinEnergy*(kinEnergy+2.*aParticle->GetPDGMass()));
  G4int PDGcode=aParticle->GetPDGEncoding();

  G4VQCrossSection* CHIPSmanagerEl=0; 
  G4VQCrossSection* CHIPSmanagerInEl=0;

  if     (PDGcode == 2212)   // Projectile is Proton
  {
   CHIPSmanagerEl=PxsManagerEl;  CHIPSmanagerInEl=PxsManagerInEl;
  } else if(PDGcode == 2112)  // Projectile is Neutron 
  {
   CHIPSmanagerEl=NxsManagerEl;  CHIPSmanagerInEl=NxsManagerInEl;   
  } else if(PDGcode == -2212) // Projectile is Anti-Proton
  {
   CHIPSmanagerEl=PBARxsManagerEl; CHIPSmanagerInEl=PBARxsManagerInEl;
  } else if(PDGcode == -2112) // Projectile is Anti-Neutron
  {
   CHIPSmanagerEl=PBARxsManagerEl; CHIPSmanagerInEl=PBARxsManagerInEl;
  } else if(PDGcode ==   211) // Projectile is Pi+ 
  {
   CHIPSmanagerEl=PIPxsManagerEl; CHIPSmanagerInEl=PIPxsManagerInEl;
  } else if(PDGcode ==  -211) // Projectile is Pi-
  {
   CHIPSmanagerEl=PIMxsManagerEl; CHIPSmanagerInEl=PIMxsManagerInEl;
  } else if(PDGcode ==  321)  // Projectile is K+ 
  {
   CHIPSmanagerEl=KPxsManagerEl; CHIPSmanagerInEl=KPxsManagerInEl; 
  } else if(PDGcode ==  -321) // Projectile is K- 
  {
   CHIPSmanagerEl=KMxsManagerEl; CHIPSmanagerInEl=KMxsManagerInEl; 
  } 

  G4double Xelastic(0.), Xinelastic(0.);

  if((CHIPSmanagerEl != 0) && (CHIPSmanagerInEl != 0))
  { 
   Xelastic   = CHIPSmanagerEl->GetCrossSection(false,momentum,Z,(G4int)N,PDGcode);
   Xinelastic = CHIPSmanagerInEl->GetCrossSection(false,momentum,Z,(G4int)N,PDGcode);
  }

  return Xelastic+Xinelastic; 
}

////////////////////////////////////////////////////////////////////////////////
G4double G4ComponentCHIPShadronNuclearXS::GetTotalIsotopeCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4int A )
{ return GetTotalElementCrossSection(aParticle, kinEnergy, Z, (G4double) A);  }

////////////////////////////////////////////////////////////////////////////////
G4double G4ComponentCHIPShadronNuclearXS::GetInelasticElementCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4double N)
{
  G4double momentum = std::sqrt(kinEnergy*(kinEnergy+2.*aParticle->GetPDGMass()));
  G4int PDGcode=aParticle->GetPDGEncoding();

  G4VQCrossSection* CHIPSmanagerInEl=0;

  if     (PDGcode == 2212)   // Projectile is Proton
  {
   CHIPSmanagerInEl=PxsManagerInEl;
  } else if(PDGcode == 2112)  // Projectile is Neutron 
  {
   CHIPSmanagerInEl=NxsManagerInEl;   
  } else if(PDGcode == -2212) // Projectile is Anti-Proton
  {
   CHIPSmanagerInEl=PBARxsManagerInEl;
  } else if(PDGcode == -2112) // Projectile is Anti-Neutron
  {
   CHIPSmanagerInEl=PBARxsManagerInEl;
  } else if(PDGcode ==   211) // Projectile is Pi+ 
  {
   CHIPSmanagerInEl=PIPxsManagerInEl;
  } else if(PDGcode ==  -211) // Projectile is Pi-
  {
   CHIPSmanagerInEl=PIMxsManagerInEl;
  } else if(PDGcode ==  321)  // Projectile is K+ 
  {
   CHIPSmanagerInEl=KPxsManagerInEl; 
  } else if(PDGcode ==  -321) // Projectile is K- 
  {
   CHIPSmanagerInEl=KMxsManagerInEl; 
  } 

  G4double Xinelastic(0.);

  if(CHIPSmanagerInEl != 0)
  { 
   Xinelastic = CHIPSmanagerInEl->GetCrossSection(false,momentum,Z,(G4int)N,PDGcode);
  }
  return Xinelastic; 
}

///////////////////////////////////////////////////////////////////////////////
G4double G4ComponentCHIPShadronNuclearXS::GetInelasticIsotopeCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4int A)
{return GetInelasticElementCrossSection(aParticle, kinEnergy, Z, (G4double) A); }
 
///////////////////////////////////////////////////////////////////////////////
G4double G4ComponentCHIPShadronNuclearXS::GetElasticElementCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4double N)
{
  G4double momentum = std::sqrt(kinEnergy*(kinEnergy+2.*aParticle->GetPDGMass()));
  G4int PDGcode=aParticle->GetPDGEncoding();

  G4VQCrossSection* CHIPSmanagerEl=0; 

  if     (PDGcode == 2212)   // Projectile is Proton
  {
   CHIPSmanagerEl=PxsManagerEl;
  } else if(PDGcode == 2112)  // Projectile is Neutron 
  {
   CHIPSmanagerEl=NxsManagerEl;   
  } else if(PDGcode == -2212) // Projectile is Anti-Proton
  {
   CHIPSmanagerEl=PBARxsManagerEl;
  } else if(PDGcode == -2112) // Projectile is Anti-Neutron
  {
   CHIPSmanagerEl=PBARxsManagerEl;
  } else if(PDGcode ==   211) // Projectile is Pi+ 
  {
   CHIPSmanagerEl=PIPxsManagerEl;
  } else if(PDGcode ==  -211) // Projectile is Pi-
  {
   CHIPSmanagerEl=PIMxsManagerEl;
  } else if(PDGcode ==  321)  // Projectile is K+ 
  {
   CHIPSmanagerEl=KPxsManagerEl;
  } else if(PDGcode ==  -321) // Projectile is K- 
  {
   CHIPSmanagerEl=KMxsManagerEl;
  } 

  G4double Xelastic(0.);

  if(CHIPSmanagerEl != 0)
  { 
   Xelastic   = CHIPSmanagerEl->GetCrossSection(false,momentum,Z,(G4int)N,PDGcode);
  }
  return Xelastic; 
}
 
///////////////////////////////////////////////////////////////////////////////
G4double G4ComponentCHIPShadronNuclearXS::GetElasticIsotopeCrossSection
(const G4ParticleDefinition* aParticle, G4double kinEnergy, G4int Z, G4int A)
{ return GetElasticElementCrossSection(aParticle, kinEnergy, Z, (G4double) A); }

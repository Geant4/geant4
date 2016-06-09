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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme).                     *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4GeneralSpaceNNCrossSection.cc
//
// Version:		B.1
// Date:		15/04/04
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4GeneralSpaceNNCrossSection.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"

#include <iomanip>


G4GeneralSpaceNNCrossSection::G4GeneralSpaceNNCrossSection ()
 : G4VCrossSectionDataSet("General Space NN")
{
  protonInelastic = new G4ProtonInelasticCrossSection();
  ionProton       = new G4IonProtonCrossSection();
  TripathiGeneral = new G4TripathiCrossSection();
  TripathiLight   = new G4TripathiLightCrossSection();
  Shen            = new G4IonsShenCrossSection();
  
  return;
}


G4GeneralSpaceNNCrossSection::~G4GeneralSpaceNNCrossSection ()
{
  delete protonInelastic;
  delete ionProton;
  delete TripathiGeneral;
  delete TripathiLight;
  delete Shen;
}
///////////////////////////////////////////////////////////////////////////////
//

G4bool G4GeneralSpaceNNCrossSection::IsIsoApplicable
 (const G4DynamicParticle* theProjectile, G4int ZZ, G4int AA)
{
  G4bool result = protonInelastic->IsIsoApplicable(theProjectile, ZZ, AA);
  if (!result)
  {
    result = ionProton->IsIsoApplicable(theProjectile, ZZ, AA);
    if (!result)
    {
      result = TripathiGeneral->IsIsoApplicable(theProjectile, ZZ, AA);
      if (!result)
        result = Shen->IsIsoApplicable(theProjectile, ZZ, AA);
    }
  }
  return result;
}

G4bool G4GeneralSpaceNNCrossSection::IsApplicable
  (const G4DynamicParticle* theProjectile, const G4Element* theTarget)
{
  G4bool result = protonInelastic->IsApplicable(theProjectile, theTarget);
  if (!result)
  {
    result = ionProton->IsApplicable(theProjectile, theTarget);
    if (!result)
    {
      result = TripathiGeneral->IsApplicable(theProjectile, theTarget);
      if (!result)
        result = Shen->IsApplicable(theProjectile, theTarget);
    }
  }
  return result;
}
///////////////////////////////////////////////////////////////////////////////
//

G4double G4GeneralSpaceNNCrossSection::GetCrossSection
  (const G4DynamicParticle* theProjectile, const G4Element* theTarget,
  G4double theTemperature)
{
  G4int nIso = theTarget->GetNumberOfIsotopes();
  G4double xsection = 0;
     
  if (nIso) {
    G4double sig;
    G4IsotopeVector* isoVector = theTarget->GetIsotopeVector();
    G4double* abundVector = theTarget->GetRelativeAbundanceVector();
    G4int ZZ;
    G4int AA;
     
    for (G4int i = 0; i < nIso; i++) {
      ZZ = (*isoVector)[i]->GetZ();
      AA = (*isoVector)[i]->GetN();
      sig = GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      xsection += sig*abundVector[i];
    }
   
  } else {
    G4int ZZ = G4lrint(theTarget->GetZ());
    G4int AA = G4lrint(theTarget->GetN());
    xsection = GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
  }
    
  return xsection;
}


G4double
G4GeneralSpaceNNCrossSection::GetZandACrossSection(const G4DynamicParticle* theProjectile,
                                                   G4int ZZ, G4int AA, G4double theTemperature)
{
  G4double result = 0.0;

  const G4double AT = AA;
  const G4double ZT = ZZ;
  const G4double AP = theProjectile->GetDefinition()->GetBaryonNumber();
  const G4double ZP = theProjectile->GetDefinition()->GetPDGCharge();

  if (verboseLevel >= 2)
  {
    G4cout <<"In G4GeneralSpaceNNCrossSection::GetCrossSection" <<G4endl;
    G4cout <<"Projectile A = " <<std::setw(8) <<AP 
           <<" Z = "           <<std::setw(8) <<ZP
           <<" Energy = "      <<theProjectile->GetKineticEnergy()/AP
           <<" MeV/nuc" <<G4endl;
    G4cout <<"Target     A = " <<std::setw(8) <<AT
           <<" Z = "           <<std::setw(8) <<ZT
           <<G4endl;
  }
  if (theProjectile->GetDefinition()==G4Proton::Proton())
  {
    if (ZT>5)
    {
      result = protonInelastic->
        GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4ProtonInelasticCrossSection" <<G4endl;
    }
    else
    {
      result = TripathiLight->
        GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
  }
  else if (AT==1 && ZT==1)
  {
    if (ZP>5)
    {
      result = ionProton->
        GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4IonProtonCrossSection" <<G4endl;
    }
    else
    {
      result = TripathiLight->
        GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
  }
  else
  {
    if (TripathiLight->IsIsoApplicable(theProjectile, ZZ, AA))
    {
      result = TripathiLight->
        GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
    else if (TripathiGeneral->IsIsoApplicable(theProjectile, ZZ, AA))
    {
      result = TripathiGeneral->
        GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiCrossSection" <<G4endl;
    }
    else if (Shen->IsIsoApplicable(theProjectile, ZZ, AA))
    {
      result = Shen->
        GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4IonsShenCrossSection" <<G4endl;
    }
  }
  if (verboseLevel >= 2)
  {
    G4cout <<"Cross-section = " <<result/millibarn <<" mbarn" <<G4endl;
    G4cout <<G4endl;
  }

  return result;
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
// * This  code  implementation  is  the  intellectual   property  of *
// * QinetiQ.   Rights  to use  this  software  are  granted  to  the *
// * European Space Agency  under  the  provisions  of  Clause 39  of *
// * the    Standard    Conditions   for  Contract.   ESA    contract *
// * 17191/03/NL/LvH  provides  the  rights to distribute source code *
// * through  the  Geant4 Collaboration Agreement for general  public *
// * use.  Some elements of this source code may be derived from code *
// * developed   by   other   member   institutes   of   the   Geant4 *
// * Collaboration, and the provision of  their code under the Geant4 *
// * Collaboration is acknowledged.                                   *
// *                                                                  *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4GeneralSpaceNNCrossSection.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include <iomanip>
////////////////////////////////////////////////////////////////////////////////
//
G4GeneralSpaceNNCrossSection::G4GeneralSpaceNNCrossSection ()
{
  protonInelastic = new G4ProtonInelasticCrossSection();
  ionProton       = new G4IonProtonCrossSection();
  TripathiGeneral = new G4TripathiCrossSection();
  TripathiLight   = new G4TripathiLightCrossSection();
  Shen            = new G4IonsShenCrossSection();
  
  return;
}
////////////////////////////////////////////////////////////////////////////////
//
G4GeneralSpaceNNCrossSection::~G4GeneralSpaceNNCrossSection ()
{
  delete protonInelastic;
  delete ionProton;
  delete TripathiGeneral;
  delete TripathiLight;
  delete Shen;
}
////////////////////////////////////////////////////////////////////////////////
//
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
////////////////////////////////////////////////////////////////////////////////
//
G4double G4GeneralSpaceNNCrossSection::GetCrossSection
  (const G4DynamicParticle* theProjectile, const G4Element* theTarget,
  G4double theTemperature)
{
  G4double result = 0.0;

  const G4double AT = theTarget->GetN();
  const G4double ZT = theTarget->GetZ();
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
        GetCrossSection(theProjectile, theTarget, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4ProtonInelasticCrossSection" <<G4endl;
    }
    else
    {
      result = TripathiLight->
        GetCrossSection(theProjectile, theTarget, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
  }
  else if (AT==1 && ZT==1)
  {
    if (ZP>5)
    {
      result = ionProton->
        GetCrossSection(theProjectile, theTarget, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4IonProtonCrossSection" <<G4endl;
    }
    else
    {
      result = TripathiLight->
        GetCrossSection(theProjectile, theTarget, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
  }
  else
  {
    if (TripathiLight->IsApplicable(theProjectile, theTarget))
    {
      result = TripathiLight->
        GetCrossSection(theProjectile, theTarget, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
    else if (TripathiGeneral->IsApplicable(theProjectile, theTarget))
    {
      result = TripathiGeneral->
        GetCrossSection(theProjectile, theTarget, theTemperature);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiCrossSection" <<G4endl;
    }
    else if (Shen->IsApplicable(theProjectile, theTarget))
    {
      result = Shen->
        GetCrossSection(theProjectile, theTarget, theTemperature);
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
////////////////////////////////////////////////////////////////////////////////
//

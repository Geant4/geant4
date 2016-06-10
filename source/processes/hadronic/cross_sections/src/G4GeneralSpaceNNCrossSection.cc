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
// 19 Aug 2011, V.Ivanchenko move to new design and make x-section per element
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include <iomanip>

#include "G4GeneralSpaceNNCrossSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4IonProtonCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadTmpUtil.hh"
#include "G4Proton.hh"

G4GeneralSpaceNNCrossSection::G4GeneralSpaceNNCrossSection ()
 : G4VCrossSectionDataSet("General Space NN")
{
  protonInelastic = new G4ProtonInelasticCrossSection();
  ionProton       = new G4IonProtonCrossSection();
  TripathiGeneral = new G4TripathiCrossSection();
  TripathiLight   = new G4TripathiLightCrossSection();
  Shen            = new G4IonsShenCrossSection();
  theProton       = G4Proton::Proton();  
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

G4bool G4GeneralSpaceNNCrossSection::IsElementApplicable
  (const G4DynamicParticle* theProjectile, G4int, const G4Material*)
{
  return (1 <= theProjectile->GetDefinition()->GetBaryonNumber());
}

///////////////////////////////////////////////////////////////////////////////
//

G4double G4GeneralSpaceNNCrossSection::GetElementCrossSection
  (const G4DynamicParticle* theProjectile, G4int ZT, const G4Material* mat)
{
  G4double result = 0.0;
  G4int ZP = G4lrint(theProjectile->GetDefinition()->GetPDGCharge()/eplus);

  if (verboseLevel >= 2)
  {
    G4int AP = theProjectile->GetDefinition()->GetBaryonNumber();
    G4cout <<"In G4GeneralSpaceNNCrossSection::GetCrossSection" <<G4endl;
    G4cout <<"Projectile A = " <<std::setw(8) <<AP 
           <<" Z = "           <<std::setw(8) <<ZP
           <<" Energy = "      <<theProjectile->GetKineticEnergy()/AP
           <<" MeV/nuc" <<G4endl;
    G4cout <<"Target     Z = " <<std::setw(8) <<ZT
           <<G4endl;
  }
  if (theProjectile->GetDefinition()==theProton)
  {
    if (ZT>5)
    {
      result = protonInelastic->
        GetElementCrossSection(theProjectile, ZT, mat);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4ProtonInelasticCrossSection" <<G4endl;
    }
    else
    {
      result = TripathiLight->
        GetElementCrossSection(theProjectile, ZT, mat);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
  }
  else if (ZT==1)
  {
    if (ZP>5)
    {
      result = ionProton->
        GetElementCrossSection(theProjectile, ZT, mat);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4IonProtonCrossSection" <<G4endl;
    }
    else
    {
      result = TripathiLight->
        GetElementCrossSection(theProjectile, ZT, mat);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
  }
  else
  {
    if (TripathiLight->IsElementApplicable(theProjectile, ZT, mat))
    {
      result = TripathiLight->
        GetElementCrossSection(theProjectile, ZT, mat);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiLightCrossSection" <<G4endl;
    }
    else if (TripathiGeneral->IsElementApplicable(theProjectile, ZT, mat))
    {
      result = TripathiGeneral->
        GetElementCrossSection(theProjectile, ZT, mat);
      if (verboseLevel >= 2)
        G4cout <<"Selecting G4TripathiCrossSection" <<G4endl;
    }
    else if (Shen->IsElementApplicable(theProjectile, ZT, mat))
    {
      result = Shen->
        GetElementCrossSection(theProjectile, ZT, mat);
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

///////////////////////////////////////////////////////////////////////////////
//

void G4GeneralSpaceNNCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
   outFile << "G4GeneralSpaceNNCrossSection calculates hadronic inelastic\n"
           << "cross sections of interest in space science, by using the\n"
           << "following cross sections:\n"
           << "- G4ProtonInelasticCrossSection : for proton projectile\n"
           << "  on targets with Z > 5;\n"
           << "- G4TripathiLightCrossSection : for proton projectile\n"
           << "  on targets with Z <= 5;\n"
           << "  for targets with Z = 1 and projectile Z <= 5;\n"
           << "  for neutron, or deuteron, or 3He, or alpha projectile\n"
           << "  with kinetic energy less than 10 GeV per nucleon,\n"
           << "  in any target;\n"
           << "  for 3He and 4He targets, for any projectile with\n"
           << "  kinetic energy less than 10 GeV per nucleon;\n"
           << "- G4IonProtonCrossSection : for projectile with Z > 5\n"
           << "  on hydrogen target;\n"
           << "- G4TripathiCrossSection : for any projectile with A >=3\n"
           << "  and kinetic energy less than 1 GeV per nucleon,\n"
           << "  for any target, if the previous cross section is\n"
           << "  not applicable;\n"
           << "- G4IonsShenCrossSection : in all remaining cases, up to\n"
           << "  projectile kinetic energy of 1 TeV per nucleon.\n";
}


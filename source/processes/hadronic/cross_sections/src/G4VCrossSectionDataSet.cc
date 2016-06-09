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
// $Id: G4VCrossSectionDataSet.cc,v 1.7.4.1 2009/03/03 11:48:00 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4VCrossSectionDataSet
//
// Author  F.W. Jones, TRIUMF, 20-JAN-97
//
// Modifications:
// 23.01.2009 V.Ivanchenko move constructor and destructor to source
//

#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

G4VCrossSectionDataSet::G4VCrossSectionDataSet() :
  verboseLevel(0)
{
  G4CrossSectionDataSetRegistry::Instance()->Register(this);
}

G4VCrossSectionDataSet::~G4VCrossSectionDataSet()
{
  G4CrossSectionDataSetRegistry::Instance()->DeRegister(this);
}

// Override this method to test for particle and isotope applicability

G4bool 
G4VCrossSectionDataSet::IsZAApplicable(const G4DynamicParticle*,
                                       G4double /*ZZ*/, G4double /*AA*/)
{
  return true;
}


G4double 
G4VCrossSectionDataSet::GetIsoCrossSection(const G4DynamicParticle* aParticle,
                                           const G4Isotope* anIsotope,
                                           G4double aTemperature)
{
  G4double ZZ = anIsotope->GetZ();
  G4double AA = anIsotope->GetN();
  return GetIsoZACrossSection(aParticle, ZZ, AA, aTemperature);
}

// Override this method to get real isotopic cross sections

G4double 
G4VCrossSectionDataSet::GetIsoZACrossSection(const G4DynamicParticle*,
                                             G4double /*ZZ*/, G4double AA,
                                             G4double /*aTemperature*/)
{
  return 62*std::pow(AA, 2./3.)*millibarn;
}

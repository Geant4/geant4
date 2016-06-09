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
// $Id: G4VCrossSectionDataSet.cc,v 1.7 2006/12/29 02:05:48 dennis Exp $
// 

#include "G4VCrossSectionDataSet.hh"

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

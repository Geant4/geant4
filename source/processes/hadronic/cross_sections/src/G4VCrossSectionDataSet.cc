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
// $Id: G4VCrossSectionDataSet.cc,v 1.12 2011-01-09 02:37:48 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4Pow.hh"

G4VCrossSectionDataSet::G4VCrossSectionDataSet(const G4String& nam) :
  verboseLevel(0),minKinEnergy(0.0),maxKinEnergy(100*TeV),name(nam) 
{
  G4CrossSectionDataSetRegistry::Instance()->Register(this);
}

G4VCrossSectionDataSet::~G4VCrossSectionDataSet()
{
  G4CrossSectionDataSetRegistry::Instance()->DeRegister(this);
}

// Override these methods to test for particle and isotope applicability
G4bool 
G4VCrossSectionDataSet::IsZAApplicable(const G4DynamicParticle*,
                                       G4double /*ZZ*/, G4double /*AA*/)
{ //obsolete method, should not be used
  static G4bool onceOnly(true);
  if ( onceOnly )
  {
      G4cerr << 
      "Warning: G4VCrossSectionDataSet::IsZAApplicable() is obsolete, invoked by " 
      <<  GetName() << G4endl;
      onceOnly=false;
  }
  return true;
}

G4bool 
G4VCrossSectionDataSet::IsIsoApplicable(const G4DynamicParticle* aPart,
					G4int Z, G4int A)
{
//  return true;
    return IsZAApplicable(aPart, Z, A);
}

G4double 
G4VCrossSectionDataSet::GetIsoCrossSection(const G4DynamicParticle* aParticle,
                                           const G4Isotope* anIsotope,
                                           G4double aTemperature)
{
  G4int ZZ = anIsotope->GetZ();
  G4int AA = anIsotope->GetN();
  return GetZandACrossSection(aParticle, ZZ, AA, aTemperature);
}

// Override this method to get real isotopic cross sections

G4double 
G4VCrossSectionDataSet::GetIsoZACrossSection(const G4DynamicParticle*,
                                             G4double /*ZZ*/, G4double AA,
                                             G4double /*aTemperature*/)
{ //obsolete method, should not be used
  static G4bool onceOnly(true);
  if ( onceOnly )
  {
      G4cerr << 
      "Warning: G4VCrossSectionDataSet::GetIsoZACrossSection() is obsolete, invoked by " 
      <<  GetName() << G4endl;
      onceOnly=false;
  }
  return 62*G4Pow::GetInstance()->A23(AA)*millibarn;
}

G4double 
G4VCrossSectionDataSet::GetZandACrossSection(const G4DynamicParticle* aPart,
                                             G4int Z, G4int N,
                                             G4double aTemperature)
{
//GF  return 62*G4Pow::GetInstance()->Z23(N)*millibarn;
  return GetIsoZACrossSection(aPart, G4double(Z), G4double(N), aTemperature);
}

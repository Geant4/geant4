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
// $Id: G4HadronCaptureDataSet.cc 66241 2012-12-13 18:34:42Z gunter $
//
//
// G4 Physics class: HadronCaptureDataSet for cross sections
// F.W. Jones, TRIUMF, 19-MAY-98
// 
// 19 Aug 2011, V.Ivanchenko move to new design and make x-section per element
//

#include "G4HadronCaptureDataSet.hh"
#include <iostream>

G4HadronCaptureDataSet::G4HadronCaptureDataSet(const G4String& nam)
 : G4VCrossSectionDataSet(nam)
{
  theHadronCrossSections = G4HadronCrossSections::Instance();
}

G4HadronCaptureDataSet::~G4HadronCaptureDataSet()
{}

G4bool
G4HadronCaptureDataSet::IsElementApplicable(const G4DynamicParticle*, 
					    G4int /*Z*/, const G4Material*)
{
  return true;
}

G4double
G4HadronCaptureDataSet::GetElementCrossSection(const G4DynamicParticle* aParticle, 
					       G4int Z, const G4Material*)
{
  return theHadronCrossSections->GetCaptureCrossSection(aParticle, Z);
}


void G4HadronCaptureDataSet::CrossSectionDescription(std::ostream& outFile) const 
{
  outFile << "G4HadronCaptureDataSet contains neutron capture cross\n"
          << "sections developed as part of the Gheisha hadronic package\n"
          << "by H. Fesefeldt.  The cross sections are valid for all\n"
          << "incident neutron energies, but they do not represent any of\n"
          << "the detailed resonances known to exist at low energies.\n"
          << "The cross sections depend only on Z and not A.\n";
}

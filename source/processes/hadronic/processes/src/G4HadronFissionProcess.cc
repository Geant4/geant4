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
// F.W. Jones, TRIUMF, 03-DEC-96
// 
// This is a prototype of a low-energy fission process.
// Currently it is based on the GHEISHA routine FISSIO,
// and conforms fairly closely to the original Fortran.
// Note: energy is in MeV and momentum is in MeV/c.
//
// 27-MAR-97 F.W.Jones: first version for Alpha release
// 20-JUN-97 F.W.Jones: added check for zero cross section
//
// 19-MAY-98 FWJ: variant G4HadronFission process for
// G4CrossSectionDataSet/DataStore class design.
// 29-JUN-98 FWJ: default data set G4HadronCrossSections
// 01-SEP-2008 V.Ivanchenko: use methods from the base class
// 14-Sep-12 M.Kelsey -- Pass subType code to base ctor
//

#include "G4HadronFissionProcess.hh"
#include "G4HadronFissionDataSet.hh"
#include "G4Neutron.hh"

G4HadronFissionProcess::G4HadronFissionProcess(const G4String& processName) : 
  G4HadronicProcess(processName,fFission)
{
  AddDataSet(new G4HadronFissionDataSet());
}


G4HadronFissionProcess::~G4HadronFissionProcess()
{}


G4bool G4HadronFissionProcess::IsApplicable(const G4ParticleDefinition& p)
{ 
  return (&p == G4Neutron::Neutron());
}


void G4HadronFissionProcess::ProcessDescription(std::ostream& outFile) const
{
  outFile << "G4HadronFissionProcess handles neutron-induced fission of nuclei\n"
          << "by invoking one or more hadronic models and one or more hadronic\n"
          << "cross sections.\n";
}

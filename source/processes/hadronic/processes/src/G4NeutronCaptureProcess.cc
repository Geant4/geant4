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
// G4 Process: Low-energy Neutron Capture
// F.W. Jones, TRIUMF, 03-DEC-96
// 
// This is a prototype of a low-energy neutron capture process.
// Currently it is based on the GHEISHA routine CAPTUR,
// and conforms fairly closely to the original Fortran.
//
// 27-MAR-97 FWJ: first version for Alpha release
// 20-JUN-97 FWJ: added check for zero cross section
//
// 19-MAY-98 FWJ: variant G4HadronCapture process for
// G4CrossSectionDataSet/DataStore class design.
// 29-JUN-98 FWJ: default data set G4HadronCrossSections
// design re-done removing most of the original parts: JPW 2003.
// 01-SEP-2008 V.Ivanchenko: use methods from the base class
// 14-Sep-12 M.Kelsey -- Pass subType code to base ctor


#include "G4NeutronCaptureProcess.hh"
#include "G4Neutron.hh"
#include "G4NeutronCaptureXS.hh"

G4NeutronCaptureProcess::G4NeutronCaptureProcess(const G4String& processName) : 
  G4HadronicProcess(processName, fCapture)
{
  AddDataSet(new G4NeutronCaptureXS());
}

G4bool
G4NeutronCaptureProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return (&aParticleType == G4Neutron::Neutron());
}

void G4NeutronCaptureProcess::ProcessDescription(std::ostream& outFile) const
{
  outFile << "G4NeutronCaptureProcess handles the capture of neutrons by nuclei\n"
	  << "following by gamma/electron de-excitation cascade. One or more\n"
          << "hadronic models and hadronic cross section sets may be invoked.\n";
}


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
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
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


#include "G4HadronCaptureProcess.hh"

G4HadronCaptureProcess::G4HadronCaptureProcess(const G4String& processName) : 
   G4HadronicProcess(processName)
{
   G4HadronicProcess::AddDataSet(new G4HadronCaptureDataSet);
}

G4HadronCaptureProcess::~G4HadronCaptureProcess()
{
}


G4double G4HadronCaptureProcess::
GetMicroscopicCrossSection(const G4DynamicParticle* aParticle,
                           const G4Element* anElement, G4double aTemp)
{
   // gives the microscopic cross section in GEANT4 internal units
   if (!G4HadronicProcess::GetCrossSectionDataStore()) {
      G4Exception("G4HadronCaptureProcess", "007", FatalException, 
                  "G4HadronCaptureProcess: no cross section data Store");
      return DBL_MIN;
   }
   return G4HadronicProcess::GetCrossSectionDataStore()
            ->GetCrossSection(aParticle, anElement, aTemp);
} 
G4bool
G4HadronCaptureProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return (aParticleType == *(G4Neutron::Neutron()));
}


void G4HadronCaptureProcess::
DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
   if (!G4HadronicProcess::GetCrossSectionDataStore()) {
      G4Exception("G4HadronCaptureProcess", "007", FatalException, 
                  "G4HadronCaptureProcess: no cross section data set");
      return;
   }
   G4HadronicProcess::GetCrossSectionDataStore()
      ->DumpPhysicsTable(aParticleType);
}

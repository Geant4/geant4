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
// G4 Process: Low-energy Elastic scattering
// F.W. Jones, TRIUMF, 04-JUN-96
// 
// This is a prototype of a low-energy elastic scattering process.
// At present it just uses some GEANT3.21/GHEISHA code, translated
// to C++ methods, to calculate cross sections and scattering events.
// Eventually these methods will be replaced by ones that are more 
// accurate in the low-energy region (< 3GeV).

// Modified by FWJ 03-DEC-96: uses G4LCrossSectionData
// Modified by FWJ    FEB-97: adapted to new tracking design
// Modified by FWJ    MAR-97: newer new tracking design
//
// Modified by FWJ 27-MAR-97: first version for Alpha release
// 20-JUN-97 FWJ: added check for zero cross section
//
// 14-APR-98 FWJ: variant G4HadronElastic process for
// G4CrossSectionDataSet/DataStore class design.
// 29-JUN-98 FWJ: default data set G4HadronCrossSections
// design re-done JPW 2003.


#include "G4HadronElasticProcess.hh"
 
G4HadronElasticProcess::
G4HadronElasticProcess(const G4String& processName) : 
      G4HadronicProcess(processName)
{
  AddDataSet(new G4HadronElasticDataSet);
}

G4HadronElasticProcess::~G4HadronElasticProcess()
{
}
 

void G4HadronElasticProcess::
BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
   if (!G4HadronicProcess::GetCrossSectionDataStore()) {
      return;
   }
   G4HadronicProcess::GetCrossSectionDataStore()->BuildPhysicsTable(aParticleType);
}


G4double G4HadronElasticProcess::
GetMicroscopicCrossSection(
      const G4DynamicParticle* aParticle, const G4Element* anElement, G4double aTemp)
{
   // gives the microscopic cross section in GEANT4 internal units
   if (!G4HadronicProcess::GetCrossSectionDataStore()) {
      G4Exception("G4HadronElasticProcess", "007", FatalException,
                  "no cross section data Store");
      return DBL_MIN;
   }
   return G4HadronicProcess::GetCrossSectionDataStore()
          ->GetCrossSection(aParticle, anElement, aTemp);
}


G4bool
G4HadronElasticProcess::
IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return (aParticleType == *(G4PionPlus::PionPlus()) ||
           aParticleType == *(G4PionZero::PionZero()) ||
           aParticleType == *(G4PionMinus::PionMinus()) ||
           aParticleType == *(G4KaonPlus::KaonPlus()) ||
           aParticleType == *(G4KaonZeroShort::KaonZeroShort()) ||
           aParticleType == *(G4KaonZeroLong::KaonZeroLong()) ||
           aParticleType == *(G4KaonMinus::KaonMinus()) ||
           aParticleType == *(G4Proton::Proton()) ||
           aParticleType == *(G4AntiProton::AntiProton()) ||
           aParticleType == *(G4Neutron::Neutron()) ||
           aParticleType == *(G4AntiNeutron::AntiNeutron()) ||
           aParticleType == *(G4Lambda::Lambda()) ||
           aParticleType == *(G4AntiLambda::AntiLambda()) ||
           aParticleType == *(G4SigmaPlus::SigmaPlus()) ||
           aParticleType == *(G4SigmaZero::SigmaZero()) ||
           aParticleType == *(G4SigmaMinus::SigmaMinus()) ||
           aParticleType == *(G4AntiSigmaPlus::AntiSigmaPlus()) ||
           aParticleType == *(G4AntiSigmaZero::AntiSigmaZero()) ||
           aParticleType == *(G4AntiSigmaMinus::AntiSigmaMinus()) ||
           aParticleType == *(G4XiZero::XiZero()) ||
           aParticleType == *(G4XiMinus::XiMinus()) ||
           aParticleType == *(G4AntiXiZero::AntiXiZero()) ||
           aParticleType == *(G4AntiXiMinus::AntiXiMinus()) ||
           aParticleType == *(G4Deuteron::Deuteron()) ||
           aParticleType == *(G4Triton::Triton()) ||
           aParticleType == *(G4Alpha::Alpha()) ||
           aParticleType == *(G4OmegaMinus::OmegaMinus()) ||
           aParticleType == *(G4AntiOmegaMinus::AntiOmegaMinus()));
}


void 
G4HadronElasticProcess::
DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
   if (!G4HadronicProcess::GetCrossSectionDataStore()) {
      G4Exception("G4HadronElasticProcess", "111", JustWarning, "G4HadronElasticProcess: no cross section data set");
      return;
   }

   G4HadronicProcess::GetCrossSectionDataStore()->DumpPhysicsTable(aParticleType);
}

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronElasticProcess.cc,v 1.2 1999-12-15 14:53:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//


#include "G4HadronElasticProcess.hh"
 
// Can this be brought into the class?
//unsigned G4HadronElasticProcessHashFun(const G4ParticleDefinition& aParticleDefinition) {
//   return aParticleDefinition.GetParticleName().hash();
//}

G4HadronElasticProcess::
G4HadronElasticProcess(const G4String& processName) : 
      G4HadronicProcess(processName),
      //      thePhysicsDictionary(G4HadronElasticProcessHashFun),
      theCrossSectionDataStore(new G4CrossSectionDataStore)
{
   theCrossSectionDataStore->AddDataSet(new G4HadronElasticDataSet);
   aParticleChange.SetNumberOfSecondaries(1);
}

G4HadronElasticProcess::~G4HadronElasticProcess()
{
   aParticleChange.Clear();
}
 

void 
G4HadronElasticProcess::
BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
   if (!theCrossSectionDataStore) {
      G4Exception("G4HadronElasticProcess: no cross section data set");
      return;
   }

   theCrossSectionDataStore->BuildPhysicsTable(aParticleType);

  //   G4double LowEdgeEnergy, Value;
  //   G4bool isOutRange;
  //
  //   static const G4ElementTable* theElementTable;
  //   theElementTable = G4Element::GetElementTable();
  //   G4int numOfElements = G4Element::GetNumberOfElements();
  //   if (verboseLevel > 1) {
  //      G4cout << "BuildPhysicsTable: numOfElements " << numOfElements << G4endl;
  //      G4cout << "BuildPhysicsTable: thePhysicsTable " << thePhysicsTable << G4endl;
  //      G4cout << "  for particle: " << aParticleType.GetParticleName() << " " <<
  //              aParticleType.GetPDGEncoding() << G4endl;
  //   }
  //
  //   RWBoolean exists = thePhysicsDictionary.contains(&aParticleType);
  //   if (verboseLevel > 1) {
  //      G4cout << "exists=" << exists << G4endl;
  //   }
  //   if (exists) {
  //      thePhysicsTable = thePhysicsDictionary.findValue(&aParticleType);
  //      if (verboseLevel > 1) {
  //         G4cout << "  Retrieve tPT=" << thePhysicsTable << G4endl;
  //      }
// Remove any previous items from table
  //      thePhysicsTable->clearAndDestroy();
  //      thePhysicsTable->resize(numOfElements);
  //   }
  //   else {
  //      thePhysicsTable = new G4PhysicsTable(numOfElements);
  //      if (verboseLevel > 1) {
  //         G4cout << "  make new tPT=" << thePhysicsTable << G4endl;
  //         G4cout << "  entries() = " << thePhysicsTable->entries() << G4endl;
  //      }
  //      thePhysicsDictionary.insertKeyAndValue(&aParticleType, thePhysicsTable);
  //   }
// Make a PhysicsVector for each element
  //   for (G4int J = 0; J < numOfElements; J++) { 
     //      (*thePhysicsTable)(J) =
     //         theCrossSectionData.MakePhysicsVector(*this,
     //                                               aParticleType,
     //                                               (*theElementTable)[J]);
  //      if (verboseLevel > 1)
  //         G4cout << "  added PhysicsVector " << (*thePhysicsTable)(J) << G4endl;
  //   }
}


G4double 
G4HadronElasticProcess::GetMicroscopicCrossSection(
      const G4DynamicParticle* aParticle, const G4Element* anElement)
{
   // gives the microscopic cross section in GEANT4 internal units

   if (!theCrossSectionDataStore) {
      G4Exception("G4HadronElasticProcess: no cross section data Store");
      return DBL_MIN;
   }
   return theCrossSectionDataStore->GetCrossSection(aParticle, anElement);

  //   G4int J ;
  //   G4bool isOutRange ;
  //   G4double MicroscopicCrossSection ;
  //
  //   RWBoolean exists = 
  //      thePhysicsDictionary.contains(aParticle->GetDefinition());
  //   if (!exists) {
  //      G4Exception("G4HadronElasticProcess::GetMicroscopicCrossSection: "
  //                  "no physics table for particle");
  //      return G4double(0);
  //   }
  //   thePhysicsTable = 
  //      thePhysicsDictionary.findValue(aParticle->GetDefinition());
  //   if (verboseLevel > 1) {
  //      G4cout << "  Retrieve tPT=" << thePhysicsTable << G4endl;
  //   }
  //
  //   J = anElement->GetIndex();
  //   MicroscopicCrossSection = (*((*thePhysicsTable)(J))).
  //                              GetValue(aParticle->GetTotalMomentum()/GeV,
  //                                       isOutRange);
  //   return MicroscopicCrossSection;
}


G4double
G4HadronElasticProcess::
GetMeanFreePathBasic(const G4DynamicParticle* aParticle,
                     const G4Material* aMaterial)
{
   // Returns the mean free path in GEANT4 internal units.
   // The cross section data is stored for individual elements, 
   // so the total macroscopic cross section is summed over the elements.
   // (Cf. G4MuPairProduction which uses the material's effective Z to
   // calculate directly the cross section.

   const G4ElementVector* theElementVector;
   const G4double* theAtomicNumDensityVector;

   G4int J = aMaterial->GetIndex();

   theElementVector = aMaterial->GetElementVector();
   theAtomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();

   G4double sigma = 0.;

   for (G4int i = 0; i < aMaterial->GetNumberOfElements(); i++) {
     sigma = sigma + theAtomicNumDensityVector[i] *
             GetMicroscopicCrossSection(aParticle, (*theElementVector)(i));
   }
   if (verboseLevel > 1)
     G4cout << "G4HadronElasticProcess::GetMeanFreePathBasic: sigma="
          << sigma << G4endl;
   if (sigma > 0.)
      return 1./sigma;
   else
      return DBL_MAX;
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
   if (!theCrossSectionDataStore) {
      G4Exception("G4HadronElasticProcess: no cross section data set");
      return;
   }

   theCrossSectionDataStore->DumpPhysicsTable(aParticleType);

  //   RWBoolean exists = thePhysicsDictionary.contains(&aParticleType);
  //   if (verboseLevel > 1) {
  //      G4cout << "DumpPhysicsTable: exists=" << exists << G4endl;
  //   }
  //   if (!exists) {
  //      G4cout << "DumpPhysicsTable: no phyics table for " << 
  //              aParticleType.GetParticleName() << G4endl;
  //      return;
  //   }
  //   G4PhysicsTable* pt = thePhysicsDictionary.findValue(&aParticleType);
  //   if (verboseLevel > 1) {
  //      G4cout << "  found PhysicsTable pt=" << pt << G4endl;
  //   }
  //
  //   const G4ElementTable* et = G4Element::GetElementTable();
  //   G4int numOfElements = G4Element::GetNumberOfElements();
  //   if (verboseLevel > 1)
  //      G4cout << "DumpPhysicsTable: numOfElements=" << numOfElements << G4endl;
  //   G4LPhysicsFreeVector* pv;
  //   for (G4int J = 0; J < numOfElements; J++) { 
  //      G4cout << G4endl << "Cross section data for " << 
  //              aParticleType.GetParticleName() << " on " <<
  //              ((*et)[J])->GetName() << G4endl << G4endl;
  //      pv = (G4LPhysicsFreeVector*)(*thePhysicsTable)(J);
  //      if (verboseLevel > 1) {
  //         G4cout << "  found PhysicsVector pv=" << pv << G4endl;
  //      }
  //      pv->DumpValues();
  //   }
}

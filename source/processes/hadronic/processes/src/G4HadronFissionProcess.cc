// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronFissionProcess.cc,v 1.1 1999-01-07 16:13:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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
//


#include "G4HadronFissionProcess.hh"

//unsigned 
//G4HadronicFissionProcessHashFun(const G4ParticleDefinition& aParticleDefinition)
//{
//   return aParticleDefinition.GetParticleName().hash();
//}

G4HadronFissionProcess::
G4HadronFissionProcess(const G4String& processName) : 
   G4HadronicProcess(processName),
   //   thePhysicsDictionary(G4HadronicFissionProcessHashFun)
   theCrossSectionDataStore(new G4CrossSectionDataStore)
{
   theCrossSectionDataStore->AddDataSet(new G4HadronFissionDataSet);
}

G4HadronFissionProcess::~G4HadronFissionProcess()
{
}
 

void 
G4HadronFissionProcess::BuildThePhysicsTable(G4ParticleDefinition& aParticleType)
{
   if (!theCrossSectionDataStore) {
      G4Exception("G4HadronFissionProcess: no cross section data set");
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
   //      G4cout << "BuildThePhysicsTable: numOfElements " << numOfElements << endl;
   //      G4cout << "BuildThePhysicsTable: thePhysicsTable " << thePhysicsTable << endl;
   //      G4cout << "  for particle: " << aParticleType.GetParticleName() << " " <<
   //              aParticleType.GetPDGEncoding() << endl;
   //   }
   //   RWBoolean exists = thePhysicsDictionary.contains(&aParticleType);
   //   if (verboseLevel > 1) {
   //      G4cout << "exists=" << exists << endl;
   //   }
   //   if (exists) {
   //      thePhysicsTable = thePhysicsDictionary.findValue(&aParticleType);
   //      if (verboseLevel > 1) {
   //         G4cout << "  Retrieve tPT=" << thePhysicsTable << endl;
   //      }
   //// Remove any previous items from table
   //      thePhysicsTable->clearAndDestroy();
   //      thePhysicsTable->resize(numOfElements);
   //   }
   //   else {
   //      thePhysicsTable = new G4PhysicsTable(numOfElements);
   //      if (verboseLevel > 1) {
   //         G4cout << "  make new tPT=" << thePhysicsTable << endl;
   //         G4cout << "  entries() = " << thePhysicsTable->entries() << endl;
   //      }
   //      thePhysicsDictionary.insertKeyAndValue(&aParticleType, thePhysicsTable);
   //   }
   //// Make a PhysicsVector for each element
   //   for (G4int J = 0; J < numOfElements; J++) { 
   //      (*thePhysicsTable)(J) =
   //         theCrossSectionData.MakePhysicsVector(*this,
   //                                               aParticleType,
   //                                               (*theElementTable)[J]);
   //   }
}


G4double 
G4HadronFissionProcess::GetMicroscopicCrossSection(const G4DynamicParticle* aParticle,
                                       const G4Element* anElement)
{
   // gives the microscopic cross section in GEANT4 internal units

   if (!theCrossSectionDataStore) {
      G4Exception("G4HadronFissionProcess: no cross section data Store");
      return DBL_MIN;
   }
   return theCrossSectionDataStore->GetCrossSection(aParticle, anElement);

   //   G4int J ;
   //   G4bool isOutRange ;
   //   G4double MicroscopicCrossSection ;
   //
   //   RWBoolean exists = 
   //     thePhysicsDictionary.contains(aParticle->GetDefinition());
   //   if (!exists) {
   //      G4Exception("G4HadronFissionProcess::GetMicroscopicCrossSection: "
   //                  "no physics table for particle");
   //      return G4double(0);
   //   }
   //   thePhysicsTable = 
   //      thePhysicsDictionary.findValue(aParticle->GetDefinition());
   //   if (verboseLevel > 1) {
   //      G4cout << "  Retrieve tPT=" << thePhysicsTable << endl;
   //   }
   //
   //   J = anElement->GetIndex();
   //   MicroscopicCrossSection = (*((*thePhysicsTable)(J))).
   //                              GetValue(aParticle->GetKineticEnergy()/GeV,
   //                                       isOutRange);
   //   return MicroscopicCrossSection;
}


G4double
G4HadronFissionProcess::GetMeanFreePathBasic(const G4DynamicParticle* aParticle,
                                 const G4Material* aMaterial)
{
   // returns the mean free path in GEANT4 internal units

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
     G4cout << "G4HadronFissionProcess::GetMeanFreePathBasic: sigma=" 
          << sigma << endl;
   if (sigma > 0.)
      return 1./sigma;
   else
      return DBL_MAX;
}
 

void
G4HadronFissionProcess::
DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
   if (!theCrossSectionDataStore) {
      G4Exception("G4HadronFissionProcess: no cross section data set");
      return;
   }

   theCrossSectionDataStore->DumpPhysicsTable(aParticleType);

   //   RWBoolean exists = thePhysicsDictionary.contains(&aParticleType);
   //   if (verboseLevel > 1) {
   //      G4cout << "DumpPhysicsTable: exists=" << exists << endl;
   //   }
   //   if (!exists) {
   //      G4Exception("G4HadronFissionProcess::DumpPhysicsTable: "
   //                  "no physics table for particle");
   //      return;
   //   }
   //   G4PhysicsTable* pt = thePhysicsDictionary.findValue(&aParticleType);
   //
   //   const G4ElementTable* et = G4Element::GetElementTable();
   //   G4int numOfElements = G4Element::GetNumberOfElements();
   //   if (verboseLevel > 1)
   //      G4cout << "DumpPhysicsTable: numOfElements=" << numOfElements << endl;
   //   G4LPhysicsFreeVector* pv;
   //   for (G4int J = 0; J < numOfElements; J++) { 
   //      G4cout << endl << "Fission cross section data for " << 
   //              aParticleType.GetParticleName() << " on " <<
   //              ((*et)[J])->GetName() << endl << endl;
   //      G4cout << endl << "Ek(GeV)    sigma(mb)" << endl << endl;
   //      pv = (G4LPhysicsFreeVector*)(*thePhysicsTable)(J);
   //      pv->DumpValues();
   //   }
}

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ElectronNucleusDataSet.hh,v 1.1 2000-09-27 07:13:06 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4ElectronNucleusDataSet -- header file
// M.V. Kossov, ITEP(Moscow), 10-SEP-00
//

#ifndef G4ElectronNucleusDataSet_h
#define G4ElectronNucleusDataSet_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4HadronCrossSections.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"


class G4ElectronNucleusDataSet : public G4VCrossSectionDataSet
{
public:

  G4ElectronNucleusDataSet()               // Constructor @@??
   {
	 //theHadronCrossSections = G4HadronCrossSections::Instance();
   }

   ~G4ElectronNucleusDataSet() {}

   G4bool IsApplicable(const G4DynamicParticle* aParticle,
                       const G4Element* anElement)
   {
	 //return theHadronCrossSections->IsApplicable(aParticle, anElement);
     // Possible prototype
     G4bool result = false;
     G4double AA=targetAtomicNumber+targetAtomicNumber;
	 G4int i=6;
	 if     (AA<ANucl(0)+ANucl(1)) i=0;
	 else if(AA<ANucl(1)+ANucl(2)) i=1;
	 else if(AA<ANucl(2)+ANucl(3)) i=2;
	 else if(AA<ANucl(3)+ANucl(4)) i=3;
	 else if(AA<ANucl(4)+ANucl(5)) i=4;
	 else if(AA<ANucl(5)+ANucl(6)) i=5;
     if( aParticle->GetDefinition()->GetPDGEncoding()==22 &&
         aParticle->GetKineticEnergy()/GeV>ThresholdEnergy(i))
       ) result = true;
   }

   G4double GetCrossSection(const G4DynamicParticle* aParticle,
                            const G4Element* anElement); //Defined in .cc
   //{
	 //return theHadronCrossSections->GetInelasticCrossSection(aParticle,
     //                                                         anElement);
   //}

   void BuildPhysicsTable(const G4ParticleDefinition&) {}

   void DumpPhysicsTable(const G4ParticleDefinition&) {}

private
  G4double GetCrossSection(const G4DynamicParticle* aPart, const G4Element* anEle);
  G4double AbsorbtionByNucleus(G4int i, G4double E);
  G4double HighEnergyOld(G4double E);
  G4double HighEnergyNew(G4double E);
  G4double LinearFit(G4double X, G4int N, const G4double& XN, const G4double& YN);
  G4double SitSint(G4double E, G4double p1, G4double p2, G4double p3, G4double p4);
  G4double CorA(G4double E, G4int i);
  G4double ANucl(G4int i);
  G4int    ZNucl(G4int i);
  G4int    NNucl(G4int i);
  G4double SumQQ(G4int i);
  G4double ThresholdEnergy(G4int i);
  G4double PbDataEnergy(G4int i);

// Body
//private:

  //G4HadronCrossSections* theHadronCrossSections;
};


// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ElectronNucleusDataSet.hh,v 1.2 2000-09-27 12:23:27 mkossov Exp $
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
     G4double targetAtomicNumber = anElement->GetN();
     //G4double nTargetProtons = anElement->GetZ();
     //G4double nTargetNeutrons = targetAtomicNumber-nTargetProtons;
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
         aParticle->GetKineticEnergy()/GeV>ThresholdEnergy(i)
       ) result = true;
   }

   G4double GetCrossSection(const G4DynamicParticle* aParticle,
                            const G4Element* anElement); //Defined in .cc
   //{
	 //return theHadronCrossSections->GetInelasticCrossSection(aParticle,
     //                                                         anElement);
   //}

  void BuildPhysicsTable(const G4ParticleDefinition&) {};

  void DumpPhysicsTable(const G4ParticleDefinition&) {};

private:
  G4double ANucl(G4int i);
  G4int    ZNucl(G4int i);
  G4int    NNucl(G4int i);
  G4double SumQQ(G4int i);
  G4double ThresholdEnergy(G4int i);

// Body
//private:
//G4HadronCrossSections* theHadronCrossSections;
};

// "A" of measured nuclei @@ Move to header
inline G4double G4ElectronNucleusDataSet::ANucl(G4int i)
{
  static const G4int n  = 7;
  static G4double AN[n] = {1.008, 2.01, 9.0122, 12.011, 26.982, 63.546, 207.2};
  return AN[i];
}

// "Z" of measured nuclei @@ Move to header
inline G4int G4ElectronNucleusDataSet::ZNucl(G4int i)
{
  static const G4int n  = 7;
  static G4int ZN[n] = {1, 1, 4, 6, 13, 29,  82};
  return ZN[i];
}

// "N" of measured nuclei @@ Move to header
inline G4int G4ElectronNucleusDataSet::NNucl(G4int i)
{
  static const G4int n  = 7;
  static G4int NN[n] = {0, 1, 5, 6, 14, 35, 125};
  return NN[i];
}

// "3*sum(q^2)" for nuclei
inline G4double G4ElectronNucleusDataSet::SumQQ(G4int i)
{
  static const G4int n  = 7;
  static G4double QQ[n] = {3., 5., 22., 30., 67., 157., 496.};
  return QQ[i];
}

// Gives the threshold energy for different nuclei @@ Move to header
inline G4double G4ElectronNucleusDataSet::ThresholdEnergy(G4int i)
{
  static const G4int n = 7;
  static G4double T[n] = {0.200, 0.0022, 0.012, 0.018, 0.013, 0.010, 0.007};
  return T[i];
}

#endif

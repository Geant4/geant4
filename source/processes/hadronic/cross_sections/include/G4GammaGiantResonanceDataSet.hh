// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GammaGiantResonanceDataSet.hh,v 1.2 2000-09-27 12:23:27 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4GammaGiantResonanceDataSet -- header file
// M.V. Kossov, ITEP(Moscow), 10-SEP-00
//

#ifndef G4GammaGiantResonanceDataSet_h
#define G4GammaGiantResonanceDataSet_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4HadronCrossSections.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"


class G4GammaGiantResonanceDataSet : public G4VCrossSectionDataSet
{
public:

  G4GammaGiantResonanceDataSet()               // Constructor @@??
   {
	 //theHadronCrossSections = G4HadronCrossSections::Instance();
   }

   ~G4GammaGiantResonanceDataSet() {}

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
         aParticle->GetKineticEnergy()/MeV<ThresholdEnergy(i)
       ) result = true;
   }

   G4double GetCrossSection(const G4DynamicParticle* aParticle,
                            const G4Element* anElement);
   //{
	 //return theHadronCrossSections->GetInelasticCrossSection(aParticle,
     //                                                         anElement);
   //}

   void BuildPhysicsTable(const G4ParticleDefinition&) {}

   void DumpPhysicsTable(const G4ParticleDefinition&) {}

private:
  G4double AbsorbtionByNucleus(G4int i, G4double E);
  G4double HighEnergyOld(G4double E);
  G4double HighEnergyNew(G4double E);
  G4double LinearFit(G4double X, G4int N, const G4double* XN, const G4double* YN);
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

// Working function @@ Move to header
inline G4double G4GammaGiantResonanceDataSet::SitSint(G4double E, G4double p1,
                                        G4double p2, G4double p3, G4double p4)
{
  G4double at=atan(pow(E-p2,p4)*p3);
  return 1.+p1*sin(at+at);
}

// Correction function for Be,C @@ Move to header
inline G4double G4GammaGiantResonanceDataSet::CorA(G4double E, G4int i)
{
  static const G4int n = 7;
  static G4double p1[n] = {0., 0., 0.107630, 0.3149300, 0.3717400, 0.21533, 0.16519};
  static G4double p2[n] = {0., 0., 0.014351,-0.0028181, 0.0084629, 0.01033, 0.0066295};
  static G4double p3[n] = {0., 0.,   1425.6,    2250.2,    229.32,  46.911, 389.67};
  static G4double p4[n] = {0., 0.,   1.7264,    2.6033,    1.6613,  1.1760, 1.7950};
  if(E<0.16&&E>ThresholdEnergy(i)) return SitSint(E, p1[i], p2[i], p3[i], p4[i]);
  else                             return 1.;
}

// "A" of measured nuclei @@ Move to header
inline G4double G4GammaGiantResonanceDataSet::ANucl(G4int i)
{
  static const G4int n  = 7;
  static G4double AN[n] = {1.008, 2.01, 9.0122, 12.011, 26.982, 63.546, 207.2};
  return AN[i];
}

// "Z" of measured nuclei @@ Move to header
inline G4int G4GammaGiantResonanceDataSet::ZNucl(G4int i)
{
  static const G4int n  = 7;
  static G4int ZN[n] = {1, 1, 4, 6, 13, 29,  82};
  return ZN[i];
}

// "N" of measured nuclei @@ Move to header
inline G4int G4GammaGiantResonanceDataSet::NNucl(G4int i)
{
  static const G4int n  = 7;
  static G4int NN[n] = {0, 1, 5, 6, 14, 35, 125};
  return NN[i];
}

// "3*sum(q^2)" for nuclei
inline G4double G4GammaGiantResonanceDataSet::SumQQ(G4int i)
{
  static const G4int n  = 7;
  static G4double QQ[n] = {3., 5., 22., 30., 67., 157., 496.};
  return QQ[i];
}

// Gives the threshold energy for different nuclei @@ Move to header
inline G4double G4GammaGiantResonanceDataSet::ThresholdEnergy(G4int i)
{
  static const G4int n = 7;
  static G4double T[n] = {0.200, 0.0022, 0.012, 0.018, 0.013, 0.010, 0.007};
  return T[i];
}

// Gives the threshold energy of Pb data use for different nuclei @@ Move to header
inline G4double G4GammaGiantResonanceDataSet::PbDataEnergy(G4int i)
{
  static const G4int n = 7;
  static G4double T[n] = {0.227, 0.0022, 0.032, 0.038, 0.035, 0.033, 0.027};
  return T[i];
}

#endif

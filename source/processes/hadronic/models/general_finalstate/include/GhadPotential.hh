#ifndef GhadPotential_h
#define GhadPotential_h

#include "G4FermiMomentum.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4AntiKaonZero.hh"
#include "G4KaonZero.hh"
#include "G4AntiProton.hh"

class GhadPotential
{
  public:
  
  GhadPotential(G4double anA, G4double aZ) : theA(anA), theZ(aZ)
  {
   if(anA < 17) theDencity = new G4NuclearShellModelDensity(anA, aZ);
   else         theDencity = new G4NuclearFermiDensity(anA, aZ);
   theFermi.Init(anA, aZ);
  }
 
  G4double GetPotential(G4ParticleDefinition * aN, G4double radius)
  {
    if(dynamic_cast<G4Neutron *>(aN)) return GetPotential(G4Neutron::Neutron(), radius);
    if(dynamic_cast<G4Proton *>(aN)) return GetPotential(G4Proton::Proton(), radius);
    if(dynamic_cast<G4PionPlus *>(aN)) return GetPotential(G4PionPlus::PionPlus(), radius);
    if(dynamic_cast<G4PionMinus *>(aN)) return GetPotential(G4PionMinus::PionMinus(), radius);
    if(dynamic_cast<G4PionZero *>(aN)) return GetPotential(G4PionZero::PionZero(), radius);
    if(dynamic_cast<G4KaonMinus *>(aN)) return GetPotential(G4KaonMinus::KaonMinus(), radius);
    if(dynamic_cast<G4KaonPlus *>(aN)) return GetPotential(G4KaonPlus::KaonPlus(), radius);
    if(dynamic_cast<G4KaonZero *>(aN)) return GetPotential(G4KaonZero::KaonZero(), radius);
    if(dynamic_cast<G4AntiKaonZero *>(aN)) return GetPotential(G4AntiKaonZero::AntiKaonZero(), radius);
    if(dynamic_cast<G4AntiProton *>(aN)) return GetPotential(G4AntiProton::AntiProton(), radius);
    return 0;
  }
  
  private:
  G4double GetPotential(G4Neutron * aN, G4double radius);
  G4double GetPotential(G4Proton * aP, G4double radius);
  G4double GetPotential(G4AntiProton * aP, G4double radius);
  G4double GetPotential(G4KaonPlus * aK, G4double radius);
  G4double GetPotential(G4KaonMinus * aK, G4double radius);
  G4double GetPotential(G4KaonZero * aK, G4double radius);
  G4double GetPotential(G4AntiKaonZero * aK, G4double radius);
  G4double GetPotential(G4PionMinus * aPi, G4double radius);
  G4double GetPotential(G4PionPlus * aPi, G4double radius);
  G4double GetPotential(G4PionZero * aPi, G4double radius);
  G4double GetKaonPotential(G4double radius);
  G4double GetPionPotential(G4double radius);
  
  G4VNuclearDensity *theDencity;
  G4FermiMomentum theFermi;
  G4double theA;
  G4double theZ;

  // Coulomb barrier
  static const G4double coulomb = 1.44 / 1.14 * MeV;
  
  // pion's potential constant (real part only)
  // 0.35 + i0.82 or 0.63 + i0.89 fermi
  static const G4double a_pion = 0.35;

  // kaon's potential constant (real part only)
  // @@@@@@@ correct @@@@@ for kaon it should differ from pions
  // 0.35 + i0.82 or 0.63 + i0.89 fermi
  static const G4double a_kaon = 0.35;

  // antiproton's potential constant (real part only)
  // 1.53 + i2.50 fermi
  static const G4double a_antiproton = 1.53;

};

#endif

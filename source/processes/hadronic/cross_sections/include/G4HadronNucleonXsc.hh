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
// Calculation of the total, elastic and inelastic cross-sections
// based on parametrisations of (proton, pion, kaon, photon) nucleon
// cross-sections and the hadron-nucleous cross-section model in 
// the framework of Glauber-Gribov approach
//
//
//
// 14.03.07 V. Grichine - first implementation
// 04.11.11 V. Grichine - update for kaon-(p,n) xsc, vector spline
// 21.02.12 V. Grichine - update for pion-(p,n) xsc, NS fit++, vector spline
//
//

#ifndef G4HadronNucleonXsc_h
#define G4HadronNucleonXsc_h

#include "globals.hh"
#include "G4Proton.hh"
#include "G4Nucleus.hh"
#include "G4LPhysicsFreeVector.hh"


class G4ParticleDefinition;
class G4DynamicParticle;

class G4HadronNucleonXsc 
{
public:

  G4HadronNucleonXsc ();
  virtual ~G4HadronNucleonXsc ();
   
  virtual
  G4bool IsApplicable(const G4DynamicParticle* aDP, const G4Element*);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle* aDP, G4int Z, G4int A);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&) 
  {G4cout << "G4HadronNucleonXsc: uses parametrisation"<<G4endl;}

  void CrossSectionDescription(std::ostream&) const;

  // Xsc parametrisations

  G4double GetHadronNucleonXscEL(const G4DynamicParticle*, const G4ParticleDefinition*);

  G4double GetHadronNucleonXscPDG(const G4DynamicParticle*, const G4ParticleDefinition*);

  G4double GetHadronNucleonXscNS(const G4DynamicParticle*, const G4ParticleDefinition*);

  G4double GetKaonNucleonXscGG(const G4DynamicParticle*, const G4ParticleDefinition*);

  G4double GetHadronNucleonXscVU(const G4DynamicParticle*, const G4ParticleDefinition*);

  // kinematics and set/get

  G4double CalculateEcmValue ( const G4double , const G4double , const G4double ); 

  G4double CalcMandelstamS( const G4double , const G4double , const G4double );

  G4double GetCoulombBarrier(const G4DynamicParticle* aParticle, const G4ParticleDefinition* nucleon );

  G4double GetTotalHadronNucleonXsc()    { return fTotalXsc;     }; 
  G4double GetElasticHadronNucleonXsc()  { return fElasticXsc;   }; 
  G4double GetInelasticHadronNucleonXsc(){ return fInelasticXsc; }; 

  /*
  void InitialiseKaonNucleonTotXsc();
  G4double GetKpProtonTotXscVector(G4double logEnergy){ return fKpProtonTotXscVector.Value(logEnergy); };
  G4double GetKpNeutronTotXscVector(G4double logEnergy){ return fKpNeutronTotXscVector.Value(logEnergy); };
  G4double GetKmProtonTotXscVector(G4double logEnergy){ return fKmProtonTotXscVector.Value(logEnergy); };
  G4double GetKmNeutronTotXscVector(G4double logEnergy){ return fKmNeutronTotXscVector.Value(logEnergy); };
  */

private:

//  const G4double fUpperLimit;
  const G4double fLowerLimit; 

  G4double fTotalXsc, fElasticXsc, fInelasticXsc;
//  G4double fHadronNucleonXsc;
  /*
  // K-nucleon tot xsc (mb) fit data, std::log(Tkin(MeV))

  static const G4double fKpProtonTotXsc[66];
  static const G4double fKpProtonTotTkin[66];

  static const G4double fKpNeutronTotXsc[75];
  static const G4double fKpNeutronTotTkin[75];

  static const G4double fKmProtonTotXsc[106];
  static const G4double fKmProtonTotTkin[106];

  static const G4double fKmNeutronTotXsc[68];
  static const G4double fKmNeutronTotTkin[68];

  G4LPhysicsFreeVector fKpProtonTotXscVector;
  G4LPhysicsFreeVector fKpNeutronTotXscVector;
  G4LPhysicsFreeVector fKmProtonTotXscVector;
  G4LPhysicsFreeVector fKmNeutronTotXscVector;
  */ 
  G4ParticleDefinition* theGamma;
  G4ParticleDefinition* theProton;
  G4ParticleDefinition* theNeutron;
  G4ParticleDefinition* theAProton;
  G4ParticleDefinition* theANeutron;
  G4ParticleDefinition* thePiPlus;
  G4ParticleDefinition* thePiMinus;
  G4ParticleDefinition* thePiZero;
  G4ParticleDefinition* theKPlus;
  G4ParticleDefinition* theKMinus;
  G4ParticleDefinition* theK0S;
  G4ParticleDefinition* theK0L;
  G4ParticleDefinition* theL;
  G4ParticleDefinition* theAntiL;
  G4ParticleDefinition* theSPlus;
  G4ParticleDefinition* theASPlus;
  G4ParticleDefinition* theSMinus;
  G4ParticleDefinition* theASMinus;
  G4ParticleDefinition* theS0;
  G4ParticleDefinition* theAS0;
  G4ParticleDefinition* theXiMinus;
  G4ParticleDefinition* theXi0;
  G4ParticleDefinition* theAXiMinus;
  G4ParticleDefinition* theAXi0;
  G4ParticleDefinition* theOmega;
  G4ParticleDefinition* theAOmega;
  G4ParticleDefinition* theD;
  G4ParticleDefinition* theT;
  G4ParticleDefinition* theA;
  G4ParticleDefinition* theHe3;

};


#endif

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
//
//
// 17.07.06 V. Grichine - first implementation
// 22.01.07 V.Ivanchenko - add interface with Z and A
// 05.03.07 V.Ivanchenko - add IfZAApplicable
// 06.03.07 V.Ivanchenko - add GetElasticGlauberGribov and GetElasticGlauberGribov
//                         for combined dataset
//
//

#ifndef G4GlauberGribovCrossSection_h
#define G4GlauberGribovCrossSection_h

#include "globals.hh"
#include "G4Proton.hh"
#include "G4Nucleus.hh"
#include "G4HadTmpUtil.hh"

#include "G4VCrossSectionDataSet.hh"

class G4ParticleDefinition;

class G4GlauberGribovCrossSection : public G4VCrossSectionDataSet
{
public:

  G4GlauberGribovCrossSection ();
  virtual ~G4GlauberGribovCrossSection ();
   
  virtual
  G4bool IsApplicable(const G4DynamicParticle* aDP, const G4Element*);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle* aDP, G4int Z, G4int A);

  virtual
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, 
			   G4double aTemperature = 0.0);

  virtual
  G4double GetZandACrossSection(const G4DynamicParticle*,
                                G4int Z, G4int A,
                                G4double aTemperature = 0.0);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&)
  {}

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&) 
  {G4cout << "G4GlauberGribovCrossSection: uses Glauber-Gribov formula"<<G4endl;}

  G4double GetRatioSD(const G4DynamicParticle*, G4int At, G4int Zt);
  G4double GetRatioQE(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double GetHadronNucleonXsc(const G4DynamicParticle*, const G4Element*);
  G4double GetHadronNucleonXsc(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double GetHadronNucleonXscPDG(const G4DynamicParticle*, const G4Element*);
  G4double GetHadronNucleonXscPDG(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double GetHadronNucleonXscNS(const G4DynamicParticle*, const G4Element*);
  G4double GetHadronNucleonXscNS(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double GetHNinelasticXsc(const G4DynamicParticle*, const G4Element*);
  G4double GetHNinelasticXsc(const G4DynamicParticle*, G4int At, G4int Zt);
  G4double GetHNinelasticXscVU(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double CalculateEcmValue ( const G4double , const G4double , const G4double ); 

  G4double CalcMandelstamS( const G4double , const G4double , const G4double );

  G4double GetElasticGlauberGribov(const G4DynamicParticle*, G4int Z, G4int A);
  G4double GetInelasticGlauberGribov(const G4DynamicParticle*, G4int Z, G4int A);

  G4double GetTotalGlauberGribovXsc()    { return fTotalXsc;     }; 
  G4double GetElasticGlauberGribovXsc()  { return fElasticXsc;   }; 
  G4double GetInelasticGlauberGribovXsc(){ return fInelasticXsc; }; 
  G4double GetProductionGlauberGribovXsc(){ return fProductionXsc; }; 
  G4double GetDiffractionGlauberGribovXsc(){ return fDiffractionXsc; }; 
  G4double GetRadiusConst()              { return fRadiusConst;  }; 

  G4double GetNucleusRadius(const G4DynamicParticle*, const G4Element*);
  G4double GetNucleusRadius(G4int At);

  inline G4double GetParticleBarCorTot(const G4ParticleDefinition* theParticle, G4int Z);
  inline G4double GetParticleBarCorIn(const G4ParticleDefinition* theParticle, G4int Z);

  inline void SetEnergyLowerLimit(G4double E ){fLowerLimit=E;};

private:

  const G4double fUpperLimit;
  G4double fLowerLimit; 
  const G4double fRadiusConst;

  static const G4double fNeutronBarCorrectionTot[93];
  static const G4double fNeutronBarCorrectionIn[93];

  static const G4double fProtonBarCorrectionTot[93];
  static const G4double fProtonBarCorrectionIn[93];

  static const G4double fPionPlusBarCorrectionTot[93];
  static const G4double fPionPlusBarCorrectionIn[93];

  static const G4double fPionMinusBarCorrectionTot[93];
  static const G4double fPionMinusBarCorrectionIn[93];

  G4double fTotalXsc, fElasticXsc, fInelasticXsc, fProductionXsc, fDiffractionXsc;
  G4double fHadronNucleonXsc;
 
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

////////////////////////////////////////////////////////////////
//
// Inlines

inline
G4double
G4GlauberGribovCrossSection::GetElasticGlauberGribov(const G4DynamicParticle* dp,
                                                     G4int Z, G4int A)
{
  GetZandACrossSection(dp, Z, A);
  return fElasticXsc;
}

/////////////////////////////////////////////////////////////////

inline
G4double
G4GlauberGribovCrossSection::GetInelasticGlauberGribov(const G4DynamicParticle* dp,
                                                       G4int Z, G4int A)
{
  GetZandACrossSection(dp, Z, A);
  return fInelasticXsc;
}

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov tot xsc, when it 
// is available, else return 1.0


inline G4double G4GlauberGribovCrossSection::GetParticleBarCorTot( 
                          const G4ParticleDefinition* theParticle, G4int Z)
{
  if(Z >= 2 && Z <= 92)
  {
    if(      theParticle == theProton ) return fProtonBarCorrectionTot[Z]; 
    else if( theParticle == theNeutron) return fNeutronBarCorrectionTot[Z]; 
    else if( theParticle == thePiPlus ) return fPionPlusBarCorrectionTot[Z];
    else if( theParticle == thePiMinus) return fPionMinusBarCorrectionTot[Z];
    else return 1.0;
  }
  else return 1.0;
}

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov in xsc, when it 
// is available, else return 1.0


inline G4double G4GlauberGribovCrossSection::GetParticleBarCorIn( 
                          const G4ParticleDefinition* theParticle, G4int Z)
{
  if(Z >= 2 && Z <= 92)
  {
    if(      theParticle == theProton ) return fProtonBarCorrectionIn[Z]; 
    else if( theParticle == theNeutron) return fNeutronBarCorrectionIn[Z]; 
    else if( theParticle == thePiPlus ) return fPionPlusBarCorrectionIn[Z];
    else if( theParticle == thePiMinus) return fPionMinusBarCorrectionIn[Z];
    else return 1.0;
  }
  else return 1.0;
}

#endif

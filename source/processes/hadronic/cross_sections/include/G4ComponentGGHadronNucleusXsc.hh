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
// 25.04.12 V. Grichine - first implementation based on G4GlauberGribovCrossSection old interface
//
//

#ifndef G4ComponentGGHadronNucleusXsc_h
#define G4ComponentGGHadronNucleusXsc_h 1

#include "globals.hh"
#include "G4Proton.hh"
#include "G4Nucleus.hh"

#include "G4VComponentCrossSection.hh"

class G4ParticleDefinition;
class G4HadronNucleonXsc;

class G4ComponentGGHadronNucleusXsc : public G4VComponentCrossSection
{
public:

  G4ComponentGGHadronNucleusXsc ();
  virtual ~G4ComponentGGHadronNucleusXsc ();

  static const char* Default_Name() { return "Glauber-Gribov"; }

  // virtual interface methods

virtual
  G4double GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy,
				       G4int Z, G4int A);

virtual
  G4double GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy, 
				       G4int Z, G4double A);

virtual
  G4double GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int A);

virtual
  G4double GetProductionIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int A);

virtual
  G4double GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double A);
virtual
  G4double GetProductionElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double A);

virtual
  G4double GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4double A);

virtual
  G4double GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int A);
 
virtual
  G4double ComputeQuasiElasticRatio(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int A);
 

   
  //  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle* aDP, G4int Z, G4int A, 
			 const G4Element* elm = 0,
			 const G4Material* mat = 0);

  //  virtual
  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,  
			      const G4Isotope* iso = 0,
			      const G4Element* elm = 0,
			      const G4Material* mat = 0);

  G4double GetRatioSD(const G4DynamicParticle*, G4int At, G4int Zt);
  G4double GetRatioQE(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double GetHadronNucleonXsc(const G4DynamicParticle*, const G4Element*);
  G4double GetHadronNucleonXsc(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double GetHadronNucleonXscPDG(const G4DynamicParticle*, const G4Element*);
  G4double GetHadronNucleonXscPDG(const G4DynamicParticle*, G4int At, G4int Zt);
  G4double GetHadronNucleonXscNS(const G4DynamicParticle*, const G4Element*);
  G4double GetHadronNucleonXscNS(const G4DynamicParticle*, G4int At, G4int Zt);

  // G4double GetKaonNucleonXscVector(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double GetHNinelasticXsc(const G4DynamicParticle*, const G4Element*);
  G4double GetHNinelasticXsc(const G4DynamicParticle*, G4int At, G4int Zt);
  G4double GetHNinelasticXscVU(const G4DynamicParticle*, G4int At, G4int Zt);

  G4double CalculateEcmValue ( const G4double , const G4double , const G4double ); 

  G4double CalcMandelstamS( const G4double , const G4double , const G4double );

  G4double GetNucleusRadius(const G4DynamicParticle*, const G4Element*);
  G4double GetNucleusRadius(G4int At);

  G4double GetAxsc2piR2(){return fAxsc2piR2;};
  G4double GetModelInLog(){return fModelInLog;};

  virtual void CrossSectionDescription(std::ostream&) const;

  inline G4double GetElasticGlauberGribov(const G4DynamicParticle*, G4int Z, G4int A);
  inline G4double GetInelasticGlauberGribov(const G4DynamicParticle*, G4int Z, G4int A);

  inline G4double GetTotalGlauberGribovXsc()    { return fTotalXsc;     }; 
  inline G4double GetElasticGlauberGribovXsc()  { return fElasticXsc;   }; 
  inline G4double GetInelasticGlauberGribovXsc(){ return fInelasticXsc; }; 
  inline G4double GetProductionGlauberGribovXsc(){ return fProductionXsc; }; 
  inline G4double GetDiffractionGlauberGribovXsc(){ return fDiffractionXsc; }; 
  inline G4double GetRadiusConst()              { return fRadiusConst;  }; 

  inline G4double GetParticleBarCorTot(const G4ParticleDefinition* theParticle, G4int Z);
  inline G4double GetParticleBarCorIn(const G4ParticleDefinition* theParticle, G4int Z);

  inline void SetEnergyLowerLimit(G4double E ){fLowerLimit=E;};

private:

//  const G4double fUpperLimit;
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
  G4double fAxsc2piR2, fModelInLog;
//  G4double fHadronNucleonXsc;
 
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

  G4HadronNucleonXsc* hnXsc;

};

////////////////////////////////////////////////////////////////
//
// Inlines

inline
G4double
G4ComponentGGHadronNucleusXsc::GetElasticGlauberGribov(const G4DynamicParticle* dp,
                                                     G4int Z, G4int A)
{
  GetIsoCrossSection(dp, Z, A);
  return fElasticXsc;
}

/////////////////////////////////////////////////////////////////

inline
G4double
G4ComponentGGHadronNucleusXsc::GetInelasticGlauberGribov(const G4DynamicParticle* dp,
                                                       G4int Z, G4int A)
{
  GetIsoCrossSection(dp, Z, A);
  return fInelasticXsc;
}

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov tot xsc, when it 
// is available, else return 1.0


inline G4double G4ComponentGGHadronNucleusXsc::GetParticleBarCorTot( 
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


inline G4double G4ComponentGGHadronNucleusXsc::GetParticleBarCorIn( 
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

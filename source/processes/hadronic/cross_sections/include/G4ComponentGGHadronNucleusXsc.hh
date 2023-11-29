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
// 25.04.12 V. Grichine - first implementation based on 
//                        G4GlauberGribovCrossSection old interface
//
// 04.09.18 V. Ivantchenko Major revision of interfaces and implementation
// 01.10.18 V. Grichine strange hyperon xsc
// 27.05.19 V. Ivantchenko Removed obsolete methods and members 

#ifndef G4ComponentGGHadronNucleusXsc_h
#define G4ComponentGGHadronNucleusXsc_h 1

#include "globals.hh"
#include "G4Proton.hh"
#include "G4Nucleus.hh"

#include "G4VComponentCrossSection.hh"

class G4ParticleDefinition;
class G4HadronNucleonXsc;
class G4Pow;

class G4ComponentGGHadronNucleusXsc final : public G4VComponentCrossSection
{
public:

  explicit G4ComponentGGHadronNucleusXsc();
  ~G4ComponentGGHadronNucleusXsc() final;

  static const char* Default_Name() { return "Glauber-Gribov"; }

  // virtual interface methods
  G4double GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy, 
				       G4int Z, G4double A) final;

  G4double GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy,
				       G4int Z, G4int A) final;

  G4double GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double A) final;

  G4double GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int A) final;

  G4double GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4double A) final;

  G4double GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int A) final;
 
  G4double ComputeQuasiElasticRatio(const G4ParticleDefinition* aParticle,
				    G4double kinEnergy, 
				    G4int Z, G4int A) final;

  // Glauber-Gribov cross section
  void ComputeCrossSections(const G4ParticleDefinition* aParticle,
			    G4double kinEnergy, G4int Z, G4int A, G4int nL = 0);

  // additional public methods
  G4double GetProductionElementCrossSection(const G4ParticleDefinition* aParticle,
					    G4double kinEnergy, 
					    G4int Z, G4double A);

  G4double GetProductionIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					    G4double kinEnergy, 
					    G4int Z, G4int A);

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

  void Description(std::ostream&) const final;

  inline G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,  
	 		             const G4Isotope* iso = nullptr,
			             const G4Element* elm = nullptr,
			             const G4Material* mat = nullptr);

  inline G4double GetElasticGlauberGribov(const G4DynamicParticle*, G4int Z, G4int A);
  inline G4double GetInelasticGlauberGribov(const G4DynamicParticle*, G4int Z, G4int A);

  inline G4double GetAxsc2piR2() const                   { return fAxsc2piR2; };
  inline G4double GetModelInLog() const                  { return fModelInLog; };
  inline G4double GetTotalGlauberGribovXsc() const       { return fTotalXsc; }; 
  inline G4double GetElasticGlauberGribovXsc() const     { return fElasticXsc; }; 
  inline G4double GetInelasticGlauberGribovXsc() const   { return fInelasticXsc; }; 
  inline G4double GetProductionGlauberGribovXsc() const  { return fProductionXsc; }; 
  inline G4double GetDiffractionGlauberGribovXsc() const { return fDiffractionXsc; }; 

  inline G4double GetParticleBarCorTot(const G4ParticleDefinition* theParticle, G4int Z);
  inline G4double GetParticleBarCorIn(const G4ParticleDefinition* theParticle, G4int Z);

private:

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
  G4double fEnergy; //Cache
 
  const G4ParticleDefinition* theGamma;
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;
  const G4ParticleDefinition* theAProton;
  const G4ParticleDefinition* theANeutron;
  const G4ParticleDefinition* thePiPlus;
  const G4ParticleDefinition* thePiMinus;
  const G4ParticleDefinition* theKPlus;
  const G4ParticleDefinition* theKMinus;
  const G4ParticleDefinition* theK0S;
  const G4ParticleDefinition* theK0L;
  const G4ParticleDefinition* theLambda;

  G4HadronNucleonXsc* hnXsc;

  // Cache
  const G4ParticleDefinition* fParticle;
  G4int fZ, fA, fL;

};

////////////////////////////////////////////////////////////////
//
// Inlines

inline G4double 
G4ComponentGGHadronNucleusXsc::GetIsoCrossSection(const G4DynamicParticle* dp, 
						  G4int Z, G4int A,  
						  const G4Isotope*,
						  const G4Element*,
						  const G4Material*)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z, A);
  return fTotalXsc;
}

inline G4double
G4ComponentGGHadronNucleusXsc::GetElasticGlauberGribov(
           const G4DynamicParticle* dp, G4int Z, G4int A)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z, A);
  return fElasticXsc;
}

/////////////////////////////////////////////////////////////////

inline G4double
G4ComponentGGHadronNucleusXsc::GetInelasticGlauberGribov(
           const G4DynamicParticle* dp, G4int Z, G4int A)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z, A);
  return fInelasticXsc;
}

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov tot xsc, when it 
// is available, else return 1.0

inline G4double G4ComponentGGHadronNucleusXsc::GetParticleBarCorTot( 
                           const G4ParticleDefinition* theParticle, G4int ZZ)
{
  G4double cor = 1.0;
  G4int z = std::min(92, std::max(ZZ, 1));
  if(      theParticle == theProton ) cor = fProtonBarCorrectionTot[z]; 
  else if( theParticle == theNeutron) cor = fNeutronBarCorrectionTot[z]; 
  else if( theParticle == thePiPlus ) cor = fPionPlusBarCorrectionTot[z];
  else if( theParticle == thePiMinus) cor = fPionMinusBarCorrectionTot[z];
  return cor;
}

/////////////////////////////////////////////////////////////////////
//
// return correction at Tkin = 90*GeV GG -> Barashenkov in xsc, when it 
// is available, else return 1.0


inline G4double G4ComponentGGHadronNucleusXsc::GetParticleBarCorIn( 
                           const G4ParticleDefinition* theParticle, G4int ZZ)
{
  G4double cor = 1.0;
  G4int z = std::min(92, std::max(ZZ, 1));
  if(      theParticle == theProton ) cor = fProtonBarCorrectionIn[z]; 
  else if( theParticle == theNeutron) cor = fNeutronBarCorrectionIn[z]; 
  else if( theParticle == thePiPlus ) cor = fPionPlusBarCorrectionIn[z];
  else if( theParticle == thePiMinus) cor = fPionMinusBarCorrectionIn[z];
  return cor;
}

#endif

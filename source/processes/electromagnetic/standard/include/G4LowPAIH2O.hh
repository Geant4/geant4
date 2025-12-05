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
// G4LowPAIxs.hh -- header file
//
// GEANT 4 class header file --- Copyright CERN 1995
// CERB Geneva Switzerland
//
// for information related to this code, please, contact
// CERN, CN Division, ASD Group
//
// Preparation of ionizing collision cross section according to Photo Absorption 
// Ionization (PAI) model for simulation of ionization energy losses in very thin
// layers of water for low energy protons and electrons.
//
// Author: Vladimir.Grichine@cern.ch
//
// History:
//
 
// 2.12.24, V. Grichine: 1st version

#ifndef G4LOWPAIH2O_HH
#define G4LOWPAIH2O_HH

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4PhysicsLogVector.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"

class G4MaterialCutsCouple;
class G4Material;
class G4SandiaTable;
class G4PhysicsLogVector;
class G4PhysicsTable;
class G4ParticleDefinition;
class G4ParticleChangeForLoss;

class G4LowPAIH2O :
   public G4VEmModel  , public G4VEmFluctuationModel
{
public:
  // Constructors

  explicit G4LowPAIH2O( const G4ParticleDefinition* p = nullptr,
	       const G4String& nam = "lowpaih2o"); 
	  
 ~G4LowPAIH2O() override;

  G4LowPAIH2O & operator=(const G4LowPAIH2O &right) = delete;
  G4LowPAIH2O(const G4LowPAIH2O&) = delete;

  // methods

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override; 

  void InitialiseLocal(const G4ParticleDefinition*,
                       G4VEmModel* masterModel) override; 

  G4double CrossSectionPerVolume(const G4Material*,
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy,
                                 G4double cutEnergy,
                                 G4double maxEnergy) override; 
    
  G4double CrossSectionPerAtom(
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy, G4double Z,
                                      G4double A,
                                 G4double cutEnergy,
				 G4double maxEnergy);

  virtual G4double ComputeCrossSectionPerElectron(
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy,
                                 G4double cutEnergy,
                                 G4double maxEnergy);
  
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
			         G4double maxEnergy) override; 
  
  G4double SampleFluctuations(const G4MaterialCutsCouple*,
                                      const G4DynamicParticle*,
                                      const G4double tcut,
                                      const G4double tmax,
                                      const G4double length,
                                      const G4double meanLoss) override;
  void CorrectionsAlongStep(const G4Material*,
			    const G4ParticleDefinition*,
			    const G4double kinEnergy,
			    const G4double cutEnergy,
			    const G4double& length,
			    G4double& eloss) override;

  G4double Dispersion(const G4Material*, const G4DynamicParticle*,
           const G4double, const G4double, const G4double) override {return 0.;};

  inline G4double ComputeMeanFreePath( const G4ParticleDefinition*,
                                     G4double kineticEnergy,
                                     const G4Material*,
                                     G4double cutEnergy = 0.0,
                                     G4double maxEnergy = DBL_MAX);
  void Initialize();
  void InitRuthELF();
  
  void BuildPhysicsTable(const G4ParticleDefinition* pd);  
  void BuildPrEnergyTable();  
  void BuildElEnergyTable();
  G4double GetPrTransfer( G4double Tkin);
  G4double CorrectPrTransfer( G4double Tkin);
  G4double CorrectElTransfer( G4double Tkin);
  
  G4double GetElTransfer( G4double Tkin);   
  G4double GetPrMFP( G4double Tkin);   
  G4double GetPrdNdx( G4double Tkin);   
  G4double GetElMFP( G4double Tkin);   
  G4double GetEldNdx( G4double Tkin);   
  G4double PrPAId2Ndxdw( G4double omega );  
  G4double ElPAId2Ndxdw( G4double omega );  

  void     SetBe2( G4double be2 ){ fBe2 = be2; };
  G4double GetBe2(){ return fBe2; };
  void     SetOmega( G4double ww ){ fOmega = ww; };
  G4double GetOmega(){ return fOmega;};
  void     SetBias( G4double bb ){ fBias = bb; };
  G4double GetBias(){ return fBias;};

  G4double GetElectronTmax( G4double Tkin );
  G4double GetProtonTmax( G4double Tkin );

  inline G4double GetSumELF( G4double energy );
  inline G4double GetSumRuth( G4double energy );
  
private:

  // Local class members
  G4int fTotBin{0};
  G4int fBinTr{0};
  G4int fBias{0};
  G4double fCof{0.0};
  G4double fBeta{0.0};
  G4double fBe2{0.0};
  G4double fTkin{0.0};

  G4double fBmin{0.0};
  G4double fBmax{0.0};
  G4double fWmin{0.0};
  G4double fWmax{0.0};
  G4double fOmega{0.0};
  G4double fElectronDensity{0.0};
  G4double fNat{0.0};
  G4double fNel{0.0};
  G4double fMass{0.0};

  static const G4int theBin;
  static const G4double theEsum[1007], theELFsum[1007], theRuthSum[1007];
  
  G4DataVector fEsum;
  G4DataVector fELFsum, fRuthSum;
  G4DataVector fPrWmaxVector;
  
  G4Material* fMat{nullptr};
  
  G4PhysicsLogVector*  fBetaVector{nullptr};
  G4PhysicsTable*      fPrEnergyTable{nullptr};
  G4PhysicsTable*      fElEnergyTable{nullptr};
  G4PhysicsLogVector*  fTransferVector{nullptr};

  G4ParticleDefinition* theElectron{nullptr};
  G4ParticleDefinition* theProton{nullptr};
  G4ParticleChangeForLoss* fParticleChange{nullptr};
};    

///////////////////////////////////////////////////////////////////
////////////////  Inline methods //////////////////////////////////
////////////////////////////////////////////////////////////////////


/////////////////// fast (STL) get ELF /////////////////

G4double G4LowPAIH2O::GetSumELF( G4double energy )
{
  G4double ee =  energy; // /CLHEP::eV;
  G4double elf(0.), y1(0.), y2(0.), x1(0.), x2(0.), aa(0.);

  std::size_t nlow = std::lower_bound( fEsum.begin(), fEsum.end(), ee ) - fEsum.begin();

  x1 = fEsum[nlow-1];
  x2 = fEsum[nlow];
  y1 = fELFsum[nlow-1];
  y2 = fELFsum[nlow];
  aa = (y2-y1)/(x2-x1);
  elf = y1 + aa*(ee-x1);  
  return elf;
}

/////////////////////// fast (STL) get Rutherford /////////////////////

G4double G4LowPAIH2O::GetSumRuth( G4double energy )
{
  G4double ee =  energy; // /CLHEP::eV;
  G4double ruth(0.), y1(0.), y2(0.), x1(0.), x2(0.), aa(0.);

  std::size_t nlow = std::lower_bound( fEsum.begin(), fEsum.end(), ee ) - fEsum.begin();

  x1 = fEsum[nlow-1];
  x2 = fEsum[nlow];
  y1 = fRuthSum[nlow-1];
  y2 = fRuthSum[nlow];
  aa = (y2-y1)/(x2-x1);
  ruth = y1 + aa*(ee-x1);
  
  return ruth;
}

///////////////////////////////////////////////

G4double G4LowPAIH2O::ComputeMeanFreePath( const G4ParticleDefinition* pd,
                                      G4double Tkin,
                                      const G4Material*,
					   G4double, // cutEnergy, // = 0.0,
					   G4double ) //maxEnergy ) //  = DBL_MAX)
{
  G4double dndx(0.), mfp(DBL_MAX);
  
  if     ( pd == theProton )    dndx = GetPrdNdx(Tkin);
  else if( pd == theElectron )  dndx = GetEldNdx(Tkin);
  else                          return DBL_MAX;

  if( dndx > 0.) mfp = 1./dndx;
  else           mfp = DBL_MAX;

  mfp *= fBias; 
  return mfp;  
}

#endif   

/////////////////   end of G4LowPAIH2O header file  //////////////////

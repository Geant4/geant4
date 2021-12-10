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
///////////////////////////////////////////////////////////////////////////
//
// base class for 'fast' parametrisation model describing X-ray transition
// created in some G4Envelope. Angular distribuiton is very rough !!! (see DoIt
// method
//
// History:
// 06.10.05  V. Grichine first step to discrete process
// 15.01.02  V. Grichine first version
// 28.07.05, P.Gumplinger add G4ProcessType to constructor
// 28.09.07, V.Ivanchenko general cleanup without change of algorithms
// 19.09.21, V. Grichine, set/get functions for angle anf energy ranges and number of bins

#ifndef G4VXTRenergyLoss_h
#define G4VXTRenergyLoss_h 1

#include "globals.hh"
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicsTable.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VDiscreteProcess.hh"

class G4SandiaTable;
class G4VParticleChange;
class G4PhysicsFreeVector;
class G4PhysicsLinearVector;
class G4PhysicsLogVector;

class G4VXTRenergyLoss : public G4VDiscreteProcess
{
 public:
  explicit G4VXTRenergyLoss(G4LogicalVolume* anEnvelope, G4Material*,
                            G4Material*, G4double, G4double, G4int,
                            const G4String& processName = "XTRenergyLoss",
                            G4ProcessType type          = fElectromagnetic);
  virtual ~G4VXTRenergyLoss();

  virtual void ProcessDescription(std::ostream&) const override;
  virtual void DumpInfo() const override { ProcessDescription(G4cout); };

  G4VXTRenergyLoss(G4VXTRenergyLoss&) = delete;
  G4VXTRenergyLoss& operator=(const G4VXTRenergyLoss& right) = delete;

  // Virtual methods to be implemented in inherited particular TR radiators
  virtual G4double GetStackFactor(G4double energy, G4double gamma,
                                  G4double varAngle);

  virtual G4bool IsApplicable(const G4ParticleDefinition&) override;

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep) override;

  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition) override;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;
  void BuildEnergyTable();
  void BuildAngleForEnergyBank();

  void BuildTable(){};
  void BuildAngleTable();
  void BuildGlobalAngleTable();

  G4complex OneInterfaceXTRdEdx(G4double energy, G4double gamma,
                                G4double varAngle);

  G4double SpectralAngleXTRdEdx(G4double varAngle);

  virtual G4double SpectralXTRdEdx(G4double energy);

  G4double AngleSpectralXTRdEdx(G4double energy);

  G4double AngleXTRdEdx(G4double varAngle);

  G4double OneBoundaryXTRNdensity(G4double energy, G4double gamma,
                                  G4double varAngle) const;

  // for photon energy distribution tables
  G4double XTRNSpectralAngleDensity(G4double varAngle);
  G4double XTRNSpectralDensity(G4double energy);

  // for photon angle distribution tables
  G4double XTRNAngleSpectralDensity(G4double energy);
  G4double XTRNAngleDensity(G4double varAngle);

  void GetNumberOfPhotons();

  // Auxiliary functions for plate/gas material parameters
  G4double GetPlateFormationZone(G4double, G4double, G4double);
  G4complex GetPlateComplexFZ(G4double, G4double, G4double);
  void ComputePlatePhotoAbsCof();
  G4double GetPlateLinearPhotoAbs(G4double);
  void GetPlateZmuProduct();
  G4double GetPlateZmuProduct(G4double, G4double, G4double);

  G4double GetGasFormationZone(G4double, G4double, G4double);
  G4complex GetGasComplexFZ(G4double, G4double, G4double);
  void ComputeGasPhotoAbsCof();
  G4double GetGasLinearPhotoAbs(G4double);
  void GetGasZmuProduct();
  G4double GetGasZmuProduct(G4double, G4double, G4double);

  G4double GetPlateCompton(G4double);
  G4double GetGasCompton(G4double);
  G4double GetComptonPerAtom(G4double, G4double);

  G4double GetXTRrandomEnergy(G4double scaledTkin, G4int iTkin);
  G4double GetXTRenergy(G4int iPlace, G4double position, G4int iTransfer);

  G4double GetRandomAngle(G4double energyXTR, G4int iTkin);
  G4double GetAngleXTR(G4int iTR, G4double position, G4int iAngle);

  // set/get methods for class fields

  void     SetGamma(G4double gamma) { fGamma = gamma; };
  G4double GetGamma() { return fGamma; };
  void     SetEnergy(G4double energy) { fEnergy = energy; };
  G4double GetEnergy() { return fEnergy; };
  void     SetVarAngle(G4double varAngle) { fVarAngle = varAngle; };
  G4double GetVarAngle() { return fVarAngle; };
  void   SetCompton(G4bool pC) { fCompton = pC; };
  G4bool GetCompton() { return fCompton; };



  void     SetAlphaGas(G4double ag){ fAlphaGas = ag;};
  G4double GetAlphaGas() { return fAlphaGas; };  
  void     SetAlphaPlate(G4double ap){ fAlphaPlate = ap;};
  G4double GetAlphaPlate() { return fAlphaPlate; };

  void     SetTheMinEnergyTR(G4double minetr){ fTheMinEnergyTR = minetr;};
  G4double GetTheMinEnergyTR() { return fTheMinEnergyTR; };  
  void     SetTheMaxEnergyTR(G4double maxetr){ fTheMaxEnergyTR = maxetr;};
  G4double GetTheMaxEnergyTR() { return fTheMaxEnergyTR; };  

  void     SetMinEnergyTR(G4double minetr){ fMinEnergyTR = minetr;};
  G4double GetMinEnergyTR() { return fMinEnergyTR; };  
  void     SetMaxEnergyTR(G4double maxetr){ fMaxEnergyTR = maxetr;};
  G4double GetMaxEnergyTR() { return fMaxEnergyTR; };  
  
  void     SetTheMinAngle(G4double minang){ fTheMinAngle = minang;};
  G4double GetTheMinAngle() { return fTheMinAngle; };  
  void     SetTheMaxAngle(G4double maxang){ fTheMaxAngle = maxang;};
  G4double GetTheMaxAngle() { return fTheMaxAngle; };

  void     SetMinThetaTR(G4double minatr){ fMinThetaTR = minatr;};
  G4double GetMinThetaTR() { return fMinThetaTR; };  
  void     SetMaxThetaTR(G4double maxatr){ fMaxThetaTR = maxatr;};
  G4double GetMaxThetaTR() { return fMaxThetaTR; };

  // modes of XTR angle distribution
  
  void   SetFastAngle(G4bool fatr){ fFastAngle = fatr;};
  G4bool GetFastAngle() { return fFastAngle; };  
  void   SetAngleRadDistr(G4bool fatr){ fAngleRadDistr = fatr;};
  G4bool GetAngleRadDistr() { return fAngleRadDistr; };  
  

  

  G4PhysicsLogVector* GetProtonVector() { return fProtonEnergyVector; };
  G4int GetTotBin() { return fTotBin; };
  G4PhysicsFreeVector* GetAngleVector(G4double energy, G4int n);

 protected:
  //   min TR energy
  G4double fTheMinEnergyTR; 
  //   max TR energy
  G4double fTheMaxEnergyTR; 
  G4double fTheMinAngle;  //  min theta of TR quanta
  G4double fTheMaxAngle;  //  1.e-4;  //  max theta of TR quanta

  // static const members
  
  // min Tkin of proton in tables
  static constexpr G4double fMinProtonTkin = 100. * CLHEP::GeV;
  // max Tkin of proton in tables
  static constexpr G4double fMaxProtonTkin = 100. * CLHEP::TeV;
  // physical constants for plasma energy
  static constexpr G4double fPlasmaCof =
    4. * CLHEP::pi * CLHEP::fine_structure_const * CLHEP::hbarc * CLHEP::hbarc *
    CLHEP::hbarc / CLHEP::electron_mass_c2;
  static constexpr G4double fCofTR = CLHEP::fine_structure_const / CLHEP::pi;

  G4int fTotBin;  //  number of bins in log-gamma scale
  G4int fBinTR;   //  number of bins in TR energy-angle vectors
  
  G4ParticleDefinition* fPtrGamma;  // pointer to TR photon

  G4double* fGammaCutInKineticEnergy;  // TR photon cut in energy array
  G4LogicalVolume* fEnvelope;
  G4PhysicsTable* fAngleDistrTable;
  G4PhysicsTable* fEnergyDistrTable;
  G4PhysicsTable* fAngleForEnergyTable;
  G4PhysicsLogVector* fProtonEnergyVector;
  G4PhysicsLogVector* fXTREnergyVector;
  G4SandiaTable* fPlatePhotoAbsCof;
  G4SandiaTable* fGasPhotoAbsCof;

  G4ParticleChange fParticleChange;
  std::vector<G4PhysicsTable*> fAngleBank;

  G4double fGammaTkinCut;  // Tkin cut of TR photon in current mat.
  G4double fMinEnergyTR;   //  min TR energy in material
  G4double fMaxEnergyTR;   //  max TR energy in material
  G4double fMinThetaTR, fMaxThetaTR;    //  min-max theta of TR quanta
  G4double fTotalDist;
  G4double fPlateThick;
  G4double fGasThick;
  G4double fAlphaPlate;
  G4double fAlphaGas;
  G4double fGamma;     // current Lorentz factor
  G4double fEnergy;    // energy and
  G4double fVarAngle;  // angle squared!
  G4double fLambda;
  G4double fSigma1;
  G4double fSigma2;  // plasma energy Sq of matter1/2

  G4int fMatIndex1;
  G4int fMatIndex2;
  G4int fPlateNumber;

  G4bool fExitFlux;
  G4bool fFastAngle, fAngleRadDistr;
  G4bool fCompton;

  G4int secID = -1;  // creator modelID
};

#endif

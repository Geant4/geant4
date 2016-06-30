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
//
// $Id: G4VXTRenergyLoss.hh 97385 2016-06-02 09:59:53Z gcosmo $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// base class for 'fast' parametrisation model describing X-ray transition
// created in some G4Envelope. Anglur distribuiton is very rough !!! (see DoIt
// method
// 
// History:
// 06.10.05 V. Grichine first step to discrete process
// 15.01.02 V. Grichine first version
// 28.07.05, P.Gumplinger add G4ProcessType to constructor
// 28.09.07, V.Ivanchenko general cleanup without change of algorithms
//

#ifndef G4VXTRenergyLoss_h
#define G4VXTRenergyLoss_h 1

#include <complex>
#include "globals.hh"
#include "Randomize.hh"

#include "G4LogicalVolume.hh"

#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Gamma.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VContinuousProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4Integrator.hh"
#include "G4ParticleChange.hh"

class G4SandiaTable;
class G4VParticleChange;
class G4PhysicsFreeVector;
class G4PhysicsLinearVector;

class G4VXTRenergyLoss : public G4VDiscreteProcess  // G4VContinuousProcess
{
public:

  explicit G4VXTRenergyLoss (G4LogicalVolume *anEnvelope,G4Material*,
                    G4Material*, G4double,G4double,G4int,
                    const G4String & processName = "XTRenergyLoss",
                    G4ProcessType type = fElectromagnetic);
  virtual  ~G4VXTRenergyLoss ();

  // These virtual has to be implemented in inherited particular TR radiators
 
  virtual  G4double GetStackFactor( G4double energy, G4double gamma,
                                                     G4double varAngle );

  virtual G4bool IsApplicable(const G4ParticleDefinition&) override;

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
				   const G4Step&  aStep) override;

  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;
  void BuildEnergyTable() ;
  void BuildAngleForEnergyBank() ;

  void BuildTable(){} ;
  void BuildAngleTable() ;
  void BuildGlobalAngleTable() ;

  G4complex OneInterfaceXTRdEdx( G4double energy, 
                                G4double gamma,
                                G4double varAngle ) ;

  G4double SpectralAngleXTRdEdx(G4double varAngle) ;

  virtual  G4double SpectralXTRdEdx(G4double energy) ;

  G4double AngleSpectralXTRdEdx(G4double energy) ;

  G4double AngleXTRdEdx(G4double varAngle) ;


  /////////////////////////////////////////////////////////////

  G4double OneBoundaryXTRNdensity( G4double energy,
                                   G4double gamma,
				   G4double varAngle ) const ;


  // for photon energy distribution tables

  G4double XTRNSpectralAngleDensity(G4double varAngle) ;
  G4double XTRNSpectralDensity(G4double energy) ;
  
  // for photon angle distribution tables

  G4double XTRNAngleSpectralDensity(G4double energy) ;
  G4double XTRNAngleDensity(G4double varAngle) ;

  void GetNumberOfPhotons() ;  

  // Auxiliary functions for plate/gas material parameters

  G4double  GetPlateFormationZone(G4double,G4double,G4double);
  G4complex GetPlateComplexFZ(G4double,G4double,G4double);
  void      ComputePlatePhotoAbsCof();
  G4double  GetPlateLinearPhotoAbs(G4double);
  void      GetPlateZmuProduct() ;
  G4double  GetPlateZmuProduct(G4double,G4double,G4double);

  G4double  GetGasFormationZone(G4double,G4double,G4double);
  G4complex GetGasComplexFZ(G4double,G4double,G4double);
  void      ComputeGasPhotoAbsCof();
  G4double  GetGasLinearPhotoAbs(G4double);
  void      GetGasZmuProduct();
  G4double  GetGasZmuProduct(G4double,G4double,G4double);

  G4double  GetPlateCompton(G4double);
  G4double  GetGasCompton(G4double);
  G4double  GetComptonPerAtom(G4double,G4double);

  G4double  GetXTRrandomEnergy( G4double scaledTkin, G4int iTkin );
  G4double  GetXTRenergy( G4int iPlace, G4double position, G4int iTransfer  );

  G4double  GetRandomAngle( G4double energyXTR, G4int iTkin );
  G4double  GetAngleXTR(G4int iTR,G4double position,G4int iAngle);

  G4double  GetGamma()   {return fGamma;}; 
  G4double  GetEnergy()  {return fEnergy;};                
  G4double  GetVarAngle(){return fVarAngle;};
               
  void SetGamma(G4double gamma)      {fGamma    = gamma;}; 
  void SetEnergy(G4double energy)    {fEnergy   = energy;};                
  void SetVarAngle(G4double varAngle){fVarAngle = varAngle;};               
  void SetAngleRadDistr(G4bool pAngleRadDistr){fAngleRadDistr=pAngleRadDistr;};               
  void SetCompton(G4bool pC){fCompton=pC;};               

  G4PhysicsLogVector* GetProtonVector(){ return fProtonEnergyVector;};
  G4int GetTotBin(){return fTotBin;};           
  G4PhysicsFreeVector* GetAngleVector(G4double energy, G4int n);

protected:

  G4ParticleDefinition* fPtrGamma ;    // pointer to TR photon

  G4double* fGammaCutInKineticEnergy ; // TR photon cut in energy array

  G4double         fGammaTkinCut ;     // Tkin cut of TR photon in current mat.
  G4LogicalVolume* fEnvelope ;
  G4PhysicsTable*  fAngleDistrTable ;
  G4PhysicsTable*  fEnergyDistrTable ;

  G4PhysicsLogVector* fProtonEnergyVector ;
  G4PhysicsLogVector* fXTREnergyVector ;

  G4double fTheMinEnergyTR;            //   min TR energy
  G4double fTheMaxEnergyTR;            //   max TR energy
  G4double fMinEnergyTR;               //  min TR energy in material
  G4double fMaxEnergyTR;               //  max TR energy in material
  G4double fTheMaxAngle;               //  max theta of TR quanta
  G4double fTheMinAngle;               //  max theta of TR quanta
  G4double fMaxThetaTR;                //  max theta of TR quanta
  G4int    fBinTR;                     //  number of bins in TR vectors

  G4double fMinProtonTkin;             // min Tkin of proton in tables
  G4double fMaxProtonTkin;             // max Tkin of proton in tables
  G4int    fTotBin;                    // number of bins in log scale
  G4double fGamma;                     // current Lorentz factor
  G4double fEnergy;                    // energy and
  G4double fVarAngle;                  // angle squared
  G4double fLambda;

  G4double fPlasmaCof ;                // physical consts for plasma energy
  G4double fCofTR ;  

  G4bool   fExitFlux;
  G4bool   fAngleRadDistr;
  G4bool   fCompton;
  G4double fSigma1; 
  G4double fSigma2;                    // plasma energy Sq of matter1/2

  G4int    fMatIndex1;
  G4int    fMatIndex2;
  G4int    fPlateNumber;

  G4double fTotalDist;
  G4double fPlateThick;
  G4double fGasThick;     
  G4double fAlphaPlate;
  G4double fAlphaGas ;

  G4SandiaTable* fPlatePhotoAbsCof;
 
  G4SandiaTable* fGasPhotoAbsCof;

  G4ParticleChange fParticleChange;

  G4PhysicsTable*                    fAngleForEnergyTable;
  std::vector<G4PhysicsTable*>       fAngleBank;

private:

  // copy constructor and hide assignment operator
  G4VXTRenergyLoss(G4VXTRenergyLoss &) = delete;
  G4VXTRenergyLoss & operator=(const G4VXTRenergyLoss &right) = delete;

};

#endif

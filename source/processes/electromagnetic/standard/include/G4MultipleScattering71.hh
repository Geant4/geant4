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
// $Id: G4MultipleScattering71.hh,v 1.6 2008/07/16 11:27:41 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
//------------- G4MultipleScattering71 physics process --------------------------
//               by Laszlo Urban, March 2001
//
// 07-08-01 new methods Store/Retrieve PhysicsTable
// 23-08-01 new angle and z distribution,energy dependence reduced,
//          Store,Retrieve methods commented out temporarily, L.Urban
// 11-09-01 G4MultipleScatteringx put as default: G4MultipleScattering
//          Store,Retrieve methods reactived (mma)
// 13-09-01 Unused TrueToGeomTransformation method deleted,
//          class description (L.Urban)
// 19-09-01 come back to previous process name msc
// 17-04-02 NEW angle distribution + boundary algorithm modified, L.Urban
// 22-04-02 boundary algorithm modified -> important improvement in timing !!!!
//          (L.Urban)
// 24-05-02 changes in data members, L.Urban
// 30-10-02 changes in data members, L.Urban
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 05-02-03 changes in data members, L.Urban
// 28-03-03 Move to model design (V.Ivanchenko)
// 18-04-03 Change name (V.Ivanchenko)
// 16-06-03: ShortLived are not applicable any more (V.Ivanchenko)
// 17-08-04 name of data member facxsi changed to factail together
//          with the corresponding set function (L.Urban)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 15-04-05 optimize internal interfaces (V.Ivanchenko)
// 03-10-05 Process is freezed with the name 71 (V.Ivanchenko)
// 07-03-06 Create G4UrbanMscModel and move there step limit calculation (V.Ivanchenko)
//
//------------------------------------------------------------------------------
//
// $Id: G4MultipleScattering71.hh,v 1.6 2008/07/16 11:27:41 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $

// class description
//
//  The class simulates the multiple scattering for any kind
//  of charged particle.
//
// class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MultipleScattering71_h
#define G4MultipleScattering71_h 1

#include "G4VMultipleScattering.hh"
#include "G4MscModel71.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MultipleScattering71 : public G4VMultipleScattering

{
public:    // with description

  G4MultipleScattering71(const G4String& processName="msc");

  virtual ~G4MultipleScattering71();

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p);

  virtual G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  G4double TruePathLengthLimit(const G4Track&  track,
			       G4double& lambda,
			       G4double  currentMinimalStep);

  // Print few lines of informations about the process: validity range,
  void PrintInfo();

  // geom. step length distribution should be sampled or not
  void Setsamplez(G4bool value);

  // activate boundary algorithm
  void SetBoundary(G4bool value);

  // to reduce the energy/step dependence
  void Setdtrl(G4double value);

  void Setfactail(G4double value);

  // Steplimit after boundary crossing = facrange*range
  // estimated nb of steps at boundary nsmallstep = 1/facrange
  void SetFacrange(G4double val);

  // corrs to transport cross section for high energy
  void SetNuclCorrPar(G4double val);
  void SetFactPar(G4double val);

protected:

  // This function initialise models
  void InitialiseProcess(const G4ParticleDefinition*);

  // This method is used for tracking, it returns step limit
  virtual G4double GetContinuousStepLimit(const G4Track& track,
                                        G4double previousStepSize,
                                        G4double currentMinimalStep,
                                        G4double& currentSafety);

private:

  //  hide assignment operator as  private
  G4MultipleScattering71 & operator = (const G4MultipleScattering71 &right);
  G4MultipleScattering71 ( const G4MultipleScattering71 &);

private:        // data members

  G4MscModel71* model;

  G4double lowKineticEnergy;
  G4double highKineticEnergy;
  G4int    totBins;

  G4double truePathLength;
  G4double geomPathLength;
  G4double trueStepLength;
  G4double range;

  G4double facrange;
  G4double tlimit;
  G4double tlimitmin;
  G4double dtrl;
  G4double NuclCorrPar;
  G4double FactPar;
  G4double factail;
  G4double cf;

  G4int    stepnolastmsc;
  G4int    nsmallstep;

  G4bool   samplez;
  G4bool   boundary;
  G4bool   isInitialized;

};

//--------------------------------------------------------------------
//  inline methods
//--------------------------------------------------------------------

inline G4bool G4MultipleScattering71::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4MultipleScattering71::AlongStepDoIt(
                                                        const G4Track&,
                                                        const G4Step& step)
{
  G4double geomStepLength = step.GetStepLength();
  if((geomStepLength == geomPathLength) && (truePathLength <= range))
     trueStepLength = truePathLength;
  else
     trueStepLength = model->TrueStepLength(geomStepLength);
  fParticleChange.ProposeTrueStepLength(trueStepLength);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4MultipleScattering71::PostStepDoIt(const G4Track& track,
							       const G4Step& step)
{
  fParticleChange.Initialize(track);
  std::vector<G4DynamicParticle*>* p = 0;
  model->SampleSecondaries(p, CurrentMaterialCutsCouple(),track.GetDynamicParticle(),
		    step.GetStepLength(),step.GetPostStepPoint()->GetSafety());
  return &fParticleChange;
}

// geom. step length distribution should be sampled or not
inline void G4MultipleScattering71::Setsamplez(G4bool value)
{
  samplez = value;
}

// activate boundary algorithm
inline void G4MultipleScattering71::SetBoundary(G4bool value)
{
  boundary = value;
}

// to reduce the energy/step dependence
inline void G4MultipleScattering71::Setdtrl(G4double value)                
{
  dtrl = value;
}

inline void G4MultipleScattering71::Setfactail(G4double value)              
{
  factail = value;
}

// Steplimit after boundary crossing = facrange*range
// estimated nb of steps at boundary nsmallstep = 1/facrange
inline void G4MultipleScattering71::SetFacrange(G4double val)              
{
  facrange = val;
  nsmallstep = G4int(std::log((cf+facrange-1.)/facrange)/std::log(cf))+1;
}

// corrs to transport cross section for high energy
inline  void G4MultipleScattering71::SetNuclCorrPar(G4double val)            
{
  NuclCorrPar = val;
}

inline  void G4MultipleScattering71::SetFactPar(G4double val)                
{
  FactPar = val;
}

#endif




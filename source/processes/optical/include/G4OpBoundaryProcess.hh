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
//
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Boundary Process Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpBoundaryProcess.hh
// Description: Discrete Process -- reflection/refraction at
//                                  optical interfaces
// Version:     1.1
// Created:     1997-06-18
// Modified:    2005-07-28 add G4ProcessType to constructor
//              1999-10-29 add method and class descriptors
//              1999-10-10 - Fill NewMomentum/NewPolarization in
//                           DoAbsorption. These members need to be
//                           filled since DoIt calls
//                           aParticleChange.SetMomentumChange etc.
//                           upon return (thanks to: Clark McGrew)
//              2006-11-04 - add capability of calculating the reflectivity
//                           off a metal surface by way of a complex index
//                           of refraction - Thanks to Sehwook Lee and John
//                           Hauptman (Dept. of Physics - Iowa State Univ.)
//              2009-11-10 - add capability of simulating surface reflections
//                           with Look-Up-Tables (LUT) containing measured
//                           optical reflectance for a variety of surface
//                           treatments - Thanks to Martin Janecek and
//                           William Moses (Lawrence Berkeley National Lab.)
//              2013-06-01 - add the capability of simulating the transmission
//                           of a dichronic filter
//              2017-02-24 - add capability of simulating surface reflections
//                           with Look-Up-Tables (LUT) developed in DAVIS
//
// Author:      Peter Gumplinger
//              adopted from work by Werner Keil - April 2/96
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpBoundaryProcess_h
#define G4OpBoundaryProcess_h 1

#include "G4OpticalPhoton.hh"
#include "G4OpticalSurface.hh"
#include "G4RandomTools.hh"
#include "G4VDiscreteProcess.hh"

enum G4OpBoundaryProcessStatus
{
  Undefined,
  Transmission,
  FresnelRefraction,
  FresnelReflection,
  TotalInternalReflection,
  LambertianReflection,
  LobeReflection,
  SpikeReflection,
  BackScattering,
  Absorption,
  Detection,
  NotAtBoundary,
  SameMaterial,
  StepTooSmall,
  NoRINDEX,
  PolishedLumirrorAirReflection,
  PolishedLumirrorGlueReflection,
  PolishedAirReflection,
  PolishedTeflonAirReflection,
  PolishedTiOAirReflection,
  PolishedTyvekAirReflection,
  PolishedVM2000AirReflection,
  PolishedVM2000GlueReflection,
  EtchedLumirrorAirReflection,
  EtchedLumirrorGlueReflection,
  EtchedAirReflection,
  EtchedTeflonAirReflection,
  EtchedTiOAirReflection,
  EtchedTyvekAirReflection,
  EtchedVM2000AirReflection,
  EtchedVM2000GlueReflection,
  GroundLumirrorAirReflection,
  GroundLumirrorGlueReflection,
  GroundAirReflection,
  GroundTeflonAirReflection,
  GroundTiOAirReflection,
  GroundTyvekAirReflection,
  GroundVM2000AirReflection,
  GroundVM2000GlueReflection,
  Dichroic,
  CoatedDielectricReflection,
  CoatedDielectricRefraction,
  CoatedDielectricFrustratedTransmission
};

class G4OpBoundaryProcess : public G4VDiscreteProcess
{
 public:
  explicit G4OpBoundaryProcess(const G4String& processName = "OpBoundary",
                               G4ProcessType type          = fOptical);
  virtual ~G4OpBoundaryProcess();

  virtual G4bool IsApplicable(
    const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable' only for an optical photon.

  virtual G4double GetMeanFreePath(const G4Track&, G4double,
                                   G4ForceCondition* condition) override;
  // Returns infinity; i. e. the process does not limit the step, but sets the
  // 'Forced' condition for the DoIt to be invoked at every step. However, only
  // at a boundary will any action be taken.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;
  // This is the method implementing boundary processes.

  virtual G4OpBoundaryProcessStatus GetStatus() const;
  // Returns the current status.

  virtual void SetInvokeSD(G4bool);
  // Set flag for call to InvokeSD method.

  virtual void PreparePhysicsTable(const G4ParticleDefinition&) override;

  virtual void Initialise();

  void SetVerboseLevel(G4int);

 private:
  G4OpBoundaryProcess(const G4OpBoundaryProcess& right) = delete;
  G4OpBoundaryProcess& operator=(const G4OpBoundaryProcess& right) = delete;

  G4bool G4BooleanRand(const G4double prob) const;

  G4ThreeVector GetFacetNormal(const G4ThreeVector& Momentum,
                               const G4ThreeVector& Normal) const;

  void DielectricMetal();
  void DielectricDielectric();

  void DielectricLUT();
  void DielectricLUTDAVIS();

  void DielectricDichroic();
  void CoatedDielectricDielectric();

  void ChooseReflection();
  void DoAbsorption();
  void DoReflection();

  G4double GetIncidentAngle();
  // Returns the incident angle of optical photon

  G4double GetReflectivity(G4double E1_perp, G4double E1_parl,
                           G4double incidentangle, G4double RealRindex,
                           G4double ImaginaryRindex);
  // Returns the Reflectivity on a metallic surface

  G4double GetReflectivityThroughThinLayer(G4double sinTL, G4double E1_perp,
                                           G4double E1_parl, G4double wavelength,
                                           G4double cost1, G4double cost2);
  // Returns the Reflectivity on a coated surface

  void CalculateReflectivity();

  void BoundaryProcessVerbose() const;

  // Invoke SD for post step point if the photon is 'detected'
  G4bool InvokeSD(const G4Step* step);

  G4ThreeVector fOldMomentum;
  G4ThreeVector fOldPolarization;

  G4ThreeVector fNewMomentum;
  G4ThreeVector fNewPolarization;

  G4ThreeVector fGlobalNormal;
  G4ThreeVector fFacetNormal;

  const G4Material* fMaterial1;
  const G4Material* fMaterial2;

  G4OpticalSurface* fOpticalSurface;

  G4MaterialPropertyVector* fRealRIndexMPV;
  G4MaterialPropertyVector* fImagRIndexMPV;
  G4Physics2DVector* fDichroicVector;

  G4double fPhotonMomentum;
  G4double fRindex1;
  G4double fRindex2;

  G4double fSint1;

  G4double fReflectivity;
  G4double fEfficiency;
  G4double fTransmittance;
  G4double fSurfaceRoughness;

  G4double fProb_sl, fProb_ss, fProb_bs;
  G4double fCarTolerance;

  // Used by CoatedDielectricDielectric()
  G4double fCoatedRindex, fCoatedThickness;

  G4OpBoundaryProcessStatus fStatus;
  G4OpticalSurfaceModel fModel;
  G4OpticalSurfaceFinish fFinish;

  G4int f_iTE, f_iTM;

  G4int fNumWarnings; // number of times small step warning printed

  size_t idx_dichroicX      = 0;
  size_t idx_dichroicY      = 0;
  size_t idx_rindex1        = 0;
  size_t idx_rindex_surface = 0;
  size_t idx_reflect        = 0;
  size_t idx_eff            = 0;
  size_t idx_trans          = 0;
  size_t idx_lobe           = 0;
  size_t idx_spike          = 0;
  size_t idx_back           = 0;
  size_t idx_rindex2        = 0;
  size_t idx_groupvel       = 0;
  size_t idx_rrindex        = 0;
  size_t idx_irindex        = 0;
  size_t idx_coatedrindex   = 0;

  // Used by CoatedDielectricDielectric()
  G4bool fCoatedFrustratedTransmission = true;

  G4bool fInvokeSD;
};

////////////////////
// Inline methods
////////////////////

inline G4bool G4OpBoundaryProcess::G4BooleanRand(const G4double prob) const
{
  // Returns a random boolean variable with the specified probability
  return (G4UniformRand() < prob);
}

inline G4bool G4OpBoundaryProcess::IsApplicable(
  const G4ParticleDefinition& aParticleType)
{
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

inline G4OpBoundaryProcessStatus G4OpBoundaryProcess::GetStatus() const
{
  return fStatus;
}

inline void G4OpBoundaryProcess::ChooseReflection()
{
  G4double rand = G4UniformRand();
  if(rand < fProb_ss)
  {
    fStatus      = SpikeReflection;
    fFacetNormal = fGlobalNormal;
  }
  else if(rand < fProb_ss + fProb_sl)
  {
    fStatus = LobeReflection;
  }
  else if(rand < fProb_ss + fProb_sl + fProb_bs)
  {
    fStatus = BackScattering;
  }
  else
  {
    fStatus = LambertianReflection;
  }
}

inline void G4OpBoundaryProcess::DoAbsorption()
{
  fStatus = Absorption;

  if(G4BooleanRand(fEfficiency))
  {
    // EnergyDeposited =/= 0 means: photon has been detected
    fStatus = Detection;
    aParticleChange.ProposeLocalEnergyDeposit(fPhotonMomentum);
  }
  else
  {
    aParticleChange.ProposeLocalEnergyDeposit(0.0);
  }

  fNewMomentum     = fOldMomentum;
  fNewPolarization = fOldPolarization;

  aParticleChange.ProposeTrackStatus(fStopAndKill);
}

inline void G4OpBoundaryProcess::DoReflection()
{
  if(fStatus == LambertianReflection)
  {
    fNewMomentum = G4LambertianRand(fGlobalNormal);
    fFacetNormal = (fNewMomentum - fOldMomentum).unit();
  }
  else if(fFinish == ground)
  {
    fStatus = LobeReflection;
    if(!fRealRIndexMPV || !fImagRIndexMPV)
    {
      fFacetNormal = GetFacetNormal(fOldMomentum, fGlobalNormal);
    }
    // else
      // complex ref. index to be implemented
    fNewMomentum =
      fOldMomentum - (2. * fOldMomentum * fFacetNormal * fFacetNormal);
  }
  else
  {
    fStatus      = SpikeReflection;
    fFacetNormal = fGlobalNormal;
    fNewMomentum =
      fOldMomentum - (2. * fOldMomentum * fFacetNormal * fFacetNormal);
  }
  fNewPolarization =
    -fOldPolarization + (2. * fOldPolarization * fFacetNormal * fFacetNormal);
}

#endif /* G4OpBoundaryProcess_h */

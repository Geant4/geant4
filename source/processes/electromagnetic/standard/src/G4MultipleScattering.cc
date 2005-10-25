//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4MultipleScattering.cc,v 1.37 2005-10-25 18:31:28 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------------------
// 16/05/01 value of cparm changed , L.Urban
// 18/05/01 V.Ivanchenko Clean up against Linux ANSI compilation
// 07/08/01 new methods Store/Retrieve PhysicsTable (mma)
// 23-08-01 new angle and z distribution,energy dependence reduced,
//          Store,Retrieve methods commented out temporarily, L.Urban
// 27-08-01 in BuildPhysicsTable:aParticleType.GetParticleName()=="mu+" (mma)
// 28-08-01 GetContinuousStepLimit and AlongStepDoIt moved from .icc file (mma)
// 03-09-01 value of data member factlim changed, L.Urban
// 10-09-01 small change in GetContinuousStepLimit, L.Urban
// 11-09-01 G4MultipleScatteringx put as default G4MultipleScattering
//          store/retrieve physics table reactivated (mma)
// 13-09-01 corr. in ComputeTransportCrossSection, L.Urban
// 14-09-01 protection in GetContinuousStepLimit, L.Urban
// 17-09-01 migration of Materials to pure STL (mma)
// 27-09-01 value of data member factlim changed, L.Urban
// 31-10-01 big fixed in PostStepDoIt,L.Urban
// 17-04-02 NEW angle distribution + boundary algorithm modified, L.Urban
// 22-04-02 boundary algorithm modified -> important improvement in timing 
// 24-04-02 some minor changes in boundary algorithm, L.Urban
// 06-05-02 bug fixed in GetContinuousStepLimit, L.Urban
// 24-05-02 changes in angle distribution and boundary algorithm, L.Urban
// 11-06-02 bug fixed in ComputeTransportCrossSection, L.Urban
// 12-08-02 bug fixed in PostStepDoIt (lateral displacement), L.Urban
// 15-08-02 new angle distribution, L.Urban
// 26-09-02 angle distribution + boundary algorithm modified, L.Urban
// 15-10-02 temporary fix for proton scattering
// 30-10-02 modified angle distribution,mods in boundary algorithm,
//          changes in data members, L.Urban
// 11-12-02 precision problem in ComputeTransportCrossSection
//          for small Tkin/for heavy particles cured from L.Urban
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 05-02-03 changes in data members, new sampling for geom.
//          path length, step dependence reduced with new
//          method
// 28-03-03 Move to model design (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 23-04-04 value of data member dtrl changed from 0.15 to 0.05 (L.Urban)
// 17-08-04 name of facxsi changed to factail (L.Urban)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 07-02-05 correction in order to have a working Setsamplez function (L.Urban)
// 15-04-05 optimize internal interface (V.Ivanchenko)
// 12-09-05 new TruePathLengthLimit - facrange works for every track from
//             start, geometry also influences the limit
// 02-10-05 conditions limiting the step are finalized + code cleaning (L.Urban)
// 03-10-05 weaker step limitation for Tkin > Tlimit (L.Urban)
// 05-10-05 value of data member tlimitmin has been changed (L.Urban)
// 06-10-05 correction in TruePathLengthLimit, timing improved.(L.Urban)
// 07-10-05 bug fixed in TruePathLengthLimit (L.Urban)
// 11-10-05 change in TruePathLengthLimit conditions,slightly better 
//          timing and much weaker cut dependence (L.Urban)
// 13-10-05 move SetFacrange(0.02) from InitialiseProcess to constructor
// 23-10-05 new Boolean data member prec (false ~ 7.1 like, true new step
//          limit in TruePathLengthLimit, L.Urban)
// 25-10-05 prec renamed to stepLimitAlgorithm, set function triggers
//          'default' facrange too, true - 0.02, false - 0.2 (L.Urban)
//
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MultipleScattering.hh"
#include "G4MscModel.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MultipleScattering::G4MultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName),
    totBins(120),
    dtrl(0.05),
    factail(1.0),
    boundary(true),
    isInitialized(false)
{
  lowKineticEnergy = 0.1*keV;
  highKineticEnergy= 100.*TeV;

  SetstepLimitAlgorithm(true);
  SetFacrange(0.02);

  Tkinlimit        = 2.*MeV;
  tlimit           = 1.e10*mm;
  tlimitmin        = facrange*1.e-3*mm;
  geombig          = 1.e50*mm;
  geommin          = 5.e-6*mm;
  facgeom          = 4.;
  safety           = 0.*mm;
  facsafety        = 0.95;
  facsafety2       = 0.10;

  SetBinning(totBins);
  SetMinKinEnergy(lowKineticEnergy);
  SetMaxKinEnergy(highKineticEnergy);

  Setsamplez(false) ;
  SetLateralDisplasmentFlag(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MultipleScattering::~G4MultipleScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4MultipleScattering::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MultipleScattering::InitialiseProcess
                                         (const G4ParticleDefinition* particle)
{
  if(isInitialized) return;

  if (particle->GetParticleType() == "nucleus") {
    boundary = false;
    SetLateralDisplasmentFlag(false);
    SetBuildLambdaTable(false);
  } else {
    SetBuildLambdaTable(true);
  }
  // compute Tlimit for particle
  Tlimit = Tkinlimit*electron_mass_c2/particle->GetPDGMass();
  G4MscModel* em = new G4MscModel(dtrl,factail,samplez);
  em->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
  em->SetLowEnergyLimit(lowKineticEnergy);
  em->SetHighEnergyLimit(highKineticEnergy);
  AddEmModel(1, em);
  isInitialized = true;
  navigator = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MultipleScattering::TruePathLengthLimit(const G4Track&  track,
                                                   G4double& lambda,
                                                   G4double  currentMinimalStep)
{
  G4double tPathLength = currentMinimalStep;
  G4double range = CurrentRange() ;

  //standard  version
  if(stepLimitAlgorithm)
  {
    safety = track.GetStep()->GetPreStepPoint()->GetSafety() ;
    if((track.GetStep()->GetPreStepPoint()->GetStepStatus() ==
                                                 fGeomBoundary)
       || (track.GetCurrentStepNumber() == 1))
    {
      // not so strong step restriction above Tlimit
      G4double facr = facrange;
      if(track.GetKineticEnergy() > Tlimit)
        facr *= track.GetKineticEnergy()/Tlimit;

      // constraint from the physics
      if (range > lambda) tlimit = facr*range;
      else                tlimit = facr*lambda;

      G4double geomlimit = GeomLimit(track);
      // constraint from the geometry (if tlimit above is too big)
      if ((geomlimit > geommin) && (tlimit > geomlimit/facgeom))
      {
        if(track.GetStep()->GetPreStepPoint()->GetStepStatus() ==
                                                    fGeomBoundary)
          tlimit = geomlimit/facgeom;
        else
          tlimit = 2.*geomlimit/facgeom;
      }

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;
    }
    if(track.GetCurrentStepNumber() == 1)
    {
      // range <= safety ---> particle is not able to leave volume
      if(range <= facsafety*safety)
        return currentMinimalStep;

      //if track starts far from boundaries increase tlimit!
      if(tlimit < facsafety2*safety)
        tlimit = facsafety2*safety ;
    }
    // range <= safety ---> particle is not able to leave volume
    if(range <= facsafety*safety)
      return currentMinimalStep;

    if(tPathLength > tlimit) tPathLength = tlimit;

    if((tPathLength < facsafety2*safety)
         && (currentMinimalStep > facsafety2*safety))
      tPathLength = facsafety2*safety;
  }
  // version similar to 7.1 (needed for some experiments)
  else
  {
    if(track.GetCurrentStepNumber() == 1)
      tlimit = geombig;
    if (track.GetStep()->GetPreStepPoint()->GetStepStatus() ==
                                                 fGeomBoundary)
    {
      if (range > lambda) tlimit = facrange*range;
      else                tlimit = facrange*lambda;

      if(tlimit < tlimitmin) tlimit = tlimitmin;
    }
    if(tPathLength > tlimit) tPathLength = tlimit;
  }

  return tPathLength ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MultipleScattering::GeomLimit(const G4Track&  track)
{
  G4double geomlimit = geombig;

  // do not call navigator for the World volume
  if((track.GetVolume() != 0)
           && (track.GetVolume()->GetName() != "World"))
  {
    const G4double cstep = geombig;
    navigator->LocateGlobalPointWithinVolume(
                  track.GetStep()->GetPreStepPoint()->GetPosition());
    geomlimit = navigator->ComputeStep(
                  track.GetStep()->GetPreStepPoint()->GetPosition(),
                  track.GetMomentumDirection(),
                  cstep,
                  safety);

    if(geomlimit < geommin) geomlimit = geommin;
  }

  return geomlimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MultipleScattering::PrintInfo()
{
  if(boundary) {
    G4cout << "      Boundary/stepping algorithm is active with facrange= "
           << facrange
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


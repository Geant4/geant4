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
// $Id: G4MultipleScattering80.cc,v 1.2 2005-09-05 12:13:46 urban Exp $
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
// 11-09-01 G4MultipleScattering80x put as default G4MultipleScattering80
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
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MultipleScattering80.hh"
#include "G4MscModel80.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MultipleScattering80::G4MultipleScattering80(const G4String& processName)
  : G4VMultipleScattering(processName),
    totBins(120),
    facrange(0.02),
    dtrl(0.05),
    NuclCorrPar (0.0615),
    FactPar(0.40),
    factail(1.0),
    boundary(true),
    isInitialized(false)
{
  lowKineticEnergy = 0.1*keV;
  highKineticEnergy= 100.*TeV;

  tlimit           = 1.e10*mm;
  tlimitmin        = facrange*1.e-3*mm;
  factlim          = 0.5;
  geombig          = 1.e50*mm;
  geommin          = 5.e-7*mm;
  facgeom          = 2.;
  tskin            = tlimit/facgeom;
  tid              = -1;
  pid              = -1;
  stepnobound      = 100000000;
  nsmallstep       = G4int(facgeom);
  safety           = 0.*mm;

  SetBinning(totBins);
  SetMinKinEnergy(lowKineticEnergy);
  SetMaxKinEnergy(highKineticEnergy);

  Setsamplez(false) ;
  SetLateralDisplasmentFlag(true);

  ipr = 0;
  npr = 0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MultipleScattering80::~G4MultipleScattering80()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4MultipleScattering80::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MultipleScattering80::InitialiseProcess
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
  G4MscModel80* em = new G4MscModel80(dtrl,NuclCorrPar,FactPar,factail,samplez);
  em->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
  em->SetLowEnergyLimit(lowKineticEnergy);
  em->SetHighEnergyLimit(highKineticEnergy);
  AddEmModel(1, em);
  isInitialized = true;
  navigator = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();
  SetFacrange(0.02); // **********************************************
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MultipleScattering80::TruePathLengthLimit(const G4Track&  track,
                                                   G4double& lambda,
                                                   G4double  currentMinimalStep)
{
  G4double tPathLength = currentMinimalStep;
  G4double range = CurrentRange() ;

  if((track.GetCurrentStepNumber() == 1) ||
     (track.GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary))
  {
    // constraint from the physics
    if (range > lambda) tlimit = facrange*range;
    else                tlimit = facrange*lambda;

    //lower limit for tlimit
    if(tlimit < tlimitmin) tlimit = tlimitmin;
    tskin = tlimit/facgeom;

    if(track.GetCurrentStepNumber() == 1)
    {
      tid = track.GetTrackID() ;
      pid = track.GetParentID() ;
      stepnobound      = 100000000;
    }
    if(track.GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary)
    {
      stepnobound = track.GetCurrentStepNumber() ;
    }
  }
  if(ipr < npr)
  {
    ipr += 1;
    G4cout.precision(6);
    G4cout << "\nparticle=" << track.GetDefinition()->GetParticleName() 
           << G4endl;
    G4cout << "tID=" << track.GetTrackID() << "  pID=" << track.GetParentID()
           << "  stepno=" << track.GetCurrentStepNumber() << G4endl;
    G4cout << "tid=" << tid << "  pid=" << pid << "  stepnobound="
           << stepnobound << G4endl;
    G4cout << "material=" << track.GetMaterial()->GetName()
           << "  Tkin=" << track.GetKineticEnergy() << "  range=" << range
           << "  lambda=" << lambda << G4endl;
    G4cout << "facrange=" << facrange << "  tlimit=" << tlimit
           << "  tskin=" << tskin << G4endl;
    G4cout << "x,y,z:" << track.GetStep()->GetPreStepPoint()->GetPosition()
           << G4endl;
    G4cout << "dir  :" << track.GetMomentumDirection() << G4endl;
  }

  G4double geomlimit = GeomLimit(track);
  if(ipr < npr) G4cout << "safety=" << safety << G4endl;
  if(ipr < npr) G4cout << "geomlimit=" << geomlimit << "  currentMinimalStep="
                       << currentMinimalStep << G4endl;

  // simple case: particle is not able to reach the boundary
  if(safety >= range)
    tPathLength = range ;
  else
  {
    // small steps just after crossing a boundary
    if((track.GetCurrentStepNumber() >= stepnobound) &&
       (track.GetCurrentStepNumber() <= stepnobound+nsmallstep) &&
       (track.GetTrackID() == tid) && (track.GetParentID() == pid))
    {
      if(track.GetCurrentStepNumber() == stepnobound)
        tPathLength = geommin ;
      else
        if(tPathLength > tskin) tPathLength = tskin;
    }
    else
    {
      if(tPathLength > tlimit) tPathLength = tlimit;
      if((factlim*safety > tlimit) && (currentMinimalStep > factlim*safety))
        tPathLength = factlim*safety;

    }
    //check geometry as well (small steps before reaching a boundary)
    if(tPathLength > geomlimit) tPathLength = geomlimit;
  }

  if(ipr < npr)
  {
    G4cout << "before return tPathLength=" << tPathLength
           << " ---------------------------------------------- " << G4endl;
  }

  return tPathLength ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MultipleScattering80::GeomLimit(const G4Track&  track)
{
  G4double geomlimit = geombig;
  safety = track.GetStep()->GetPreStepPoint()->GetSafety() ;
  // do not call navigator for big geommin and for World
  if((geommin < geombig) && (track.GetVolume() != 0)
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

    if(geomlimit > facgeom*tskin) geomlimit -= facgeom*tskin;
    else if(geomlimit > tskin) geomlimit = tskin;
    if(geomlimit < geommin) geomlimit = geommin;
  }

  return geomlimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MultipleScattering80::PrintInfo()
{
  if(boundary) {
    G4cout << "      Boundary algorithm is active with facrange= "
           << facrange
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


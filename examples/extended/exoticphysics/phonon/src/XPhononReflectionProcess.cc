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
/// \file exoticphysics/phonon/src/XPhononReflectionProcess.cc
/// \brief Implementation of the XPhononReflectionProcess class
//
// $Id$
//

#include "XPhononReflectionProcess.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4RandomTools.hh"

#include "XTPhononFast.hh"
#include "XTPhononSlow.hh"
#include "XLPhonon.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"

#include "XPhononTrackInformation.hh"
#include "XLatticeManager3.hh"

#include "G4SystemOfUnits.hh"


XPhononReflectionProcess::XPhononReflectionProcess(const G4String& aName)
:G4VDiscreteProcess(aName)
{
   fAlminum = NULL;
   kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
   G4cout<<"\n XPhononReflectionProcess::Constructor: Geometry surface tolerance is: " << kCarTolerance /mm << " mm";
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhononReflectionProcess::~XPhononReflectionProcess()
{;}

XPhononReflectionProcess::XPhononReflectionProcess(XPhononReflectionProcess& right)
: G4VDiscreteProcess(right)
{;}
 
G4double 
  XPhononReflectionProcess::GetMeanFreePath( 
       const G4Track&, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{
// Always return DBL_MAX and Forced
// This ensures that the process is called
// at the end of every step. In 
// PostStepDoIt the process decides whether
// the step encountered a volume boundary 
// and a reflection should be applied

   *condition = Forced;

   return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



G4VParticleChange*
  XPhononReflectionProcess::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step& aStep )
{ 

  //This process handles the interaction of phonons with
  //boundaries. Implementation of this class is highly 
  //geometry dependent.Currently, phonons are killed when
  //they reach a boundary. If the other side of the 
  //boundary was Al, a hit is registered.
  
  aParticleChange.Initialize(aTrack);
   
  //Check if current step is limited by a volume boundary
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if(postStepPoint->GetStepStatus()!=fGeomBoundary)
   {

     //make sure that correct phonon velocity is used after the step

     XLatticeManager3* LM = XLatticeManager3::GetXLatticeManager();
     XPhysicalLattice* Lattice = LM->GetXPhysicalLattice(aTrack.GetVolume());

     int pol = 0;
     if (aTrack.GetDefinition() == XLPhonon::PhononDefinition()){
       pol=0;
     }
     else if (aTrack.GetDefinition() == XTPhononSlow::Definition()){
       pol=1;
     }
     else if (aTrack.GetDefinition() == XTPhononFast::Definition()){
       pol=2;
     }      
     
     //Since step was not a volume boundary, just set correct phonon velocity and return
     aParticleChange.ProposeVelocity(Lattice->MapKtoV(pol, aTrack.GetMomentumDirection())*m/s);
     return &aParticleChange;
   }
  
 
 //do nothing but return is the step is too short
 //This is to allow actual reflection where after
 //the first boundary crossing a second, infinitesimal
 //step occurs crossing back into the original volume
 if(aTrack.GetStepLength()<=kCarTolerance/2)
   { 

     return &aParticleChange;
   }
  
  G4double eKin = aTrack.GetKineticEnergy();     
  aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  return &aParticleChange; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool XPhononReflectionProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
  return ((&aPD==XTPhononFast::PhononDefinition())||(&aPD==XLPhonon::PhononDefinition())||(&aPD==XTPhononSlow::PhononDefinition()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhononReflectionProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{
   if(!fAlminum)
   { fAlminum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



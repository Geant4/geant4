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
/// \file exoticphysics/phonon/src/XPhononScatteringProcess.cc
/// \brief Implementation of the XPhononScatteringProcess class
//
// $Id$
//

#include "XPhononScatteringProcess.hh"

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "XLPhonon.hh"
#include "XTPhononFast.hh"
#include "XTPhononSlow.hh"
#include "XPhononTrackInformation.hh"
#include "XPhysicalLattice.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"


XPhononScatteringProcess::XPhononScatteringProcess(const G4String& aName)
: G4VDiscreteProcess(aName)
{

   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhononScatteringProcess::~XPhononScatteringProcess()
{;}

XPhononScatteringProcess::XPhononScatteringProcess(XPhononScatteringProcess& right)
: G4VDiscreteProcess(right)
{;}
 
G4double 
  XPhononScatteringProcess::GetMeanFreePath( 
       const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{

  //Use LatticeManager3 to get physical lattice of current volume
  XLatticeManager3* LM = XLatticeManager3::GetXLatticeManager();
  XPhysicalLattice* Lattice = LM->GetXPhysicalLattice(aTrack.GetVolume());
  if(Lattice==0) G4cout<<"\n\nXPhononScatteringProcess::GetMeanFreePath: WARNING!! PHYSICAL LATTICE POINTER IS NULL!!!\n\n";

  //Dynamical constants retrieved from PhysicalLattice
  G4double B=Lattice->GetScatteringConstant();
  G4double h=6.626e-34*m2*kg/s;
  G4double E= aTrack.GetKineticEnergy();

  //Calculate mean free path
  G4double mfp = 1/((E/h)*(E/h)*(E/h)*(E/h)*B)*aTrack.GetVelocity();

   *condition = NotForced;
 
   return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4VParticleChange*
  XPhononScatteringProcess::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step& aStep)
{

  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if(postStepPoint->GetStepStatus()==fGeomBoundary)
   { return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);}

    //Initialize particle change
    aParticleChange.Initialize(aTrack);

    // Kill current track, then create scattered
    // phonon as a secondary
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.SetNumberOfSecondaries(1);

    //randomly generate a new direction
    //modeMixer determines what the new 
    //polarization type will be
    G4Track* sec;
    G4ThreeVector vgroup;  
    G4ThreeVector newDir = G4RandomDirection();
    G4double modeMixer = G4UniformRand();
    
    //Use LattaiceManager3 to obtains PhysicalLattice of current volume
    XLatticeManager3* LM = XLatticeManager3::GetXLatticeManager();
    XPhysicalLattice* Lattice = LM->GetXPhysicalLattice(aTrack.GetVolume());
    double cProbST=Lattice->GetSTDOS();
    double cProbFT=Lattice->GetFTDOS()+cProbST;

    //Generate the new track after scattering
    //the probabilities for the different po-
    //larization types depends on the DOS
    int polarization;
    if(modeMixer<cProbST){  
      polarization = 1;
      vgroup=Lattice->MapKtoVDir(1, newDir);
      vgroup=Lattice->fLocalToGlobal.TransformAxis(vgroup);
      sec=new G4Track(new G4DynamicParticle(XTPhononSlow::PhononDefinition(),vgroup.unit(), aTrack.GetKineticEnergy()), aTrack.GetGlobalTime(), aTrack.GetPosition());

    }else if(modeMixer<cProbFT){
      polarization = 2;
      vgroup=Lattice->MapKtoVDir(2, newDir);
      vgroup=Lattice->fLocalToGlobal.TransformAxis(vgroup);
      sec=new G4Track(new G4DynamicParticle(XTPhononFast::PhononDefinition(),vgroup.unit(), aTrack.GetKineticEnergy()), aTrack.GetGlobalTime(), aTrack.GetPosition());

    } else {
      polarization = 0;
      vgroup=Lattice->MapKtoVDir(0, newDir);
      vgroup=Lattice->fLocalToGlobal.TransformAxis(vgroup);
      sec=new G4Track(new G4DynamicParticle(XLPhonon::PhononDefinition(),vgroup.unit(), aTrack.GetKineticEnergy()), aTrack.GetGlobalTime(), aTrack.GetPosition());

    }

    sec->SetUserInformation(new XPhononTrackInformation(Lattice->fLocalToGlobal.TransformAxis(newDir)));
    sec->SetVelocity(Lattice->MapKtoV(polarization, newDir)*m/s);    
    sec->UseGivenVelocity(true);
    aParticleChange.AddSecondary(sec);

    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4bool XPhononScatteringProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
  return ((&aPD==XLPhonon::PhononDefinition())|(&aPD==XTPhononFast::PhononDefinition())|(&aPD==XTPhononSlow::PhononDefinition()));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


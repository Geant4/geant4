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
// Rich advanced example for Geant4
// RichTbSteppingAction.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "globals.hh"
#include "RichTbSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "RichTbAnalysisManager.hh"
#include "RichTbMaterial.hh"
#include "RichTbGeometryParameters.hh"
#include "RichTbMaterialParameters.hh"
#include "RichTbRunConfig.hh"
#include "RichTbPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4ThreeVector.hh"
#include "G4OpticalPhoton.hh"
#include "G4PionMinus.hh"

#ifdef G4ANALYSIS_USE
  #include "AIDA/AIDA.h"
#endif

RichTbSteppingAction::
RichTbSteppingAction(RichTbRunConfig* rConfig ,
                     RichTbPrimaryGeneratorAction* RPrimGenAction)
{
  richtbRunConfig= rConfig;
  rPrimGenAction = RPrimGenAction;
  HpdPhElectronKE=rConfig->getHpdPhElectronEnergy();
  uParticleChange=new G4VParticleChange();  
}

RichTbSteppingAction::~RichTbSteppingAction()
{
}

void RichTbSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  RichTbGenericHisto(aStep);
}

void RichTbSteppingAction:: RichTbDebugHisto(const G4Step*)
{
}

void RichTbSteppingAction::RichTbGenericHisto(const G4Step* aStep)
{
   G4StepPoint* pPreStepPoint  = aStep ->GetPreStepPoint();
   G4StepPoint* pPostStepPoint = aStep ->GetPostStepPoint();
   const G4ThreeVector prePos= pPreStepPoint->GetPosition();
   const G4ThreeVector postPos= pPostStepPoint->GetPosition();

   // In the following 1000 mm is the Z coord of a point
   // between the mirror and the dowstream end of the vessel.

   if(prePos.z()<1000*mm && prePos.z() > 0.0*mm )
   {
     //check to see if we are at a boundary

     if (pPostStepPoint->GetStepStatus() == fGeomBoundary)
     {

        G4Track* aPhotTrack = aStep -> GetTrack();
        const G4DynamicParticle* aParticle = aPhotTrack->GetDynamicParticle();
        const G4double PhotonEnergy = aParticle->GetKineticEnergy();

        G4String VolNameD="Agel";
        G4String VolNameE="VesselEnclosure";
        G4String VolNameF="MirrorSphe";
        G4String VolNameG="RadFrame";
        G4String VolNameH="FilterBox";
        G4String VolNameQ="HpdQuartzWindow";
        G4String VolNameP="HpdMaster";
     
        if( aParticle->GetDefinition() == G4OpticalPhoton::OpticalPhoton())
	{
	   if( pPreStepPoint -> GetPhysicalVolume() && 
		pPostStepPoint -> GetPhysicalVolume() )
	   {
	     G4String tpreVol = pPreStepPoint->GetPhysicalVolume()->GetName();
	     G4String tpostVol = pPostStepPoint->GetPhysicalVolume()->GetName();

#ifdef G4ANALISYS_USE
	     if(richtbRunConfig-> GetRichTbParticleEnergyCode() == 1 ||
	       richtbRunConfig-> GetRichTbParticleEnergyCode() == 2 )
	     {
               // Optical photons are generated as beam particle.
               // count the photons entering the aerogel volume.
               // this is essentially same as the photons generated.
	       
               if(( tpreVol ==  VolNameG || tpreVol == VolNameE )  && 
                  ( tpostVol ==  VolNameD ) )
	       {
                  RichTbAnalysisManager * analysis =
		    RichTbAnalysisManager::getInstance();
                  analysis->bumpNumPhotonsBeforeAerogel();
               }
             }
#endif

      if(PhotonEnergy > 0.0 )
      {

        const G4double PhotonWavelength =  
                       PhotMomWaveConv*1.0*eV/ PhotonEnergy;    

        G4String tpreVol = pPreStepPoint -> GetPhysicalVolume()->GetName();
        G4String tpostVol = pPostStepPoint -> GetPhysicalVolume()->GetName();

        // First for the mirror volume

#ifdef G4ANALYSIS_USE
	RichTbAnalysisManager * analysis =
	   RichTbAnalysisManager::getInstance();
	if(tpreVol ==  VolNameE && tpostVol ==  VolNameF )
	{
	  analysis->getfhistoWBeforeMirror()->fill(PhotonWavelength);
	  analysis->bumpNumPhotonsBeforeMirror();
	}
	if(tpreVol ==  VolNameF && tpostVol ==  VolNameE )
	{
          analysis->getfhistoWAfterMirror()->fill(PhotonWavelength);
          analysis->bumpNumPhotonsAfterMirror();
	}

        //Now for the Aerogel Volume
        if(richtbRunConfig-> GetRichTbParticleEnergyCode() == 0 )
	{
          if(( tpreVol ==  VolNameG || tpreVol == VolNameE )  && 
             ( tpostVol ==  VolNameD ) )
	  {
            if(prePos.z() < postPos.z() )
	    {
	      RichTbAnalysisManager * analysis = RichTbAnalysisManager::getInstance();
	      analysis->bumpNumPhotonsBeforeAerogel();
            }
	  }
        }
#endif

        if( ( tpreVol ==  VolNameD ) && 
            ( tpostVol ==  VolNameG || tpostVol == VolNameE ) )
	{

          // G4double XatAgelExit=postPos.x();
          // G4double YatAgelExit=postPos.y();
          // G4double ZatAgelExit=postPos.z();
          const G4ThreeVector PhotCurMom =  
                aPhotTrack->GetMomentumDirection();
          // G4double CurExitangle= std::acos(PhotCurMom.z());
 
          // Plot the Angle of emission of the photon.
          // When there is no Raylegh scattering this is the
          // Cherenkov angle. For now we only consider the charged
          // track to be of direction 001. Later this may be changed.
          // So the angle considered is just the angle of the photon
          // track.

#ifdef G4ANALYSIS_USE
          RichTbAnalysisManager * analysis =
	    RichTbAnalysisManager::getInstance();

          const G4ThreeVector PhotOrgUnitMom =  
            aPhotTrack->GetVertexMomentumDirection();
          G4double Ckv_angle= std::acos(PhotOrgUnitMom.z());

          analysis->getfhistoCkvProdSmall()->fill(Ckv_angle);
      
          // Now for the photon emission point in aerogel
      
          const G4ThreeVector PhotEmisPt = aPhotTrack->GetVertexPosition();
          analysis->getfhistoEmisZ()->fill( PhotEmisPt.z());         
#endif
        }
      }
     }
    }
   }
  }
}

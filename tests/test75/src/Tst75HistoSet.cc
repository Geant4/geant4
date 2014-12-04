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
#include "globals.hh"
#include "G4PhysicalConstants.hh"

#include "Tst75HistoSet.hh"

#include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"


Tst75HistoSet::Tst75HistoSet( G4String htitle, G4float maxke )
   : TstHistoSet( htitle ), fMaxKE( maxke )
{

   TString title( htitle.c_str() ) ;    
   //
   // plots for secondary proton
   //
   fProton45  = new TH1F("p45deg", title+"+ p (45 deg)",100,0.0,fMaxKE);
   fProton60  = new TH1F("p60deg", title+"+ p (60 deg)",100,0.0,fMaxKE);
   fProton72  = new TH1F("p72deg", title+"+ p (72 deg)",100,0.0,fMaxKE);
   fProton84  = new TH1F("p84deg", title+"+ p (84 deg)",100,0.0,fMaxKE);
   fProton90  = new TH1F("p90deg", title+"+ p (90 deg)",100,0.0,fMaxKE);
   fProton135 = new TH1F("p135deg",title+"+ p (135 deg)",100,0.0,fMaxKE);
   fProKE     = new TH1F("proKE",  title+"proton; Kinetic Energy (MeV)",100,0.0,fMaxKE);
   fProCosTh  = new TH1F("pCosTheta", title+"p;cos #Theta",  100,-1.,1.);

    //
    // plots for secondary pi-
    //    
    fPim28deg = new TH1F("pim28deg", title+"+ pi- (28.4 deg)", 100, 0., fMaxKE );
    fPim28deg->GetXaxis()->SetTitle("momentum of secondary pi- (MeV/c)");
    fPim44deg = new TH1F("pim44deg", title+"+ pi- (44.2 deg)", 100, 0., fMaxKE );
    fPim44deg->GetXaxis()->SetTitle("momentum of secondary pi- (MeV/c)");

    //
    // plots for secondary pi+
    //
    fPip28deg = new TH1F("pip28deg", title+"+ pi+ (28.4 deg)", 100, 0., fMaxKE );
    fPip28deg->GetXaxis()->SetTitle("momentum of secondary pi+ (MeV/c)");
    fPip44deg = new TH1F("pip44deg", title+"+ pi+ (44.2 deg)", 100, 0., fMaxKE );
    fPip44deg->GetXaxis()->SetTitle("momentum of secondary pi+ (MeV/c)");
  
   // fInteraction = new G4VParticleChange();


}

Tst75HistoSet::~Tst75HistoSet()
{

   if (fProton45)  delete fProton45;
   if (fProton60)  delete fProton60;
   if (fProton72)  delete fProton72;   
   if (fProton84)  delete fProton84;
   if (fProton90)  delete fProton90;
   if (fProton135) delete fProton135;
   
   if (fPim28deg)  delete fPim28deg;
   if (fPim44deg)  delete fPim44deg;
   if (fPip28deg)  delete fPip28deg;
   if (fPip44deg)  delete fPip44deg;

}

void Tst75HistoSet::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& ) 
{

      // loop over secondaries and cleanup
      G4int nsecnd = aChange->GetNumberOfSecondaries();

      G4float wei = 1.;
      G4float binFactor = 1.;
      
      for (G4int i=0; i<nsecnd; ++i) 
      {
        G4Track*           sTrack = aChange->GetSecondary(i);
// -->        G4TrackStatus sts =  	sTrack->GetTrackStatus ();
        const G4DynamicParticle* sdp = sTrack->GetDynamicParticle();
        const G4ParticleDefinition* sdpd = sdp->GetDefinition();
        G4String sName = sdpd->GetParticleName();
	
	G4float secKE = sTrack->GetKineticEnergy();
	G4float secP  = sTrack->GetMomentum().mag();
	G4float cosTheta = sTrack->GetMomentumDirection().z();
	
	// secondary proton
	//
	if ( sName == "proton" )
	{
	   fProKE->Fill(secKE,1.0);
	   fProCosTh->Fill(cosTheta,1.0);
	   if (cosTheta > 0.6 && cosTheta <= 0.8) 
	   {
	      binFactor = 0.2*twopi * fProton45->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fProton45->Fill(secKE,wei);
	   } 
	   else if (cosTheta > 0.4 && cosTheta <= 0.6) 
	   {
	      binFactor = 0.2*twopi * fProton60->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fProton60->Fill(secKE,wei);
	   } 
	   else if (cosTheta > 0.2 && cosTheta <= 0.4) 
	   {
	      binFactor = 0.2*twopi * fProton72->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fProton72->Fill(secKE,wei);
	   } 
	   else if (cosTheta > 0.0 && cosTheta <= 0.2) 
	   {
	      binFactor = 0.2*twopi * fProton84->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fProton84->Fill(secKE,wei);
	   } 
	   else if (cosTheta > -0.1 && cosTheta < 0.1) 
	   {
	      binFactor = 0.2*twopi * fProton90->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fProton90->Fill(secKE,wei);
	   } 
	   else if (cosTheta > -0.6 && cosTheta <= 0.8) 
	   {
	      binFactor = 0.2*twopi * fProton135->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fProton135->Fill(secKE,wei);
	   }
	}
	//
	// secondary pi-
	//
	else if ( sName == "pi-" )
	{	   
	   if ( cosTheta > 0.7 && cosTheta <= 0.75 )
	   {
	      // pi- out at 44.2 degree (cos(44.2)=0.7169)
	      //
	      binFactor = 0.05*twopi * fPim44deg->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fPim44deg->Fill( secP, wei );
	   }
	   else if ( cosTheta > 0.85 && cosTheta <= 0.9 )
	   {
	      // pi- out at 28.4 (cos(28.4)=0.8796)
	      //
	      binFactor = 0.05*twopi * fPim28deg->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fPim28deg->Fill( secP, wei );
	   }
	}
	//
	// secondary pi+
	//
	else if ( sName == "pi+" )
	{
	   if ( cosTheta > 0.7 && cosTheta <= 0.75 )
	   {
	      // pi- out at 44.2 degree (cos(44.2)=0.7169)
	      //
	      binFactor = 0.05*twopi * fPip44deg->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fPip44deg->Fill( secP, wei );
	   }
	   else if ( cosTheta > 0.85 && cosTheta <= 0.9 )
	   {
	      // pi- out at 28.4 (cos(28.4)=0.8796)
	      //
	      binFactor = 0.05*twopi * fPip28deg->GetBinWidth(1);
	      wei = 1. / binFactor;
	      fPip28deg->Fill( secP, wei );
	   }
	}

      } // end loop over secondaries

      
   return;
   
}

void Tst75HistoSet::Write( G4int stat, G4double weight )
{

   float scale = weight / ((float)stat);
      
   fProton45->Scale( scale );
   fProton45->Write();
   fProton60->Scale( scale );
   fProton60->Write();
   fProton72->Scale( scale );
   fProton72->Write();
   fProton84->Scale( scale );
   fProton84->Write();
   fProton90->Scale( scale );
   fProton90->Write();
   fProton135->Scale( scale );
   fProton135->Write();
   
   fProKE->Write();
   fProCosTh->Write();

   fPim28deg->Scale( scale );
   fPim28deg->Write();
   fPim44deg->Scale( scale );
   fPim44deg->Write();
   
   fPip28deg->Scale( scale );
   fPip28deg->Write();
   fPip44deg->Scale( scale );
   fPip44deg->Write();

   return;

}

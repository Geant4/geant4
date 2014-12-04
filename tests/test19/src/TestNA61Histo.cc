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

#include "TestNA61Histo.hh"

#include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"

TestNA61Histo::TestNA61Histo( G4String htitle )
   : TstHistoSet(htitle)
{

  G4String title;
  G4String outcome;
  
  outcome = " -> X + proton ";
       
  title = htitle + outcome + " (0<theta<20 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_0_20",title.c_str(),250,0.,25.) ); 
  title = htitle + outcome + " (20<theta<40 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_20_40",title.c_str(),250,0.,25.) ); 
  title = htitle + outcome + " (40<theta<60 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_40_60",title.c_str(),250,0.,25.) ); 
  title = htitle + outcome + " (60<theta<100 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_60_100",title.c_str(),250,0.,25.) ); 
  title = htitle + outcome + " (100<theta<140 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_100_140",title.c_str(),150,0.,15.) ); 
  title = htitle + outcome + " (140<theta<180 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_140_180",title.c_str(),150,0.,15.) ); 
  title = htitle + outcome + " (180<theta<240 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_180_240",title.c_str(),150,0.,15.) ); 
  title = htitle + outcome + " (240<theta<300 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_240_300",title.c_str(),150,0.,15.) ); 
      
  outcome = " -> X + pi- ";

  title = htitle + outcome + " (0<theta<20 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_0_20", title.c_str(), 250, 0., 25. ) );
  title = htitle + outcome + " (20<theta<40 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_20_40", title.c_str(), 250, 0., 25. ) );
  title = htitle + outcome + " (40<theta<60 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_40_60", title.c_str(), 250, 0., 25. ) );
  title = htitle + outcome + " (60<theta<100 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_60_100", title.c_str(), 250, 0., 25. ) );
  title = htitle + outcome + " (100<theta<140 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_100_140", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (140<theta<180 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_140_180", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (180<theta<240 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_180_240", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (240<theta<300 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_240_300", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (300<theta<360 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_300_360", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (360<theta<420 (mrad)) ";
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_360_420", title.c_str(), 150, 0., 15. ) );

  outcome = " -> X + pi+ ";

  title = htitle + outcome + " (0<theta<20 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_0_20", title.c_str(), 250, 0., 25. ) );
  title = htitle + outcome + " (20<theta<40 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_20_40", title.c_str(), 250, 0., 25. ) );
  title = htitle + outcome + " (40<theta<60 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_40_60", title.c_str(), 250, 0., 25. ) );
  title = htitle + outcome + " (60<theta<100 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_60_100", title.c_str(), 250, 0., 25. ) );
  title = htitle + outcome + " (100<theta<140 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_100_140", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (140<theta<180 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_140_180", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (180<theta<240 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_180_240", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (240<theta<300 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_240_300", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (300<theta<360 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_300_360", title.c_str(), 150, 0., 15. ) );
  title = htitle + outcome + " (360<theta<420 (mrad)) ";
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_360_420", title.c_str(), 150, 0., 15. ) );
  
  outcome = " -> X + K+ ";
  
  title = htitle + outcome + " (20<theta<140 (mrad)) ";
  fHistoSecKPlus.push_back( new TH1F( "kplusMult_20_140",     title.c_str(), 100, 0., 10. ) );
  title = htitle + outcome + " (140<theta<240 (mrad)) ";
  fHistoSecKPlus.push_back( new TH1F( "kplusMult_140_240",    title.c_str(), 100, 0., 10. ) );
  title = htitle + ", K+/pi+ (20<theta<140 (mrad)) ";
  fHistoSecKPlus.push_back( new TH1F( "kplus2piplus_20_140",  title.c_str(), 100, 0., 10. ) );
  title = htitle + ", K+/pi+ (140<theta<240 (mrad)) ";
  fHistoSecKPlus.push_back( new TH1F( "kplus2piplus_140_240",  title.c_str(), 100, 0., 10. ) );
  
  fHistoSecPiPlus2.push_back( new TH1F( "piplus_20_140",  " ", 100, 0., 10. ) );
  fHistoSecPiPlus2.push_back( new TH1F( "piplus_140_240", " ", 100, 0., 10. ) );

}

TestNA61Histo::~TestNA61Histo()
{

   for (size_t i=0; i<fHistoSecProton.size(); i++ )
   {
      delete fHistoSecProton[i];
   }
   for (size_t i=0; i<fHistoSecPiMinus.size(); i++ )
   {
      delete fHistoSecPiMinus[i];
   }
   for (size_t i=0; i<fHistoSecPiPlus.size(); i++ )
   {
      delete fHistoSecPiPlus[i];
   }
   for ( size_t i=0; i<fHistoSecKPlus.size(); i++ )
   {
      delete fHistoSecKPlus[i];
   }
   for ( size_t i=0; i<fHistoSecPiPlus2.size(); i++ )
   {
      delete fHistoSecPiPlus2[i];
   }

}

void TestNA61Histo::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& ) 
{

   G4int NSec = aChange->GetNumberOfSecondaries();
     
   const G4DynamicParticle* sec = 0;
      
   for (G4int i=0; i<NSec; i++) 
   {

      sec = aChange->GetSecondary(i)->GetDynamicParticle();
			
      const G4String& pname = sec->GetDefinition()->GetParticleName();
	
      // G4ThreeVector mom = sec->GetMomentum();
	
      double theta = sec->GetMomentum().theta();
      theta /= mrad;
	
      double pmom = sec->GetTotalMomentum() / GeV ;
	
      int id = -1;
	
      int NHProSize = fHistoSecProton.size();
      int NHPimSize = fHistoSecPiMinus.size();
      int NHPipSize = fHistoSecPiPlus.size();
      
      if ( pname == "proton" || pname == "pi+" || pname == "pi-" )
      {
	if ( theta >= 0. && theta < 20. )
	{
	   id = 0;
	}
	else if ( theta >= 20. && theta < 40. )
	{
	   id = 1;
	}
	else if ( theta >= 40. && theta < 60. )
	{
	   id = 2;
	}
	else if ( theta >= 60. && theta < 100. )
	{
	   id = 3;
	}
	else if ( theta >= 100. && theta < 140. )
	{
	   id = 4;
	}
	else if ( theta >= 140. && theta < 180. )
	{
	   id = 5;
	}
	else if ( theta > 180. && theta < 240. )
	{
	   id = 6;
	}
	else if ( theta >= 240. && theta < 300. )
	{
	   id = 7;
	} 
	else if ( theta >= 300. && theta < 360. )
	{
	   id = 8;
	}
	else if ( theta >= 360. && theta < 420. )
	{
	   id = 9;
	}
	
	if ( pname == "proton" )
	{	   
	   if ( id >=0 && id < NHProSize ) fHistoSecProton[id]->Fill( pmom );
	}
	if ( pname == "pi-" )
	{
	   if ( id >=0 && id < NHPimSize ) fHistoSecPiMinus[id]->Fill( pmom );
	}
	if ( pname == "pi+" )
	{
	   if ( id >= 0 && id < NHPipSize ) fHistoSecPiPlus[id]->Fill( pmom );
	}
      }
      
      if ( pname == "proton" || pname == "pi-" ) continue;

      id = -1;

      if ( pname == "kaon+" || pname == "pi+" )
      {
         if ( theta >=20. && theta < 140. )
	 {
	    id = 0;
	 }
	 else if ( theta >= 140. && theta < 240. )
	 {
	    id = 1;
	 }
	 if ( pname == "kaon+" )
	 {
	    if ( id >= 0 )
	    {
	       fHistoSecKPlus[id]->Fill( pmom ) ;
	       fHistoSecKPlus[id+2]->Fill( pmom ) ;
	    }
	 }
	 else if ( pname == "pi+" )
	 {
	    if ( id >= 0 ) fHistoSecPiPlus2[id]->Fill( pmom ) ;
	 }	 
      }
         
   }
      
   return;
   
}

void TestNA61Histo::Write( G4int stat, G4double )
{

   for ( size_t i=0; i<fHistoSecProton.size(); ++i )
   {
      double xbin = fHistoSecProton[i]->GetBinWidth(1);
      double scale = 1. / ((double)stat * xbin) ;
      fHistoSecProton[i]->Scale(scale);
      fHistoSecProton[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiMinus.size(); ++i )
   {
      double xbin = fHistoSecPiMinus[i]->GetBinWidth(1);
      double scale = 1. / ((double)stat * xbin) ;
      fHistoSecPiMinus[i]->Scale(scale);
      fHistoSecPiMinus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiPlus.size(); ++i )
   {
      double xbin = fHistoSecPiPlus[i]->GetBinWidth(1);
      double scale = 1. / ((double)stat * xbin) ;
      fHistoSecPiPlus[i]->Scale(scale);
      fHistoSecPiPlus[i]->Write();
   }
   
   for ( size_t i=0; i<fHistoSecKPlus.size()-2; i++ )
   {
      double xbin = fHistoSecKPlus[i]->GetBinWidth(1);
      double scale = 1. / ((double)stat * xbin) ;
      fHistoSecKPlus[i]->Scale(scale);
      fHistoSecKPlus[i]->Write();
      fHistoSecKPlus[i+2]->Divide( fHistoSecPiPlus2[i] );
      fHistoSecKPlus[i+2]->Write();
   }
   
   return;

}

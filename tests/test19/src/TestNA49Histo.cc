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

#include "TestNA49Histo.hh"

#include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"


TestNA49Histo::TestNA49Histo( G4String htitle )
   : TstHistoSet( htitle )
{

  std::string title;
  std::string outcome;
    
  fHistoNSec = new TH1D( "NSec", htitle.c_str(), 50, 0., 50. );
  
  outcome = " -> X + proton ";       
  title = htitle + outcome;

  // in principle, the exp.data are binned by 0.05 from -0.8 to 0.95, 
  // except in xF=(-0.1,0.1) where the binning is finer (0.025) 
  // so for now we do 0.025 everywhere, for simplicity
  //
  fHistoSecProton.push_back( new TH1D( "proton_dXSecdxF", title.c_str(), 81, -1.05, 1.05 ) );	
  fHistoSecProton.push_back( new TH1D( "proton_dNdxF",    title.c_str(), 81, -1.05, 1.05 ) );	
  fHistoSecProton.push_back( new TProfile( "proton_pT",   title.c_str(), 81, -1.05, 1.05, 0., 10. ) );	
  fHistoSecProton.push_back( new TProfile( "proton_pT2",  title.c_str(), 81, -1.05, 1.05, 0., 10. ) );	
       
  outcome = " -> X + antiproton ";
  title = htitle + outcome;
  
//  double pbarbins[] = { -0.3, -0.2, -0.15, -0.1, -0.075, -0.05, -0.025, 
//                         0.0,
//			 0.025, 0.05, 0.1, 0.15, 0.2, 0.3 }
  // we specify bin left edge, to make the bin center match the exp.data
  //
  double pbarbins[] = { -0.55, -0.45, -0.35, -0.25, -0.175, -0.125, -0.0875, -0.0625, -0.0375, 
                         -0.0125,
			 0.0125, 0.0375, 0.075, 0.125, 0.175, 0.25, 0.35, 0.45, 0.55 };
  int npbarbins = sizeof(pbarbins) / sizeof(double) - 1;
  
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dXSecdxF", title.c_str(), 80, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dNdxF",    title.c_str(), npbarbins, pbarbins ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT",   title.c_str(), npbarbins, pbarbins, 0., 10. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT2",  title.c_str(), npbarbins, pbarbins, 0., 10. ) );	
  //fHistoSecAntiProton.push_back( new TH1D( "antiproton_dNdxF",    title.c_str(), 80, -1., 1. ) );	
  //fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT",   title.c_str(), 80, -1., 1., 0., 10. ) );	
  //fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT2",  title.c_str(), 80, -1., 1., 0., 10. ) );	
  
  outcome = " -> X + pi- ";  

  title = htitle + outcome;
//  double pionbins[] = { -0.500, -0.400, -0.300, -0.250, -0.200, -0.150, -0.125, -0.100, -0.075, 
//                     -0.050, -0.040, -0.030, -0.020, -0.010,  
//		      0.0,
//		      0.010,  0.020,  0.030,  0.040,  0.050,
//		      0.075,  0.100,  0.125,  0.150,  0.200,  0.250,  0.300,  0.400,  0.500 }; 
		       
  // the idea is to make the bin center correspond to the number in the NA49 data file(s)
  //
  double pibins[] = { -0.550, -0.450, -0.350, -0.275, -0.225, -0.175, -0.1375, -0.1125, -0.0875, 
                     -0.0625, -0.045, -0.035, -0.025, -0.015,  
		      -0.005,
		      0.005,  0.015,  0.025,  0.035,  0.045,
		      0.0625,  0.0875,  0.1125,  0.1375,  0.175, 0.225,  0.275,  0.350,  0.450, 0.55 };  
  int npibins = sizeof(pibins) / sizeof(double) - 1;

  fHistoSecPiMinus.push_back( new TH1D( "pimimus_dXSecdxF", title.c_str(), 200, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TH1D( "piminus_dNdxF",    title.c_str(), npibins, pibins ) );
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT",   title.c_str(), npibins, pibins, 0., 10. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT2",  title.c_str(), npibins, pibins, 0., 10. ) );	
  // fHistoSecPiMinus.push_back( new TH1D( "piminus_dNdxF",    title.c_str(), 200, -1., 1. ) );	  
  // fHistoSecPiMinus.push_back( new TProfile( "piminus_pT",   title.c_str(), 200, -1., 1., 0., 10. ) );	
  // fHistoSecPiMinus.push_back( new TProfile( "piminus_pT2",  title.c_str(), 200, -1., 1., 0., 10. ) );	


  outcome = " -> X + pi+ ";
  title = htitle + outcome;

  fHistoSecPiPlus.push_back( new TH1D( "piplus_dXSecdxF", title.c_str(), 200, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TH1D( "piplus_dNdxF",    title.c_str(), npibins, pibins ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT",   title.c_str(), npibins, pibins, 0., 10. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT2",  title.c_str(), npibins, pibins, 0., 10. ) );	
  // fHistoSecPiPlus.push_back( new TH1D( "piplus_dNdxF",    title.c_str(), 200, -1., 1. ) );	
  // fHistoSecPiPlus.push_back( new TProfile( "piplus_pT",   title.c_str(), 200, -1., 1., 0., 10. ) );	
  // fHistoSecPiPlus.push_back( new TProfile( "piplus_pT2",  title.c_str(), 200, -1., 1., 0., 10. ) );	

  outcome = " -> X + neutron ";
  title = htitle + outcome;

  // this binning is reasonable because the data go by 0.1 in (0.1-0.6), then just 2 bins of 0.15
  //
  fHistoSecNeutron.push_back( new TH1D( "neutron_dNdxF",    title.c_str(), 10, 0.05, 1.05 ) );	
  
  fInteraction = new G4VParticleChange();


}

TestNA49Histo::~TestNA49Histo()
{

   if (fHistoNSec) delete fHistoNSec;
   
   for (size_t i=0; i<fHistoSecProton.size(); i++ )
   {
      delete fHistoSecProton[i];
   }
   for (size_t i=0; i<fHistoSecAntiProton.size(); i++ )
   {
      delete fHistoSecAntiProton[i];
   }
   for (size_t i=0; i<fHistoSecPiMinus.size(); i++ )
   {
      delete fHistoSecPiMinus[i];
   }
   for (size_t i=0; i<fHistoSecPiPlus.size(); i++ )
   {
      delete fHistoSecPiPlus[i];
   }
   for ( size_t i=0; i<fHistoSecNeutron.size(); i++ )
   {
      delete fHistoSecNeutron[i];
   }
   
   if ( fInteraction )
   {
      int nsc = fInteraction->GetNumberOfSecondaries();
      if ( nsc > 0  )
      {
         for ( int i=0; i<nsc; i++) 
         {   
            delete fInteraction->GetSecondary(i);
         } 
      }
      fInteraction->Clear();
      delete fInteraction;    
   }
      
}

void TestNA49Histo::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& labp ) 
{


   G4int NSec = fInteraction->GetNumberOfSecondaries();
   if ( NSec > 0 )
   {
      for ( int i=0; i<NSec; i++ )
      {
         delete fInteraction->GetSecondary(i);
      }
      fInteraction->Clear();
   }
   
   NSec = aChange->GetNumberOfSecondaries();
   
   for ( int i=0; i<NSec; i++ )
   {
      G4Track* trk = new G4Track( *(aChange->GetSecondary(i)) );
      fInteraction->AddSecondary( trk );
   }

   if ( fDoResDecay ) AccountForResDecay( fInteraction );
   
   double SQRT_S = 17.2;
   
   G4ThreeVector boostLabp = labp.boostVector();
   //G4cout << " boostLabp = " << boostLabp << G4endl;
   
   
   NSec = fInteraction->GetNumberOfSecondaries();
   
   fHistoNSec->Fill( (double)NSec );
     
   const G4DynamicParticle* sec = 0;
       
   for (G4int i=0; i<NSec; i++) 
   {
        sec = fInteraction->GetSecondary(i)->GetDynamicParticle();			
	const G4String& pname = sec->GetDefinition()->GetParticleName();
	
	G4ThreeVector mom = sec->GetMomentum() / GeV ;
	double mass = sec->GetDefinition()->GetPDGMass() / GeV;
	double ekin = sec->GetKineticEnergy() / GeV ;
	
	G4LorentzVector boostSec( mom, ekin+mass );
	//G4cout << " boostSec = " << boostSec << G4endl;
	boostSec.boost(-boostLabp);
	//G4cout << " back-boosted = " << boostSec << G4endl;
	
	double xF  = 2 * (boostSec.z()) / SQRT_S;
	
	double pT  = mom.perp() ;
	double pT2 = pT * pT ;
        
	if ( pname == "neutron" )
	{
	   fHistoSecNeutron[0]->Fill( xF );
	}
	else if ( pname == "pi-" )
	{
	   fHistoSecPiMinus[1]->Fill( xF );
	   fHistoSecPiMinus[2]->Fill( xF, pT );
	   fHistoSecPiMinus[3]->Fill( xF, pT2 );
	}
	else if ( pname == "pi+" )
	{
	   fHistoSecPiPlus[1]->Fill( xF );
	   fHistoSecPiPlus[2]->Fill( xF, pT );
	   fHistoSecPiPlus[3]->Fill( xF, pT2 );
	}
	else if ( pname == "proton" )
	{
	   fHistoSecProton[1]->Fill( xF );
	   fHistoSecProton[2]->Fill( xF, pT );
	   fHistoSecProton[3]->Fill( xF, pT2 );
	}
	else if ( pname == "anti_proton" )
	{
	   fHistoSecAntiProton[1]->Fill( xF );
	   fHistoSecAntiProton[2]->Fill( xF, pT );
	   fHistoSecAntiProton[3]->Fill( xF, pT2 );
	}	
   }
      
   return;
   
}

void TestNA49Histo::Write( G4int stat, G4double )
{

   double xbin = 1.;
   double scale = 1.;
   
   xbin = fHistoNSec->GetBinWidth(1);
   scale = 1. / ((double)stat * xbin);
   fHistoNSec->Scale( scale ) ;
   fHistoNSec->Write();
   
   for ( size_t i=0; i<fHistoSecProton.size(); ++i )   
   {
      xbin = fHistoSecProton[i]->GetBinWidth(1);
      scale = 1. / ((double)stat * xbin) ;
      if ( i <2 ) fHistoSecProton[i]->Scale(scale);
      fHistoSecProton[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecAntiProton.size(); ++i )   
   {
      if ( i < 2 )
      {
         int nbins = fHistoSecAntiProton[i]->GetXaxis()->GetNbins();
         for ( int ib = 1; ib <= nbins; ib++ )
         { 
            xbin = fHistoSecAntiProton[i]->GetBinWidth(ib);
	    scale = 1. / ((double)stat * xbin) ;
	    double yvalue = fHistoSecAntiProton[i]->GetBinContent(ib);
	    yvalue *= scale;
	    fHistoSecAntiProton[i]->SetBinContent( ib, yvalue );
         }
      }
      //scale = 1. / ((double)stat * xbin) ;
      //if ( i <2 ) fHistoSecAntiProton[i]->Scale(scale);
      fHistoSecAntiProton[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiMinus.size(); ++i )
   {
      if ( i < 2 )
      {
         int nbins = fHistoSecPiMinus[i]->GetXaxis()->GetNbins();
         for ( int ib = 1; ib <= nbins; ib++ )
         { 
            xbin = fHistoSecPiMinus[i]->GetBinWidth(ib);
	    scale = 1. / ((double)stat * xbin) ;
	    double yvalue = fHistoSecPiMinus[i]->GetBinContent(ib);
	    yvalue *= scale;
	    // double yerror = fHistoSecPiMinus[i]->GetBinError(ib);
	    // yerror *= scale;
	    fHistoSecPiMinus[i]->SetBinContent( ib, yvalue );
	    // fHistoSecPiMinus[i]->SetBinError( ib, yerror );
         }
      }
      //scale = 1. / ((double)stat * xbin) ;
      //if ( i < 2 ) fHistoSecPiMinus[i]->Scale(scale);
      fHistoSecPiMinus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiPlus.size(); ++i )
   {
      if ( i < 2 )
      {
         int nbins = fHistoSecPiPlus[i]->GetXaxis()->GetNbins();
         for ( int ib = 1; ib <= nbins; ib++ )
         { 
            xbin = fHistoSecPiPlus[i]->GetBinWidth(ib);
	    scale = 1. / ((double)stat * xbin) ;
	    double yvalue = fHistoSecPiPlus[i]->GetBinContent(ib);
	    yvalue *= scale;
	    // double yerror = fHistoSecPiPlus[i]->GetBinError(ib);
	    // yerror *= scale;
	    fHistoSecPiPlus[i]->SetBinContent( ib, yvalue );
	    // fHistoSecPiPlus[i]->SetBinError( ib, yerror );
         }
      }
      //scale = 1. / ((double)stat * xbin) ;
      //if ( i < 2 ) fHistoSecPiPlus[i]->Scale(scale);
      fHistoSecPiPlus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecNeutron.size(); ++i )
   {
     xbin = fHistoSecNeutron[i]->GetBinWidth(1);
     // double scale = 1. / ((double)stat * xbin);
     scale = 1. / ((double)stat * xbin);
     if ( i < 2 ) fHistoSecNeutron[i]->Scale(scale);
     fHistoSecNeutron[i]->Write();
   }

   return;

}

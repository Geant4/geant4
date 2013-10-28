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


#include "Tst23NA49Histo.hh"
#include "Tst23ParticleChange.hh"

// #include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"

Tst23NA49Histo::Tst23NA49Histo( G4String htitle )
   : TstHistoSet( htitle )
{

  G4String title;
  G4String outcome;
  
  fHistoNSec = new TH1D( "NSec", htitle.c_str(), 50, 0., 50. );
  
  outcome = " -> X + proton ";       
  title = htitle + outcome;

  fHistoSecProton.push_back( new TH1D( "proton_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecProton.push_back( new TH1D( "proton_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecProton.push_back( new TProfile( "proton_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecProton.push_back( new TProfile( "proton_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  fHistoSecProtonTot.push_back( new TH1D( "proton_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecProtonTot.push_back( new TH1D( "proton_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecProtonTot.push_back( new TProfile( "proton_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecProtonTot.push_back( new TProfile( "proton_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	
       
  outcome = " -> X + antiproton ";
  title = htitle + outcome;
  
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  fHistoSecAntiProtonTot.push_back( new TH1D( "antiproton_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProtonTot.push_back( new TH1D( "antiproton_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProtonTot.push_back( new TProfile( "antiproton_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecAntiProtonTot.push_back( new TProfile( "antiproton_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	
  
  outcome = " -> X + pi- ";  
  title = htitle + outcome;

  fHistoSecPiMinus.push_back( new TH1D( "pimimus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TH1D( "piminus_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  fHistoSecPiMinusTot.push_back( new TH1D( "pimimus_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinusTot.push_back( new TH1D( "piminus_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinusTot.push_back( new TProfile( "piminus_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiMinusTot.push_back( new TProfile( "piminus_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  outcome = " -> X + pi+ ";
  title = htitle + outcome;

  fHistoSecPiPlus.push_back( new TH1D( "piplus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TH1D( "piplus_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  fHistoSecPiPlusTot.push_back( new TH1D( "piplus_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlusTot.push_back( new TH1D( "piplus_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlusTot.push_back( new TProfile( "piplus_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiPlusTot.push_back( new TProfile( "piplus_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  outcome = " -> X + neutron ";
  title = htitle + outcome;

  fHistoSecNeutron.push_back( new TH1D( "neutron_dNdxF", title.c_str(), 10, 0.05, 1.05 ) );	

  fHistoSecNeutronTot.push_back( new TH1D( "neutron_dNdxF_total", title.c_str(), 10, 0.05, 1.05 ) );	

}

Tst23NA49Histo::~Tst23NA49Histo()
{

   if ( fHistoNSec ) delete fHistoNSec;
   
   for (size_t i=0; i<fHistoSecProton.size(); i++ )
   {
      delete fHistoSecProton[i];
   }
   for (size_t i=0; i<fHistoSecProtonTot.size(); i++ )
   {
      delete fHistoSecProtonTot[i];
   }
   
   for (size_t i=0; i<fHistoSecAntiProton.size(); i++ )
   {
      delete fHistoSecAntiProton[i];
   }
   for (size_t i=0; i<fHistoSecAntiProtonTot.size(); i++ )
   {
      delete fHistoSecAntiProtonTot[i];
   }
   
   for (size_t i=0; i<fHistoSecPiMinus.size(); i++ )
   {
      delete fHistoSecPiMinus[i];
   }
   for (size_t i=0; i<fHistoSecPiMinusTot.size(); i++ )
   {
      delete fHistoSecPiMinusTot[i];
   }
   
   for (size_t i=0; i<fHistoSecPiPlus.size(); i++ )
   {
      delete fHistoSecPiPlus[i];
   }
   for (size_t i=0; i<fHistoSecPiPlusTot.size(); i++ )
   {
      delete fHistoSecPiPlusTot[i];
   }
   
   for ( size_t i=0; i<fHistoSecNeutron.size(); i++ )
   {
      delete fHistoSecNeutron[i];
   }
   for ( size_t i=0; i<fHistoSecNeutronTot.size(); i++ )
   {
      delete fHistoSecNeutronTot[i];
   }

}

/*
void Tst23NA49Histo::InitHisto()
{
     
  std::string title;
  std::string outcome;
  
  fHistoNSec = new TH1D( "NSec", fHistoTitle.c_str(), 50, 0., 50. );
  
  outcome = " -> X + proton ";       
  title = fHistoTitle + outcome;

  fHistoSecProton.push_back( new TH1D( "proton_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecProton.push_back( new TH1D( "proton_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecProton.push_back( new TProfile( "proton_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecProton.push_back( new TProfile( "proton_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  fHistoSecProtonTot.push_back( new TH1D( "proton_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecProtonTot.push_back( new TH1D( "proton_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecProtonTot.push_back( new TProfile( "proton_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecProtonTot.push_back( new TProfile( "proton_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	
       
  outcome = " -> X + antiproton ";
  title = fHistoTitle + outcome;
  
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  fHistoSecAntiProtonTot.push_back( new TH1D( "antiproton_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProtonTot.push_back( new TH1D( "antiproton_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProtonTot.push_back( new TProfile( "antiproton_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecAntiProtonTot.push_back( new TProfile( "antiproton_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	
  
  outcome = " -> X + pi- ";  
  title = fHistoTitle + outcome;

  fHistoSecPiMinus.push_back( new TH1D( "pimimus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TH1D( "piminus_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  fHistoSecPiMinusTot.push_back( new TH1D( "pimimus_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinusTot.push_back( new TH1D( "piminus_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinusTot.push_back( new TProfile( "piminus_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiMinusTot.push_back( new TProfile( "piminus_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  outcome = " -> X + pi+ ";
  title = fHistoTitle + outcome;

  fHistoSecPiPlus.push_back( new TH1D( "piplus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TH1D( "piplus_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  fHistoSecPiPlusTot.push_back( new TH1D( "piplus_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlusTot.push_back( new TH1D( "piplus_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlusTot.push_back( new TProfile( "piplus_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiPlusTot.push_back( new TProfile( "piplus_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  outcome = " -> X + neutron ";
  title = fHistoTitle + outcome;

  fHistoSecNeutron.push_back( new TH1D( "neutron_dNdxF",    title.c_str(), 10, 0.05, 1.05 ) );	

  fHistoSecNeutronTot.push_back( new TH1D( "neutron_dNdxF_total",    title.c_str(), 10, 0.05, 1.05 ) );	

  return;

}
*/

void Tst23NA49Histo::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& labp ) 
{

   
   bool isFirst = (dynamic_cast<Tst23ParticleChange*>(aChange))->IsFisrtInteraction();
   
   // if ( !isFirst ) return;
   
   double SQRT_S = 17.2;
   
   G4ThreeVector boostLabp = labp.boostVector();
   
   G4int NSec = aChange->GetNumberOfSecondaries();

   if ( isFirst && NSec > 0 ) fHistoNSec->Fill( (double)NSec );
     
   const G4DynamicParticle* sec = 0;
       
   for (G4int i=0; i<NSec; i++) 
   {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();			
	const G4String& pname = sec->GetDefinition()->GetParticleName();
	
	G4ThreeVector mom = sec->GetMomentum() / GeV ;
	double mass = sec->GetDefinition()->GetPDGMass() / GeV;
	double ekin = sec->GetKineticEnergy() / GeV ;
	
	G4LorentzVector boostSec( mom, ekin+mass );
	boostSec.boost(-boostLabp);
	
	double xF  = 2 * (boostSec.z()) / SQRT_S;
	
	double pT  = mom.perp() ;
	double pT2 = pT * pT ;
        
	if ( pname == "neutron" )
	{
	   if ( isFirst ) fHistoSecNeutron[0]->Fill( xF );
	   fHistoSecNeutronTot[0]->Fill( xF );
	}
	else if ( pname == "pi-" )
	{
	   if ( isFirst )
	   {
	      fHistoSecPiMinus[1]->Fill( xF );
	      fHistoSecPiMinus[2]->Fill( xF, pT );
	      fHistoSecPiMinus[3]->Fill( xF, pT2 );
	   }
	   fHistoSecPiMinusTot[1]->Fill( xF );
	   fHistoSecPiMinusTot[2]->Fill( xF, pT );
	   fHistoSecPiMinusTot[3]->Fill( xF, pT2 );
	}
	else if ( pname == "pi+" )
	{
	   if ( isFirst )
	   {
	      fHistoSecPiPlus[1]->Fill( xF );
	      fHistoSecPiPlus[2]->Fill( xF, pT );
	      fHistoSecPiPlus[3]->Fill( xF, pT2 );
	   }
	   fHistoSecPiPlusTot[1]->Fill( xF );
	   fHistoSecPiPlusTot[2]->Fill( xF, pT );
	   fHistoSecPiPlusTot[3]->Fill( xF, pT2 );
	}
	else if ( pname == "proton" )
	{
	   if ( isFirst )
	   {
	      fHistoSecProton[1]->Fill( xF );
	      fHistoSecProton[2]->Fill( xF, pT );
	      fHistoSecProton[3]->Fill( xF, pT2 );
	   }
	   fHistoSecProtonTot[1]->Fill( xF );
	   fHistoSecProtonTot[2]->Fill( xF, pT );
	   fHistoSecProtonTot[3]->Fill( xF, pT2 );
	}
	else if ( pname == "anti_proton" )
	{
	   if ( isFirst )
	   {
	      fHistoSecAntiProton[1]->Fill( xF );
	      fHistoSecAntiProton[2]->Fill( xF, pT );
	      fHistoSecAntiProton[3]->Fill( xF, pT2 );
	   }
	   fHistoSecAntiProtonTot[1]->Fill( xF );
	   fHistoSecAntiProtonTot[2]->Fill( xF, pT );
	   fHistoSecAntiProtonTot[3]->Fill( xF, pT2 );
	}	
   }
      
   return;
   
}

void Tst23NA49Histo::Write( G4int, G4double )
{

   // fHistoFile->cd();

   double xbin = 1.;
   double norm = 1.;
   double scale = 1.;
   
   norm = fHistoNSec->Integral();
   
   xbin = (double)fHistoNSec->GetBinWidth(1);
   scale = 1. / (xbin*norm);
   fHistoNSec->Scale( scale ) ;
   fHistoNSec->Write();
   
   // secondary proton
   //
   for ( size_t i=0; i<fHistoSecProton.size(); ++i )   
   {
      xbin = fHistoSecProton[i]->GetBinWidth(1);
      // scale = 1. / ((double)stat * xbin) ;
      scale = 1. / (xbin*norm);
      if ( i <2 ) fHistoSecProton[i]->Scale(scale);
      fHistoSecProton[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecProtonTot.size(); ++i )   
   {
      xbin = fHistoSecProtonTot[i]->GetBinWidth(1);
      // scale = 1. / ((double)stat * xbin) ;
      scale = 1. / (xbin*norm);
      if ( i <2 ) fHistoSecProtonTot[i]->Scale(scale);
      fHistoSecProtonTot[i]->Write();
   }

   // secondary pbar
   //
   for ( size_t i=0; i<fHistoSecAntiProton.size(); ++i )   
   {
      xbin = fHistoSecAntiProton[i]->GetBinWidth(1);
      // scale = 1. / ((double)stat * xbin) ;
      scale = 1. / (xbin*norm);
      if ( i <2 ) fHistoSecAntiProton[i]->Scale(scale);
      fHistoSecAntiProton[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecAntiProtonTot.size(); ++i )   
   {
      xbin = fHistoSecAntiProtonTot[i]->GetBinWidth(1);
      // scale = 1. / ((double)stat * xbin) ;
      scale = 1. / (xbin*norm);
      if ( i <2 ) fHistoSecAntiProtonTot[i]->Scale(scale);
      fHistoSecAntiProtonTot[i]->Write();
   }

   // secondary pi-
   //
   for ( size_t i=0; i<fHistoSecPiMinus.size(); ++i )
   {
      xbin = fHistoSecPiMinus[i]->GetBinWidth(1);
      // scale = 1. / ((double)stat * xbin) ;
      scale = 1. / (xbin*norm);
      if ( i < 2 ) fHistoSecPiMinus[i]->Scale(scale);
      fHistoSecPiMinus[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecPiMinusTot.size(); ++i )
   {
      xbin = fHistoSecPiMinusTot[i]->GetBinWidth(1);
      // scale = 1. / ((double)stat * xbin) ;
      scale = 1. / (xbin*norm);
      if ( i < 2 ) fHistoSecPiMinusTot[i]->Scale(scale);
      fHistoSecPiMinusTot[i]->Write();
   }

   // secondary pi+
   //
   for ( size_t i=0; i<fHistoSecPiPlus.size(); ++i )
   {
      xbin = fHistoSecPiPlus[i]->GetBinWidth(1);
      // scale = 1. / ((double)stat * xbin) ;
      scale = 1. / (xbin*norm);
      if ( i < 2 ) fHistoSecPiPlus[i]->Scale(scale);
      fHistoSecPiPlus[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecPiPlusTot.size(); ++i )
   {
      xbin = fHistoSecPiPlusTot[i]->GetBinWidth(1);
      // scale = 1. / ((double)stat * xbin) ;
      scale = 1. / (xbin*norm);
      if ( i < 2 ) fHistoSecPiPlusTot[i]->Scale(scale);
      fHistoSecPiPlusTot[i]->Write();
   }

   // secondary neutron
   //
   for ( size_t i=0; i<fHistoSecNeutron.size(); ++i )
   {
      xbin = fHistoSecNeutron[i]->GetBinWidth(1);
      // double scale = 1. / ((double)stat * xbin);
      // scale = 1. / ((double)stat * xbin);
      scale = 1. / (xbin*norm);
      if ( i < 2 ) fHistoSecNeutron[i]->Scale(scale);
      fHistoSecNeutron[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecNeutronTot.size(); ++i )
   {
      xbin = fHistoSecNeutronTot[i]->GetBinWidth(1);
      // double scale = 1. / ((double)stat * xbin);
      // scale = 1. / ((double)stat * xbin);
      scale = 1. / (xbin*norm);
      if ( i < 2 ) fHistoSecNeutronTot[i]->Scale(scale);
      fHistoSecNeutronTot[i]->Write();
   }

   return;

}

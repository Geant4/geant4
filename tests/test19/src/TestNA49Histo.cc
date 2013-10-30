
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

  fHistoSecProton.push_back( new TH1D( "proton_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecProton.push_back( new TH1D( "proton_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecProton.push_back( new TProfile( "proton_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecProton.push_back( new TProfile( "proton_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	
       
  outcome = " -> X + antiproton ";
  title = htitle + outcome;
  
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	
  
  outcome = " -> X + pi- ";  
  title = htitle + outcome;

  fHistoSecPiMinus.push_back( new TH1D( "pimimus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TH1D( "piminus_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT2",  title.c_str(), 100, -1., 1., 0., 10. ) );	


  outcome = " -> X + pi+ ";
  title = htitle + outcome;

  fHistoSecPiPlus.push_back( new TH1D( "piplus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TH1D( "piplus_dNdxF",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT2",      title.c_str(), 100, -1., 1., 0., 10. ) );	

  outcome = " -> X + neutron ";
  title = htitle + outcome;

  fHistoSecNeutron.push_back( new TH1D( "neutron_dNdxF",    title.c_str(), 10, 0.05, 1.05 ) );	


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

}

void TestNA49Histo::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& labp ) 
{

   double SQRT_S = 17.2;
   
   G4ThreeVector boostLabp = labp.boostVector();
   
   G4int NSec = aChange->GetNumberOfSecondaries();
   
   fHistoNSec->Fill( (double)NSec );
     
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
      xbin = fHistoSecAntiProton[i]->GetBinWidth(1);
      scale = 1. / ((double)stat * xbin) ;
      if ( i <2 ) fHistoSecAntiProton[i]->Scale(scale);
      fHistoSecAntiProton[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiMinus.size(); ++i )
   {
      xbin = fHistoSecPiMinus[i]->GetBinWidth(1);
      scale = 1. / ((double)stat * xbin) ;
      if ( i < 2 ) fHistoSecPiMinus[i]->Scale(scale);
      fHistoSecPiMinus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiPlus.size(); ++i )
   {
      xbin = fHistoSecPiPlus[i]->GetBinWidth(1);
      scale = 1. / ((double)stat * xbin) ;
      if ( i < 2 ) fHistoSecPiPlus[i]->Scale(scale);
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

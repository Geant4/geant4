
#include "TestStoppingHisto.hh"

#include "G4VParticleChange.hh"

TestStoppingHisto::~TestStoppingHisto()
{

   for (size_t i=0; i<fHisto.size(); i++ )
   {
      delete fHisto[i];
   }
   fHisto.clear();

}

void TestStoppingHisto::Init()
{
   
   fHistoTitle = fBeam + " on " + fTarget;
   
   InitHistoGeneral();
   
   if ( fBeam == "pi-" )
   {
      InitPionMinus();
   }
   
   return;

}

void TestStoppingHisto::InitPionMinus()
{

   std::string title;
   title = fHistoTitle + "; Neutron Yield vs Kinetic Energy (MeV)";
   
   fHisto.push_back(new TH1F( "NvsT", title.c_str(), 600, 1., 151.));

   return;

}

void TestStoppingHisto::InitHistoGeneral()
{
   
   fHisto.push_back( new TH1F( "NSecondaries", fHistoTitle.c_str(), 25, 0., 25. ) );
   fHisto.push_back( new TH1F( "NChargedSecondaries", fHistoTitle.c_str(), 25, 0., 25. ) );
   fHisto.push_back( new TH1F( "NPions", fHistoTitle.c_str(), 15, 0., 15. ) );
   fHisto.push_back( new TH1F( "NChargesPions", fHistoTitle.c_str(), 15, 0., 15. ) ) ;
   fHisto.push_back( new TH1F( "NPi0s", fHistoTitle.c_str(), 15, 0., 15. ) );
   fHisto.push_back( new TH1F( "NGammas", fHistoTitle.c_str(), 15, 0., 15. ) );
   fHisto.push_back( new TH1F( "NNeutrons", fHistoTitle.c_str(), 25, 0., 25. ) );
   fHisto.push_back( new TH1F( "ChargedSecondaryMomentum", fHistoTitle.c_str(), 50, 0., 1. ) );
   fHisto.push_back( new TH1F( "ChargedPionMomentum", fHistoTitle.c_str(), 50, 0., 1. ) );
   fHisto.push_back( new TH1F( "PDGIDSecondaries", fHistoTitle.c_str(), 5000, 0., 5000. ) );
   fHisto.push_back( new TH1F( "ChargeOfSecondary", fHistoTitle.c_str(), 30, -15., 15. ) ) ;
   
   return;
   
}

void TestStoppingHisto::FillEvt( G4VParticleChange* aChange )
{

   G4int NSec = aChange->GetNumberOfSecondaries();
   
   const G4DynamicParticle* sec = 0;
   // G4ParticleDefinition* pd = 0;
   
   int NChSec    = 0;
   int NPions    = 0;
   int NChPions  = 0;
   int NGammas   = 0;
   int NKaons    = 0;
   int NNeutrons = 0;
      
   fHisto[0]->Fill( (double)NSec );
   
   for (G4int i=0; i<NSec; i++) 
   {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();
	
	fHisto[10]->Fill( (float)(sec->GetCharge()) );
	
	if ( sec->GetCharge() != 0 ) 
	{
	   NChSec++;
	   fHisto[7]->Fill( (sec->GetTotalMomentum()/GeV) );
	}
	
	const G4String& pname = sec->GetDefinition()->GetParticleName();

	float pdgid = (float)(sec->GetPDGcode());
	if ( pdgid > 4999. ) pdgid = 4999.;
	fHisto[9]->Fill( pdgid );
	
	
	if ( pname == "pi-" || pname == "pi+" || pname == "pi0" )
	{
	   NPions++;
	   if ( pname != "pi0" )
	   {
	      fHisto[8]->Fill( (sec->GetTotalMomentum()/GeV) );
	   }
	}
	
	if ( pname == "gamma" || pname == "photon" )
	{
	   NGammas++;
	}
	
	if ( pname == "neutron" )
	{
	   NNeutrons++;
	   if ( fBeam == "pi-" )
	   {
	      fHisto[11]->Fill( (sec->GetKineticEnergy()/MeV) );
	   }
	}
	
	if ( pname == "kaon+" || pname == "kaon-" || pname == "kaon0S" || pname == "kaon0L" )
	{
	   NKaons++;
	   if ( pname == "kaon+" || pname == "kaon-" )
	   {
	      fHisto[10]->Fill( (sec->GetTotalMomentum()/GeV) );
	   }
	}	
   }
   
   if ( NChSec > 0 )            fHisto[1]->Fill( (double)NChSec );
   if ( NPions > 0 )            fHisto[2]->Fill( (double)NPions );
   if ( NChPions > 0 )          fHisto[3]->Fill( (double)NChPions );
   if ( (NPions-NChPions) > 0 ) fHisto[4]->Fill( (double)(NPions-NChPions) );
   if ( NGammas > 0 )           fHisto[5]->Fill( (double)NGammas ); 
   if ( NNeutrons > 0 )         fHisto[6]->Fill( (double)NNeutrons );
   
   return;
   
}

void TestStoppingHisto::Write( int stat )
{

   std::string ptag;
   if ( fBeam == "pi-" )
   {
      ptag = "piminus";
   }
   else if ( fBeam == "kaon-" || fBeam == "Kaon-" )
   {
      ptag = "kminus";
   }
   else if ( fBeam == "anti_proton" )
   {
      ptag = "antiproton";
   }
   
   std::string fname = ptag + fTarget + fModel;
   fname += ".root";
   TFile f( fname.c_str(), "recreate" );
   for ( size_t i=0; i<fHisto.size(); i++ )
   {
      double xbin = fHisto[i]->GetBinWidth(1);
      double scale = 1. / ((double)stat * xbin) ;
      fHisto[i]->Scale(scale);
      fHisto[i]->Write();
   }
   f.Close();
   return;

}

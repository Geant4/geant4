
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
   
   if ( fBeam == "pi-" )
   {
      InitPionMinus();
   }
   else if ( fBeam == "anti_proton" )
   {
      InitAntiProton();
   }
   
   return;

}

void TestStoppingHisto::InitPionMinus()
{

   std::string title = fBeam + " on " + fTarget;
   title += "; Neutron Yield vs Kinetic Energy";
   
   fHisto.push_back(new TH1F( "NvsT", title.c_str(), 600, 1., 151.));

   return;

}

void TestStoppingHisto::InitAntiProton()
{

   std::string title = fBeam + " on " + fTarget;
   
   fHisto.push_back( new TH1F( "ChargedPionMomentum", title.c_str(), 50, 0., 1.) );
   fHisto.push_back( new TH1F( "PionMultiplicity", title.c_str(), 10, 0., 10. ) );

   return;

}


void TestStoppingHisto::FillEvt( G4VParticleChange* aChange )
{

   G4int NSec = aChange->GetNumberOfSecondaries();
   
   const G4DynamicParticle* sec = 0;
   G4ParticleDefinition* pd = 0;
   
   int NPions = 0;
      
   for (G4int i=0; i<NSec; i++) 
   {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();
	const G4String& pname = sec->GetDefinition()->GetParticleName();
	if ( fBeam == "pi-" )
	{
	   if ( pname == "neutron" )
	   {
	      fHisto[0]->Fill( (sec->GetKineticEnergy()/MeV) );
	   }
	}
	else if ( fBeam == "anti_proton" )
	{
	   if ( pname == "pi-" || pname == "pi+" )
	   {
	   // FillEvtAntiProton( sec );
	      fHisto[0]->Fill( (sec->GetTotalMomentum()/GeV) );
	      NPions++;
	   }
	   else if ( pname == "pi0" )
	   {
	      NPions++;
	   }
	}
   }
   
   if ( fBeam == "anti_proton" )
   {
      fHisto[1]->Fill( (double)NPions );
   }
   
   return;
   
}

/*
void TestStoppingHisto::FillEvtAntiProton( const G4DynamicParticle* sec )
//void TestStoppingHisto::FillEvtAntiProton( G4VParticleChange* aChange )
{
   
   const G4String& pname = sec->GetDefinition()->GetParticleName();
   if ( pname == "pi-" || pname == "pi+" )
   {
      // G4ThreeVector pmom = sec->GetMomentum();
      double ptot = sec->GetTotalMomentum()/GeV;
      fHisto[0]->Fill(ptot);
   }

   return;
}
*/

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

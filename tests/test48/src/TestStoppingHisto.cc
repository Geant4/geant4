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

#include "TestStoppingHisto.hh"

#include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"

TestStoppingHisto::~TestStoppingHisto()
{

   for (size_t i=0; i<fHisto.size(); i++ )
   {
      delete fHisto[i];
   }

   for (size_t i=0; i<fMuHisto.size(); i++ )
   {
      delete fMuHisto[i];
   }

}

void TestStoppingHisto::Init()
{
  
  TH1::SetDefaultSumw2();

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
  else if ( fBeam == "mu-" )
    {
      ptag = "muminus";
    } 
  else {
    G4cout << "You may need to extend TestStoppingHisto::Write to include beam of: " << fBeam
           << G4endl;
    ptag = fBeam;

  }
   
  fHistoTitle = fBeam + " on " + fTarget;

  G4cout << "Initializing histograms for: " << fHistoTitle
         << G4endl;

  if ( fBeam == "sigma-" )
    {
      InitSigmaMinus();
      return;
    }
   
  if ( fBeam == "kaon-" )
    {
      InitKaonMinus();
      return;
    }   
  if ( fBeam == "mu-" )
    {
      InitMuonMinus();
      return;
    }
   
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

void TestStoppingHisto::InitSigmaMinus()
{

   fHisto.push_back( new TH1F( "NSecondaries", fHistoTitle.c_str(), 50, 0., 50. ) );
   fHisto.push_back( new TH1F( "PDGIDSecondaries", fHistoTitle.c_str(), 10000, -5000., 5000. ) );
   // 0 - other
   // 1 - gamma
   // 2 - pi0
   // 3 - pi+
   // 4 - pi-
   // 5 - p
   // 6 - n
   // 7 - Lambda
   // 8 - Sigma-
   // 9 - Sigma0
   // 10 - Sigma+
   //
   TH1F* ptmult = new TH1F( "PartTypeMult", fHistoTitle.c_str(), 15, 0., 15. );
   TAxis* xaxis1 = ptmult->GetXaxis(); 
   xaxis1->SetBinLabel( 1, "other" );
   xaxis1->SetBinLabel( 2, "#gamma" );
   xaxis1->SetBinLabel( 3, "#pi^{0}" );
   xaxis1->SetBinLabel( 4, "#pi^{+}" );
   xaxis1->SetBinLabel( 5, "#pi^{-}" );
   xaxis1->SetBinLabel( 6, "p" );
   xaxis1->SetBinLabel( 7, "n" );
   xaxis1->SetBinLabel( 8, "#Lambda" );
   xaxis1->SetBinLabel( 9, "#Sigma^{-}" );
   xaxis1->SetBinLabel( 10, "#Sigma^{0}" );
   xaxis1->SetBinLabel( 11, "#Sigma^{+}" );
   fHisto.push_back( ptmult );
   
   return;

}

void TestStoppingHisto::InitKaonMinus()
{

   fHisto.push_back( new TH1F( "NSecondaries", fHistoTitle.c_str(), 50, 0., 50. ) );
   TH1F* topology = new TH1F( "Topology", fHistoTitle.c_str(), 10, 0., 10. );
   // 0 - unknown (other)
   // 1 - Lambda pi0
   // 2 - Lambda gamma pi0 (also as Sigma pi0 ???)
   // 3 - Sigma+ pi-
   // 4 - Sigma- pi+
   // 5 - Lambda pi+ pi-
   // 6 - Lambda pi0 pi0
   // 7 - Sigma0 gamma
   // 8 - Lambda gamma 
   // 9 - Sigma pi0
   TAxis* xaxis = topology->GetXaxis();
   xaxis->SetBinLabel( 1, "other" );
   xaxis->SetBinLabel( 2, "#Lambda #pi^{0}" );
   xaxis->SetBinLabel( 3, "#Lambda #gamma #pi^{0}" );
   xaxis->SetBinLabel( 4, "#Sigma^{+} #pi^{-}" );
   xaxis->SetBinLabel( 5, "#Sigma^{-} #pi^{+}" );
   xaxis->SetBinLabel( 6, "#Lambda #pi^{+} #pi^{-}" );
   xaxis->SetBinLabel( 7, "#Lambda #pi^{0} #pi^{0}" );
   xaxis->SetBinLabel( 8, "#Sigma^{0} #gamma" );
   xaxis->SetBinLabel( 9, "#Lambda #gamma" );
   xaxis->SetBinLabel( 10, "#Sigma^{0} #pi^{0}" );
   fHisto.push_back( topology );
   
   fHisto.push_back( new TH1F( "PDGIDSecondaries", fHistoTitle.c_str(), 10000, -5000., 5000. ) );
    
   std::string title;
   title = fHistoTitle + "; Energy Of Secondary Pi0 (MeV)";   
   fHisto.push_back( new TH1F( "EnergyPi0", title.c_str(), 125, 100., 350. ) );
   title = fHistoTitle + "; Energy of Secondary Pi-/+ (MeV)";
   fHisto.push_back( new TH1F( "EnergyChargedPions", title.c_str(), 125, 100., 350. ) );
   title = fHistoTitle +"; Energy of Secondary Photon (MeV)";
   fHisto.push_back( new TH1F( "EnergyPhoton", title.c_str(), 500, 0., 500. ) );
   
   // 0 - other
   // 1 - gamma
   // 2 - pi0
   // 3 - pi+
   // 4 - pi-
   // 5 - p
   // 6 - n
   // 7 - Lambda
   // 8 - Sigma-
   // 9 - Sigma0
   // 10 - Sigma+
   //
   TH1F* ptmult = new TH1F( "PartTypeMult", fHistoTitle.c_str(), 15, 0., 15. );
   TAxis* xaxis1 = ptmult->GetXaxis(); 
   xaxis1->SetBinLabel( 1, "other" );
   xaxis1->SetBinLabel( 2, "#gamma" );
   xaxis1->SetBinLabel( 3, "#pi^{0}" );
   xaxis1->SetBinLabel( 4, "#pi^{+}" );
   xaxis1->SetBinLabel( 5, "#pi^{-}" );
   xaxis1->SetBinLabel( 6, "p" );
   xaxis1->SetBinLabel( 7, "n" );
   xaxis1->SetBinLabel( 8, "#Lambda" );
   xaxis1->SetBinLabel( 9, "#Sigma^{-}" );
   xaxis1->SetBinLabel( 10, "#Sigma^{0}" );
   xaxis1->SetBinLabel( 11, "#Sigma^{+}" );
   fHisto.push_back( ptmult );
   
   //fHisto.push_back( new TH1F( "NPi0s", fHistoTitle.c_str(), 15, 0., 15. ) );
   //fHisto.push_back( new TH1F( "NGammas", fHistoTitle.c_str(), 15, 0., 15. ) );

   return;

}

void TestStoppingHisto::InitHistoGeneral()
{
   
  fHisto.push_back( new TH1F( "NSecondaries", fHistoTitle.c_str(), 25, 0., 25. ) );       // 0
  fHisto.push_back( new TH1F( "NChargedSecondaries", fHistoTitle.c_str(), 25, 0., 25. ) );// 1

   double bins_npions[] = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5 };
   int    nbins_npions = sizeof(bins_npions) / sizeof(double) - 1;
    
//  fHisto.push_back( new TH1F( "NPions", fHistoTitle.c_str(), 15, 0., 15. ) );             // 2
  fHisto.push_back( new TH1F( "NPions", fHistoTitle.c_str(), nbins_npions, bins_npions ) );             // 2
    
  fHisto.push_back( new TH1F( "NChargesPions", fHistoTitle.c_str(), 15, 0., 15. ) ) ;     // 3
  fHisto.push_back( new TH1F( "NPi0s", fHistoTitle.c_str(), 15, 0., 15. ) );              // 4
  fHisto.push_back( new TH1F( "NGammas", fHistoTitle.c_str(), 15, 0., 15. ) );            // 5
  fHisto.push_back( new TH1F( "NNeutrons", fHistoTitle.c_str(), 10, -1., 9. ) );          // 6
  fHisto.push_back( new TH1F( "ChargedSecondaryMomentum", fHistoTitle.c_str(), 50,0.,1.));// 7

  // this is actually the binning for the pbar+H
  // we'll probably need to factor out the pbar part into a separate utility
  //
  //... although it should be the same as 50 bins from 0. to 1., the "standard way"
  //
  int nbins_mom_pions = 50;
  double bins_mom_pions[51] ;
  for ( int i=0; i<51; ++i )
  {
      bins_mom_pions[i] = 0.02*i;
  }
  // fHisto.push_back( new TH1F( "ChargedPionMomentum", fHistoTitle.c_str(), 50, 0., 1. ) ); // 8
  fHisto.push_back( new TH1F( "ChargedPionMomentum", fHistoTitle.c_str(), nbins_mom_pions, bins_mom_pions ) ); // 8
  
  fHisto.push_back( new TH1F( "PDGIDSecondaries", fHistoTitle.c_str(), 5000, 0., 5000.)); // 9;
  fHisto.push_back( new TH1F( "ChargeOfSecondary", fHistoTitle.c_str(), 30, -15., 15.));  //10;
  fHisto.push_back( new TH1F( "NElectrons",  fHistoTitle.c_str(), 15, 0., 15. ) );        //11
  fHisto.push_back( new TH1F( "NMuons",      fHistoTitle.c_str(), 15, 0., 15. ) );        //12
  fHisto.push_back( new TH1F( "NNeutrinoEs", fHistoTitle.c_str(), 15, 0., 15. ) );        //13
  fHisto.push_back( new TH1F( "NNeutrinoMs", fHistoTitle.c_str(), 15, 0., 15. ) );        //14
  fHisto.push_back( new TH1F( "NeutronMomentum", fHistoTitle.c_str(), 50, 0., 1. ) );     //15
  // note that InitPionMinus is called afterwards, so ... one needs to update filling neutron kin en
  // FIXME

  return;   

}

void TestStoppingHisto::InitMuonMinus() {
   
  fHisto.push_back( new TH1F( "NSecondaries", fHistoTitle.c_str(), 25, 0., 25. ) );        // 0
  fHisto.push_back( new TH1F( "NChargedSecondaries", fHistoTitle.c_str(), 25, 0., 25. ) ); // 1
  fHisto.push_back( new TH1F( "NPions", fHistoTitle.c_str(), 15, 0., 15. ) );              // 2
  fHisto.push_back( new TH1F( "NChargesPions", fHistoTitle.c_str(), 15, 0., 15. ) );       // 3
  fHisto.push_back( new TH1F( "NPi0s", fHistoTitle.c_str(), 15, 0., 15. ) );               // 4
  fHisto.push_back( new TH1F( "NGammas", fHistoTitle.c_str(), 15, 0., 15. ) );             // 5
  fHisto.push_back( new TH1F( "NNeutrons", fHistoTitle.c_str(), 25, 0., 25. ) );           // 6
  fHisto.push_back( new TH1F( "ChargedSecondaryMomentum", fHistoTitle.c_str(), 50,0.,1.)); // 7
  fHisto.push_back( new TH1F( "ChargedPionMomentum", fHistoTitle.c_str(), 50, 0., 1. ) );  // 8
  fHisto.push_back( new TH1F( "PDGIDSecondaries", fHistoTitle.c_str(), 5000, 0., 5000. ) );// 9
  fHisto.push_back( new TH1F( "ChargeOfSecondary", fHistoTitle.c_str(), 30, -15., 15. ) ) ;//10
  fHisto.push_back( new TH1F( "NElectrons",  fHistoTitle.c_str(), 15, 0., 15. ) );         //11
  fHisto.push_back( new TH1F( "NMuons",      fHistoTitle.c_str(), 15, 0., 15. ) );         //12
  fHisto.push_back( new TH1F( "NNeutrinoEs", fHistoTitle.c_str(), 15, 0., 15. ) );         //13
  fHisto.push_back( new TH1F( "NNeutrinoMs", fHistoTitle.c_str(), 15, 0., 15. ) );         //14
  fHisto.push_back( new TH1F( "NeutronMomentum", fHistoTitle.c_str(), 50, 0., 1. ) );      //15
  fHisto.push_back( new TH1F( "NeutronKineticEnergy", fHistoTitle.c_str(), 50, 0., 100.)); //16
  fHisto.push_back( new TH1F( "NAlphas", fHistoTitle.c_str(), 25, 0., 25. ) );             //17
  fHisto.push_back( new TH1F( "AlphaMomentum", fHistoTitle.c_str(), 50, 0., 1. ) );        //18
   
  return;
   
}

void TestStoppingHisto::FillEvt( G4VParticleChange* aChange ) 
{

  G4int NSec = aChange->GetNumberOfSecondaries();
  
   if ( fBeam == "sigma-" )
   {
      FillEvtSigmaBeam( aChange );
      return;
   }
   
   if ( fBeam == "kaon-" )
   {
      FillEvtKaonBeam( aChange );
      return;
   }
   
   if ( fBeam == "mu-" )
   {
     G4TrackVector secondaries;
     for (G4int i=0; i<NSec; i++) {
       secondaries.push_back(aChange->GetSecondary(i));
     }
     FillEvtMuonMinusBeam( &secondaries );
     return;
   }
   
   const G4DynamicParticle* sec = 0;
   
   int NChSec     = 0;
   int NPions     = 0;
   int NChPions   = 0;
   int NPi0s      = 0;
   int NGammas    = 0;
   int NKaons     = 0;
   int NNeutrons  = 0;

   int NElectrons = 0;
   int NMuons     = 0;

   int NNeutrinoEs = 0;
   int NNeutrinoMs = 0;
      
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

        G4int pdgiid = sec->GetPDGcode();

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
	   else
	   {
	      NPi0s++;
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
	      fHisto[16]->Fill( (sec->GetKineticEnergy()/MeV) );
	   }
           fHisto[15]->Fill( (sec->GetTotalMomentum()/GeV) );
	}
	
	if ( pname == "kaon+" || pname == "kaon-" || pname == "kaon0S" || pname == "kaon0L" )
	{
	   NKaons++;
	   if ( pname == "kaon+" || pname == "kaon-" )
	   {
	      fHisto[10]->Fill( (sec->GetTotalMomentum()/GeV) );
	   }
	}	

	if ( std::abs(pdgiid) == 11 ) 
	{
	   NElectrons++;
	}

 	if ( std::abs(pdgiid) == 13 ) 
	{
	   NMuons++;
	}

	if ( std::abs(pdgiid) == 12 ) 
	{
	   NNeutrinoEs++;
	}

	if ( std::abs(pdgiid) == 14 ) 
	{
	   NNeutrinoMs++;
	}

   }
   
   if ( NChSec > 0 )            fHisto[1]->Fill( (double)NChSec );
   if ( NPions > 0 )            fHisto[2]->Fill( (double)NPions );
   if ( NChPions > 0 )          fHisto[3]->Fill( (double)NChPions );
   if ( (NPions-NChPions) > 0 ) fHisto[4]->Fill( (double)(NPions-NChPions) );
   if ( NGammas > 0 )           fHisto[5]->Fill( (double)NGammas ); 
   if ( NNeutrons > 0 )         fHisto[6]->Fill( (double)NNeutrons );
   if ( NElectrons  > 0 )        fHisto[11]->Fill( (double)NElectrons );
   if ( NMuons      > 0 )        fHisto[12]->Fill( (double)NMuons );
   if ( NNeutrinoEs > 0 )        fHisto[13]->Fill( (double)NNeutrinoEs );
   if ( NNeutrinoMs > 0 )        fHisto[14]->Fill( (double)NNeutrinoMs );
   
   return;
   
}

void TestStoppingHisto::FillEvtSigmaBeam( G4VParticleChange* aChange )
{

   G4int NSec = aChange->GetNumberOfSecondaries();
   fHisto[0]->Fill( (double)NSec );
   
   std::vector<int> SpeciesMult;
   for ( int i=0; i<11; i++ )
   {
      SpeciesMult.push_back(0);
   }

   const G4DynamicParticle* sec = 0;
   
   for ( G4int ip=0; ip<NSec; ip++ )
   {
      sec = aChange->GetSecondary(ip)->GetDynamicParticle();
      int pdgid = sec->GetPDGcode();
      fHisto[1]->Fill( (float)pdgid );
      if ( pdgid == 22 )
      {
//	 for ( G4int i1=0; i1<NSec; i1++ )
//	 {
//	    const G4DynamicParticle* sec1 = aChange->GetSecondary(i1)->GetDynamicParticle();
//	    sec1->DumpInfo();
//	 }
         SpeciesMult[1]++;
      }
      else if ( pdgid == 111 )
      {
         SpeciesMult[2]++;
      }
      else if ( pdgid == 211 )
      {
         SpeciesMult[3]++;
      }
      else if ( pdgid == -211 )
      {
         SpeciesMult[4]++;
      }
      else if ( pdgid == 2212 )
      {
         SpeciesMult[5]++;
      }
      else if ( pdgid == 2112 )
      {
         SpeciesMult[6]++;
      }
      else if ( pdgid == 3122 )
      {
         SpeciesMult[7]++;
      }
      else if ( pdgid == 3112 )
      {
         SpeciesMult[8]++;
      }
      else if ( pdgid == 3212 )
      {
         SpeciesMult[9]++;
      }
      else if ( pdgid == 3222 )
      {
         SpeciesMult[10]++;
      }
      else
      {
         SpeciesMult[0]++;
      }
   }

   for ( unsigned int it=0; it<SpeciesMult.size(); it++ )
   {
      fHisto[2]->Fill( (double)it, SpeciesMult[it] );
   }

   return;

}


void TestStoppingHisto::FillEvtKaonBeam( G4VParticleChange* aChange )
{

   G4int NSec = aChange->GetNumberOfSecondaries();
   fHisto[0]->Fill( (double)NSec );
   
   // Topology channels (maybe not all...)
   //
   // these topologies below seem to be valid for K- and K- on D
   //
   // 1 - Lambda pi0
   // 2 - Lambda gamma pi0 (also as Sigma pi0 ???)
   // 3 - Sigma+ pi-
   // 4 - Sigma- pi+
   // 5 - Lambda pi+ pi-
   // 6 - Lambda pi0 pi0
   // 7 - Sigma0 gamma
   // 8 - Lambda gamma 
   // 9 - Sigma pi0
   //
   // however, starting already He, topologies become mostly p/n and sometimes a pi-
   // (very popular: 2p+2n+pi- or p+3n)
   //
   // 11 - p+3n
   // 12 - 2p+2n+pi-
   
   // PDG ID's: 22 - gamma, 111 - pi0, 211 - pi+ 
   //           2212 - p, 2112 - n, 3122 - Lambda, 3112 - Sigma-, 3212 - Sigma0, 3222 - Sigma+ 
   
   // 0 - other
   // 1 - gamma
   // 2 - pi0
   // 3 - pi+
   // 4 - pi-
   // 5 - p
   // 6 - n
   // 7 - Lambda
   // 8 - Sigma-
   // 9 - Sigma0
   // 10 - Sigma+
   //
   std::vector<int> SpeciesMult;
   for ( int i=0; i<11; i++ )
   {
      SpeciesMult.push_back(0);
   }
   
   
   int channel = Topology( aChange );      
   
   fHisto[1]->Fill( (float)channel ); 

   const G4DynamicParticle* sec = 0;
   for ( G4int ip=0; ip<NSec; ip++ )
   {
      sec = aChange->GetSecondary(ip)->GetDynamicParticle();
      
//      if ( channel == 0 )
//      {
//         sec->DumpInfo();
//      }
      
      int pdgid = sec->GetPDGcode();
      if ( pdgid == 22 )
      {
//	 for ( G4int i1=0; i1<NSec; i1++ )
//	 {
//	    const G4DynamicParticle* sec1 = aChange->GetSecondary(i1)->GetDynamicParticle();
//	    sec1->DumpInfo();
//	 }
         SpeciesMult[1]++;
      }
      else if ( pdgid == 111 )
      {
         SpeciesMult[2]++;
      }
      else if ( pdgid == 211 )
      {
         SpeciesMult[3]++;
      }
      else if ( pdgid == -211 )
      {
         SpeciesMult[4]++;
      }
      else if ( pdgid == 2212 )
      {
         SpeciesMult[5]++;
      }
      else if ( pdgid == 2112 )
      {
         SpeciesMult[6]++;
      }
      else if ( pdgid == 3122 )
      {
         SpeciesMult[7]++;
      }
      else if ( pdgid == 3112 )
      {
         SpeciesMult[8]++;
      }
      else if ( pdgid == 3212 )
      {
         SpeciesMult[9]++;
      }
      else if ( pdgid == 3222 )
      {
         SpeciesMult[10]++;
      }
      else
      {
         SpeciesMult[0]++;
      }
           
      if ( pdgid == 111 )
      {
         fHisto[3]->Fill( sec->GetTotalEnergy()/MeV ) ;
      }
      else if ( std::abs(pdgid) == 211 )
      {
         fHisto[4]->Fill( sec->GetTotalEnergy()/MeV );
      }
      else if ( pdgid == 22 )
      {
         fHisto[5]->Fill( sec->GetTotalEnergy()/MeV );
      }
   }
   
   for ( unsigned int it=0; it<SpeciesMult.size(); it++ )
   {
      fHisto[6]->Fill( (double)it, SpeciesMult[it] );
   }
   
   return;

}

int TestStoppingHisto::Topology( G4VParticleChange* aChange )
{
   
   int channel = 0;
   
   G4int NSec = aChange->GetNumberOfSecondaries();
   
   const G4DynamicParticle* sec = 0;

   // find 1st hadron (since there might be leptons from EMCascade
   int firstHadronPos = -1;
   if ( NSec != 2 ) // this is an ugly hack, to deal with the fact that "traditional" code
                    // will only capture as 2-body, no leptons/gammas from EM cascade, but
		    // one of the 2-body modes will be Lambda gamma, so the "lepton/gamma skipp"
		    // will chock on that 
   {
      for ( G4int is=0; is<NSec; is++ )
      {      
         const G4String& ptype = aChange->GetSecondary(is)->GetParticleDefinition()->GetParticleType();
         if ( ptype == "meson" || ptype == "baryon" ) 
         {
            firstHadronPos = is;
            break;
         }     
      }
   }
   else
   {
      firstHadronPos = 0;
   }
   
   if ( firstHadronPos < 0 ) return channel; // nothing sensible
   
   if ( (NSec - firstHadronPos) < 2 ) return channel; // it should be at least 2 secondaries lest at the end of the list  
      
   for ( G4int ip=firstHadronPos; ip<NSec; ip++ )
   {
      sec = aChange->GetSecondary(ip)->GetDynamicParticle();
      int pdgid = sec->GetPDGcode();
      fHisto[2]->Fill( (float)pdgid );
      if ( pdgid == 3122 ) // found Lambda !!!                    
      {
	 if ( (NSec-firstHadronPos) == 2 )
	 {
	    // it's either Lambda gamma or Lambda pi0 )
	    for ( G4int ip1 = firstHadronPos; ip1<NSec; ip1++ )
	    {
	       if ( ip1 == ip ) continue;
	       int pdgid1 = aChange->GetSecondary(ip1)->GetDynamicParticle()->GetPDGcode();
	       if ( pdgid1 == 22 )
	       {
	          channel = 8;
		  break;
	       }
	       else if ( pdgid1 == 111 )
	       {
	          channel = 1;
		  break;
	       }
	    }
	 }
	 
	 if ( channel > 0 ) break;
	 
	 if ( (NSec-firstHadronPos) == 3 )
	 {
	    // it's either Lambda pi- pi+, or Lambda pi0 pi0, or Lambda pi0 gamma
	    for ( G4int ip2=firstHadronPos; ip2<NSec; ip2++ )
	    {
	       if ( ip2 == ip ) continue;
	       int pdgid2 = aChange->GetSecondary(ip2)->GetDynamicParticle()->GetPDGcode();
	       if ( std::abs(pdgid2) == 211 )
	       {
	          channel = 5; // it's Lambda pi- pi+
		  break;
	       }
	       else if ( pdgid2 == 111 )
	       {
	          for ( G4int ip3=firstHadronPos; ip3<NSec; ip3++ )
		  {
		     if ( ip3 == ip || ip3 == ip2 ) continue;
		     int pdgid3 = aChange->GetSecondary(ip3)->GetDynamicParticle()->GetPDGcode();
		     if ( pdgid3 == 22 )
		     {
		        channel = 2;
			break;
		     }
		     else if ( pdgid3 == 111 )
		     {
		        channel = 6;
			break;
		     }
		  }
		  if ( channel > 0 ) break;
	       }
	    }
	    
	    if ( channel > 0 ) break;
	    
	 } 
      }
      
      if ( channel > 0 ) break;
      
      // 3112 - Sigma-, 3212 - Sigma0, 3222 - Sigma+
      if ( pdgid == 3112 ) // Sigma-
      {
         if ( (NSec-firstHadronPos) == 2 )
	 {
	    channel = 4;
	 }
      }
      else if ( pdgid == 3222 ) // Sigma+
      {
         if ( (NSec-firstHadronPos) == 2 )
	 {
	    channel = 3;
	 }
      }
      else if ( pdgid == 3212 ) // Sigma0
      {
         for ( G4int ip4=firstHadronPos; ip4<NSec; ip4++ )
	 {
	    if ( ip4 == ip ) continue;
	    int pdgid4 = aChange->GetSecondary(ip4)->GetDynamicParticle()->GetPDGcode();
	    if ( pdgid4 == 111 )
	    {
	       channel = 9;
	       break;
	    }
	    else if ( pdgid4 == 22 )
	    {
	       channel = 7;
	       break;
	    }
	 }
      }
       
      if ( channel > 0 ) break; 
            
   } 
   
   if ( channel > 0 ) return channel;
   
   // OK, now try for "unknown" topology and/or for heavier targets
     
   
   return channel;
   
}


void TestStoppingHisto::FillEvtMuonMinusBeam(G4TrackVector *secondaries )
{

  // for muon beam only for now
   
  G4int NSec = secondaries->size();
   
  int NChSec     = 0;
  int NPions     = 0;
  int NChPions   = 0;
  int NPi0s      = 0;
  int NGammas    = 0;
  int NKaons     = 0;
  int NNeutrons  = 0;

  int NElectrons = 0;
  int NMuons     = 0;

  int NNeutrinoEs = 0;
  int NNeutrinoMs = 0;

  int NAlphas     = 0;
      
  fHisto[0]->Fill( (double)NSec );
   
  for (G4int i=0; i<NSec; i++) {
    const G4DynamicParticle *sec = ((*secondaries)[i])->GetDynamicParticle();
	
    fHisto[10]->Fill( (float)(sec->GetCharge()) );
	
    if ( sec->GetCharge() != 0 ) 
      {
        NChSec++;
        fHisto[7]->Fill( (sec->GetTotalMomentum()/GeV) );
      }
	
    const G4String& pname = sec->GetDefinition()->GetParticleName();

    G4int pdgiid = sec->GetPDGcode();

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
        else
          {
            NPi0s++;
          }
      }
	
    if ( pname == "gamma" || pname == "photon" )
      {
        NGammas++;
      }
	
    if ( pname == "neutron" )
      {
        NNeutrons++;
        fHisto[15]->Fill( (sec->GetTotalMomentum()/GeV) );
        fHisto[16]->Fill( (sec->GetKineticEnergy()/MeV) );
      }
	
    if ( pname == "alpha" )
      {
        NAlphas++;
        fHisto[18]->Fill( (sec->GetTotalMomentum()/GeV) );
      }
	
    if ( pname == "kaon+" || pname == "kaon-" || pname == "kaon0S" || pname == "kaon0L" )
      {
        NKaons++;
        if ( pname == "kaon+" || pname == "kaon-" )
          {
            fHisto[10]->Fill( (sec->GetTotalMomentum()/GeV) );
          }
      }	

    if ( std::abs(pdgiid) == 11 ) 
      {
        NElectrons++;
      }

    if ( std::abs(pdgiid) == 13 ) 
      {
        NMuons++;
      }

    if ( std::abs(pdgiid) == 12 ) 
      {
        NNeutrinoEs++;
      }

    if ( std::abs(pdgiid) == 14 ) 
      {
        NNeutrinoMs++;
      }

  }
   
  fHisto[ 1]->Fill( (double)NChSec );
  fHisto[ 2]->Fill( (double)NPions );
  fHisto[ 3]->Fill( (double)NChPions );
  fHisto[ 4]->Fill( (double)(NPions-NChPions) );
  fHisto[ 5]->Fill( (double)NGammas ); 
  fHisto[ 6]->Fill( (double)NNeutrons );
  fHisto[11]->Fill( (double)NElectrons );
  fHisto[12]->Fill( (double)NMuons );
  fHisto[13]->Fill( (double)NNeutrinoEs );
  fHisto[14]->Fill( (double)NNeutrinoMs );
  fHisto[17]->Fill( (double)NAlphas );
   
  return;
   
}

void TestStoppingHisto::Write( int stat )
{

   std::string fname = ptag + fTarget + fModel;
   if ( fJobID > -1 )
   {
      char buf[5];
      sprintf( buf, "%i", fJobID );
      fname += "-";
      fname.append( buf ); 
   }  
   fname += ".root";

   G4cout << "Writing histogram file: " << fname
          << G4endl;

   TFile f( fname.c_str(), "recreate" );
   for ( size_t i=0; i<fHisto.size(); ++i )
   {
      double xbin = fHisto[i]->GetBinWidth(1);
      double scale = 1. / ((double)stat * xbin) ;
      fHisto[i]->Scale(scale);
      fHisto[i]->Write();
   }
   for ( size_t i=0; i<fMuHisto.size(); ++i )
   {
      double xbin = fMuHisto[i]->GetBinWidth(1);
      double scale = 1. / ((double)stat * xbin) ;
      fMuHisto[i]->Scale(scale);
      fMuHisto[i]->Write();
   }
   f.Close();
   return;

}

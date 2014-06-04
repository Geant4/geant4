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

#include "G4TrackVector.hh"

#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"

// #include "G4Proton.hh" // tmp needed for Ecm/beta/gamma cross-checks

#include "G4SystemOfUnits.hh"

Tst23NA49Histo::Tst23NA49Histo( G4String htitle )
   : TstHistoSet( htitle ), TstHistoSetForNu()
{

  G4String title;
  G4String outcome;
  G4String hname;
  G4String pim_names[] = { "piminus_dNdxF_", "piminus_pT_" };
  G4String pip_names[] = { "piplus_dNdxF_",  "piplus_pT_" };
  
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
  
// FIXME !!!
//
  fHistoPTvsXFProton = 0; // for now...
  
       
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
  
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dNdxF",    title.c_str(), npbarbins, pbarbins ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT",   title.c_str(), npbarbins, pbarbins, 0., 10. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT2",  title.c_str(), npbarbins, pbarbins, 0., 10. ) );
  
  fHistoPTvsXFAntiProton = 0;	

  fHistoSecAntiProtonTot.push_back( new TH1D( "antiproton_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProtonTot.push_back( new TH1D( "antiproton_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecAntiProtonTot.push_back( new TProfile( "antiproton_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecAntiProtonTot.push_back( new TProfile( "antiproton_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	
  
  outcome = " -> X + pi- ";  
  title = htitle + outcome;

//  double pionbins[] = { -0.500, -0.400, -0.300, -0.250, -0.200, -0.150, -0.125, -0.100, -0.075, 
//                     -0.050, -0.040, -0.030, -0.020, -0.010,  
//		      0.0,
//		      0.010,  0.020,  0.030,  0.040,  0.050,
//		      0.075,  0.100,  0.125,  0.150,  0.200,  0.250,  0.300,  0.400,  0.500 }; 
		       
  // the idea is to make the bin center correspond to the number in the NA49 data file(s)
  //
  fNPiBinsXF = 30;
  fPiBinsXF = new double[fNPiBinsXF];  
  fPiBinsXF[0] = -0.550;
  fPiBinsXF[1] = -0.450;
  fPiBinsXF[2] = -0.350; 
  fPiBinsXF[3] = -0.275;
  fPiBinsXF[4] = -0.225;
  fPiBinsXF[5] = -0.175;
  fPiBinsXF[6] = -0.1375;
  fPiBinsXF[7] = -0.1125;
  fPiBinsXF[8] = -0.0875; 
  fPiBinsXF[9] = -0.0625;
  fPiBinsXF[10] = -0.045;
  fPiBinsXF[11] = -0.035;
  fPiBinsXF[12] = -0.025;
  fPiBinsXF[13] = -0.015; 
  fPiBinsXF[14] = -0.005;
  fPiBinsXF[15] =  0.005;  
  fPiBinsXF[16] =  0.015;
  fPiBinsXF[17] =  0.025;
  fPiBinsXF[18] =  0.035;  
  fPiBinsXF[19] =  0.045;
  fPiBinsXF[20] =  0.0625;  
  fPiBinsXF[21] =  0.0875;
  fPiBinsXF[22] =  0.1125;
  fPiBinsXF[23] =  0.1375;
  fPiBinsXF[24] =  0.175; 
  fPiBinsXF[25] =  0.225;
  fPiBinsXF[26] =  0.275;
  fPiBinsXF[27] =  0.350;
  fPiBinsXF[28] =  0.450;
  fPiBinsXF[29] =  0.55;  
  
  fNPiBinsPT = 17;
  fPiBinsPT = new double[fNPiBinsPT];
  fPiBinsPT[0] = 0.025;
  fPiBinsPT[1] = 0.075;
  fPiBinsPT[2] = 0.125;
  fPiBinsPT[3] = 0.175;
  fPiBinsPT[4] = 0.225;
  fPiBinsPT[5] = 0.275;
  fPiBinsPT[6] = 0.35;
  fPiBinsPT[7] = 0.45;
  fPiBinsPT[8] = 0.55;
  fPiBinsPT[9] = 0.65;
  fPiBinsPT[10] = 0.75;
  fPiBinsPT[11] = 0.85;
  fPiBinsPT[12] = 0.95; 
  fPiBinsPT[13] = 1.1;
  fPiBinsPT[14] = 1.3;
  fPiBinsPT[15] = 1.5;
  fPiBinsPT[16] = 1.7;
  
  double pibins[] = { -0.550, -0.450, -0.350, -0.275, -0.225, -0.175, -0.1375, -0.1125, -0.0875, 
                      -0.0625, -0.045, -0.035, -0.025, -0.015,  
		      -0.005,
		      0.005,  0.015,  0.025,  0.035,  0.045,
		      0.0625,  0.0875,  0.1125,  0.1375,  0.175, 0.225,  0.275,  0.350,  0.450, 0.55 };  
  int npibins = sizeof(pibins) / sizeof(double) - 1;
//  double pibins_pt = { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8 }
  double pibins_pt[] = { 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 
                         0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 
		         1.1, 1.3, 1.5, 1.7 };
  int npibins_pt = sizeof(pibins_pt)/sizeof(double) - 1; 

  fHistoSecPiMinus.push_back( new TH1D( "pimimus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TH1D( "piminus_dNdxF",    title.c_str(), fNPiBinsXF-1, fPiBinsXF ) );
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT",   title.c_str(), fNPiBinsXF-1, fPiBinsXF, 0., 10. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT2",  title.c_str(), npibins, pibins, 0., 10. ) );
  
  fHistoPTvsXFPiMinus = new TH2D( "piminis_pTvsxF", title.c_str(), npibins, pibins, npibins_pt, pibins_pt );
  
  for ( int nb=0; nb<fNPiBinsXF; ++nb )
  {
     double xF = (fPiBinsXF[nb]+fPiBinsXF[nb+1])/2.;
     std::ostringstream osCount1;
     osCount1 << xF;
     std::ostringstream osCount2;
     osCount2 << nb;
     G4String subtitle = ", xF = " + osCount1.str();
     hname = "pTpim" + osCount2.str();
     fHistoPTPiMinus.push_back( new TH1D( hname.c_str(), (title+subtitle).c_str(), npibins_pt, pibins_pt) );     
  }
  
  fHistoSecPiMinusTot.push_back( new TH1D( "pimimus_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinusTot.push_back( new TH1D( "piminus_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinusTot.push_back( new TProfile( "piminus_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiMinusTot.push_back( new TProfile( "piminus_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  for ( NuERange r=kR_0_2; r<=kR_20_50; r=NuERange(r+1) )
  {

     std::cout << " r(ange) = " << r << std::endl;
     
     std::ostringstream osCount;
     // osCount.clear();
     osCount << r;
     
     std::cout << " osCount = " << osCount.str() << std::endl;
     
     hname = pim_names[0] + osCount.str();
     fHistoPiMinusFrom1stInt[r].push_back( new TH1D( hname.c_str(), title.c_str(), npibins, pibins ) );
     hname = pim_names[1] + osCount.str(); 	
     fHistoPiMinusFrom1stInt[r].push_back( new TProfile( hname.c_str(),   title.c_str(), npibins, pibins, 0., 10. ) );	

     hname = "sec_proton_energy_to_pim_" + osCount.str();
     fHistoPiMinusFromReint[kFromP][r].push_back( new TH1D( hname.c_str(), title.c_str(), 300, 0., 150. ) );	

     hname = "sec_pim_energy_to_pim_" + osCount.str();
     fHistoPiMinusFromReint[kFromPim][r].push_back( new TH1D( hname.c_str(), title.c_str(), 300, 0., 150. ) );	

     hname = "sec_pip_energy_to_pim_" + osCount.str();
     fHistoPiMinusFromReint[kFromPip][r].push_back( new TH1D( hname.c_str(), title.c_str(), 300, 0., 150. ) );	

  }

  outcome = " -> X + pi+ ";
  title = htitle + outcome;
  
  // NOTE: use the same binning as for pi-

  fHistoSecPiPlus.push_back( new TH1D( "piplus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TH1D( "piplus_dNdxF",    title.c_str(), fNPiBinsXF-1, fPiBinsXF ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT",   title.c_str(), fNPiBinsXF-1, fPiBinsXF, 0., 10. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT2",  title.c_str(), npibins, pibins, 0., 10. ) );	

  fHistoPTvsXFPiPlus = new TH2D( "piplus_pTvsxF", title.c_str(), npibins, pibins, npibins_pt, pibins_pt );	

  for ( int nb=0; nb<fNPiBinsXF; ++nb )
  {
     double xF = (fPiBinsXF[nb]+fPiBinsXF[nb+1])/2.;
     std::ostringstream osCount1;
     osCount1 << xF;
     std::ostringstream osCount2;
     osCount2 << nb;
     G4String subtitle = ", xF = " + osCount1.str();
     hname = "pTpip" + osCount2.str();
     fHistoPTPiPlus.push_back( new TH1D( hname.c_str(), (title+subtitle).c_str(), npibins_pt, pibins_pt ) );     
  }

  fHistoSecPiPlusTot.push_back( new TH1D( "piplus_dXSecdxF_total", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlusTot.push_back( new TH1D( "piplus_dNdxF_total",    title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlusTot.push_back( new TProfile( "piplus_pT_total",   title.c_str(), 100, -1., 1., 0., 10. ) );	
  fHistoSecPiPlusTot.push_back( new TProfile( "piplus_pT2_total",  title.c_str(), 100, -1., 1., 0., 10. ) );	

  for ( NuERange r=kR_0_2; r<=kR_20_50; r=NuERange(r+1) )
  {

     std::ostringstream osCount;
     // osCount.clear();
     osCount << r;
     hname = pip_names[0] + osCount.str();
     fHistoPiPlusFrom1stInt[r].push_back( new TH1D( hname.c_str(), title.c_str(), npibins, pibins ) );
     hname = pip_names[1] + osCount.str(); 	
     fHistoPiPlusFrom1stInt[r].push_back( new TProfile( hname.c_str(),   title.c_str(), npibins, pibins, 0., 10. ) );	

     hname = "sec_proton_energy_to_pip_" + osCount.str();
     fHistoPiPlusFromReint[kFromP][r].push_back( new TH1D( hname.c_str(), title.c_str(), 300, 0., 150. ) );	
//     hname = "from_p_" + pip_names[0] + osCount.str();
//     fHistoPiPlusFromReint[kFromP][r].push_back( new TH1D( hname.c_str(), title.c_str(), 40, -1., 1. ) );	
//     hname = "from_p_" + pip_names[1] + osCount.str();
//     fHistoPiPlusFromReint[kFromP][r].push_back( new TProfile( hname.c_str(),   title.c_str(), 100, -1., 1., 0., 10. ) );	

     hname = "sec_pim_energy_to_pip_" + osCount.str();
     fHistoPiPlusFromReint[kFromPim][r].push_back( new TH1D( hname.c_str(), title.c_str(), 300, 0., 150. ) );	
//     hname = "from_pim_" + pip_names[0] + osCount.str();
//     fHistoPiPlusFromReint[kFromPip][r].push_back( new TH1D( hname.c_str(), title.c_str(), 40, -1., 1. ) );	
//     hname = "from_pim_" + pip_names[1] + osCount.str();
//     fHistoPiPlusFromReint[kFromPip][r].push_back( new TProfile( hname.c_str(),   title.c_str(), 100, -1., 1., 0., 10. ) );	

     hname = "sec_pip_energy_to_pip_" + osCount.str();
     fHistoPiPlusFromReint[kFromPip][r].push_back( new TH1D( hname.c_str(), title.c_str(), 300, 0., 150. ) );	
//     hname = "from_pip_" + pip_names[0] + osCount.str();
//     fHistoPiPlusFromReint[kFromPip][r].push_back( new TH1D( hname.c_str(), title.c_str(), 40, -1., 1. ) );	
//     hname = "from_pip_" + pip_names[1] + osCount.str();
//     fHistoPiPlusFromReint[kFromPip][r].push_back( new TProfile( hname.c_str(),   title.c_str(), 100, -1., 1., 0., 10. ) );	

  }

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
   
   // delete fHistoPTvsXFProton;
   // delete fHistoPTvsXFAntiProton;
   delete fHistoPTvsXFPiMinus;
   delete fHistoPTvsXFPiPlus;
   // delete  fHistoPTvsXFNeutron;
   
   for ( size_t i=0; i<fHistoPTPiMinus.size(); ++i )
   {
      delete fHistoPTPiMinus[i];
   }
   for ( size_t i=0; i<fHistoPTPiPlus.size(); ++i )
   {
      delete fHistoPTPiPlus[i];
   }
   
   delete [] fPiBinsXF;
   
   for ( NuERange r=kR_0_2; r<=kR_20_50; r=NuERange(r+1) )
   {

      int s1 = fHistoPiMinusFrom1stInt[r].size();
      for ( int i=0; i<s1; ++i )
      {
         delete fHistoPiMinusFrom1stInt[r][i];
	 fHistoPiMinusFrom1stInt[r][i] = 0;
         delete fHistoPiPlusFrom1stInt[r][i];
	 fHistoPiPlusFrom1stInt[r][i] = 0;
      }
      fHistoPiMinusFrom1stInt[r].clear();
      fHistoPiPlusFrom1stInt[r].clear();
      
      int s2 = fHistoPiMinusFromReint[kFromP][r].size();
      for ( int j=0; j<s2; ++j )
      {
         delete fHistoPiMinusFromReint[kFromP][r][j];
	 fHistoPiMinusFromReint[kFromP][r][j] = 0;
      }
      fHistoPiMinusFromReint[kFromP][r].clear();
      s2 = fHistoPiMinusFromReint[kFromPim][r].size();
      for ( int j=0; j<s2; ++j )
      {
         delete fHistoPiMinusFromReint[kFromPim][r][j];
	 fHistoPiMinusFromReint[kFromPim][r][j] = 0;
      }
      fHistoPiMinusFromReint[kFromPim][r].clear();
      s2 = fHistoPiMinusFromReint[kFromPip][r].size();
      for ( int j=0; j<s2; ++j )
      {
         delete fHistoPiMinusFromReint[kFromPip][r][j];
	 fHistoPiMinusFromReint[kFromPip][r][j] = 0;
      }
      fHistoPiMinusFromReint[kFromPip][r].clear();
      
      s2 = fHistoPiPlusFromReint[kFromP][r].size();
      for ( int j=0; j<s2; ++j )
      {
         delete fHistoPiPlusFromReint[kFromP][r][j];
	 fHistoPiPlusFromReint[kFromP][r][j] = 0;
      }
      fHistoPiPlusFromReint[kFromP][r].clear();
      s2 = fHistoPiPlusFromReint[kFromPim][r].size();
      for ( int j=0; j<s2; ++j )
      {
         delete fHistoPiPlusFromReint[kFromPim][r][j];
	 fHistoPiPlusFromReint[kFromPim][r][j] = 0;
      }
      fHistoPiPlusFromReint[kFromPim][r].clear();
      s2 = fHistoPiPlusFromReint[kFromPip][r].size();
      for ( int j=0; j<s2; ++j )
      {
         delete fHistoPiPlusFromReint[kFromPip][r][j];
	 fHistoPiPlusFromReint[kFromPip][r][j] = 0;
      }
      fHistoPiPlusFromReint[kFromPip][r].clear();

   }

}

void Tst23NA49Histo::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& labp ) 
{

   G4int NSec = aChange->GetNumberOfSecondaries();
 
   if ( NSec <= 0 ) return;

   fNuERange    = kNone;
   fNubarERange = kNone;
  
   bool isFirst = (dynamic_cast<Tst23ParticleChange*>(aChange))->IsFisrtInteraction();
   
   if ( isFirst && fDoResDecay ) AccountForResDecay( aChange );
      
   double SQRT_S = 17.2;
   
   G4ThreeVector boostLabp = labp.boostVector();
   
/*
   double BeamEnergy = labp.t() - 2.*G4Proton::Proton()->GetPDGMass()/GeV;
    
   // well, this way sqrt(s) will actually end up at 17.27 which is a bit more than 17.2 (taken from SPS hadrons site)
   // need to figure out and refine
   double Ecm   = sqrt(2.0*pow(0.938,2.0)+2.*BeamEnergy*0.938);
   double beta  = sqrt(pow(BeamEnergy,2.0)-pow(0.938,2.0))/(BeamEnergy+0.938);
   Ecm   = sqrt(2.0*pow(G4Proton::Proton()->GetPDGMass()/GeV,2.0)+2.*BeamEnergy*G4Proton::Proton()->GetPDGMass()/GeV);
   beta  = sqrt(pow(BeamEnergy,2.0)-pow(G4Proton::Proton()->GetPDGMass()/GeV,2.0))/(BeamEnergy+G4Proton::Proton()->GetPDGMass()/GeV);
   double gamma = 1./sqrt(1.-pow(beta,2.0));
    
   double beta1 = boostLabp.beta();
   double gamma1 = boostLabp.gamma();
*/   
   if ( isFirst && NSec > 0 ) fHistoNSec->Fill( (double)NSec );
     
   const G4Track*           trk = 0;  
   const G4DynamicParticle* sec = 0;
       
   const G4Track*           incomingtrk = (dynamic_cast<Tst23ParticleChange*>(aChange))->GetIncomingTrack();
   const G4DynamicParticle* incoming    = incomingtrk->GetDynamicParticle();
   InteractionType IntType = kOther;

   if ( incoming->GetDefinition()->GetParticleName() == "proton" )
   {
      IntType = kFromP;
   }
   else if ( incoming->GetDefinition()->GetParticleName() == "pi-" )
   {
      IntType = kFromPim;
   }
   else if ( incoming->GetDefinition()->GetParticleName() == "pi+" )
   {
      IntType = kFromPip;
   }
   
   for (G4int i=0; i<NSec; i++) 
   {
        trk = aChange->GetSecondary(i);
        sec = trk->GetDynamicParticle();			
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
	      int nb = -1;
	      for ( int ib=0; ib<fNPiBinsXF; ++ib )
	      {
	         if ( xF >= fPiBinsXF[ib] && xF < fPiBinsXF[ib+1] )
		 {
		    nb = ib;
		    break;
		 }
	      }
	      if ( nb == -1 ) continue;
	      // calculate weight as Epart/dP3 (see NA49 papers on the SPS hadrons site at CERN) 
	      double wei = CalculateBinWeight( labp, ekin, mass, pT, nb, SQRT_S );
	      fHistoPTPiMinus[nb]->Fill( pT, wei ); 
	      fHistoPTvsXFPiMinus->Fill( xF, pT, wei ); // 2D plot
	   }
	   fHistoSecPiMinusTot[1]->Fill( xF );
	   fHistoSecPiMinusTot[2]->Fill( xF, pT );
	   fHistoSecPiMinusTot[3]->Fill( xF, pT2 );
	   //
	   AccountForPionDecay( trk );
	   if ( fNubarERange == kNone ) continue;	   
	   if ( isFirst )
	   {
	      fHistoPiMinusFrom1stInt[fNubarERange][0]->Fill( xF );
	      fHistoPiMinusFrom1stInt[fNubarERange][1]->Fill( xF, pT );
	   }
	   else // reinteraction 
	   {
	      if ( IntType == kOther ) continue;
	      fHistoPiMinusFromReint[IntType][fNubarERange][0]->Fill( incomingtrk->GetVertexKineticEnergy()/GeV );
           }
	}
	else if ( pname == "pi+" )
	{
	   if ( isFirst )
	   {
	      fHistoSecPiPlus[1]->Fill( xF );
	      fHistoSecPiPlus[2]->Fill( xF, pT );
	      fHistoSecPiPlus[3]->Fill( xF, pT2 );
	      int nb = -1;
	      for ( int ib=0; ib<fNPiBinsXF; ++ib )
	      {
	         if ( xF > fPiBinsXF[ib] && xF <= fPiBinsXF[ib+1] )
		 {
		    nb = ib;
		    break;
		 }
	      }
	      if ( nb == -1 ) continue;
	      // calculate weight as Epart/dP3 (see NA49 papers on the SPS hadrons site at CERN) 
	      double wei = CalculateBinWeight( labp, ekin, mass, pT, nb, SQRT_S );
	      fHistoPTPiPlus[nb]->Fill( pT, wei ); 
	      fHistoPTvsXFPiPlus->Fill( xF, pT, wei );
	   }
	   fHistoSecPiPlusTot[1]->Fill( xF );
	   fHistoSecPiPlusTot[2]->Fill( xF, pT );
	   fHistoSecPiPlusTot[3]->Fill( xF, pT2 );
	   //
	   AccountForPionDecay( trk );
	   if ( fNuERange == kNone ) continue;
	   if ( isFirst )
	   {
	      fHistoPiPlusFrom1stInt[fNuERange][0]->Fill( xF );
	      fHistoPiPlusFrom1stInt[fNuERange][1]->Fill( xF, pT );
	   }
	   else // reiteraction
	   {
	      if ( IntType == kOther ) continue;
	      fHistoPiPlusFromReint[IntType][fNuERange][0]->Fill( incomingtrk->GetVertexKineticEnergy()/GeV );
	   }
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

   double norm = 1.;
   
   norm = fHistoNSec->Integral();
   
   fHistoNSec->Scale( 1./norm, "width" );
   fHistoNSec->Write();
   
   // secondary proton
   //
   for ( size_t i=0; i<fHistoSecProton.size(); ++i )   
   {
      if ( i <2 ) fHistoSecProton[i]->Scale(1./norm, "width");
      fHistoSecProton[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecProtonTot.size(); ++i )   
   {
      if ( i <2 ) fHistoSecProtonTot[i]->Scale(1./norm, "width");
      fHistoSecProtonTot[i]->Write();
   }
   // secondary pbar
   //
   for ( size_t i=0; i<fHistoSecAntiProton.size(); ++i )   
   {
      if ( i<2 ) fHistoSecAntiProton[i]->Scale(1./norm,"width");      
      fHistoSecAntiProton[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecAntiProtonTot.size(); ++i )   
   {
      if ( i <2 ) fHistoSecAntiProtonTot[i]->Scale(1./norm,"width");
      fHistoSecAntiProtonTot[i]->Write();
   }

   // secondary pi-
   //
   for ( size_t i=0; i<fHistoSecPiMinus.size(); ++i )
   {
      if ( i < 2 ) fHistoSecPiMinus[i]->Scale(1./norm,"width");
      fHistoSecPiMinus[i]->Write();
   }

   fHistoPTvsXFPiMinus->Scale( 1./norm, "width" );
   fHistoPTvsXFPiMinus->Write();

   for ( size_t i=0; i<fHistoPTPiMinus.size(); ++i )
   {
      fHistoPTPiMinus[i]->Scale( 1./norm ); // Note: NO scaling with "width" because it had to be E/dP3
                                            //       which is already taken into account as a weight 
      fHistoPTPiMinus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiMinusTot.size(); ++i )
   {
      if ( i < 2 ) fHistoSecPiMinusTot[i]->Scale(1./norm,"width");
      fHistoSecPiMinusTot[i]->Write();
   }

   // secondary pi+
   //
   for ( size_t i=0; i<fHistoSecPiPlus.size(); ++i )
   {
      if ( i < 2 )
      {
/*
         int nbins = fHistoSecPiPlus[i]->GetXaxis()->GetNbins();
         for ( int ib = 1; ib <= nbins; ib++ )
         { 
            xbin = fHistoSecPiPlus[i]->GetBinWidth(ib);
	    scale = 1. / (norm*xbin) ;
	    double yvalue = fHistoSecPiPlus[i]->GetBinContent(ib);
	    yvalue *= scale;
	    // double yerror = fHistoSecPiPlus[i]->GetBinError(ib);
	    // yerror *= scale;
	    fHistoSecPiPlus[i]->SetBinContent( ib, yvalue );
	    // fHistoSecPiPlus[i]->SetBinError( ib, yerror );
         }
*/
         fHistoSecPiPlus[i]->Scale( 1./norm, "width" );
      }
      fHistoSecPiPlus[i]->Write();
   }
   for ( size_t i=0; i<fHistoPTPiPlus.size(); ++i )
   {
      fHistoPTPiPlus[i]->Scale( 1./norm ); // NO scaling with "width" - see earlier comment for fHistoPTPiMinus
      fHistoPTPiPlus[i]->Write();
   }
//   fHistoPTvsXFPiPlus->Scale( 1./norm, "width" );
   int nbx = fHistoPTvsXFPiPlus->GetXaxis()->GetNbins();
   int nby = fHistoPTvsXFPiPlus->GetYaxis()->GetNbins();
   for ( int ibx=1; ibx<=nbx; ++ibx )
   {
      double xb = fHistoPTvsXFPiPlus->GetXaxis()->GetBinWidth(ibx);
      for ( int iby=1; iby<=nby; ++iby )
      {
         double yb = fHistoPTvsXFPiPlus->GetYaxis()->GetBinWidth(iby);
	 double content = fHistoPTvsXFPiPlus->GetBinContent(ibx,iby);
	 double error   = fHistoPTvsXFPiPlus->GetBinError(ibx,iby);
	 double scale = 1./(norm*xb*yb);
	 fHistoPTvsXFPiPlus->SetBinContent(ibx,iby,content*scale);
	 fHistoPTvsXFPiPlus->SetBinError(ibx,iby,error*scale);
      }
   }
   fHistoPTvsXFPiPlus->Write();
   for ( size_t i=0; i<fHistoSecPiPlusTot.size(); ++i )
   {
      if ( i < 2 ) fHistoSecPiPlusTot[i]->Scale(1./norm,"width");
      fHistoSecPiPlusTot[i]->Write();
   }

   // secondary neutron
   //
   for ( size_t i=0; i<fHistoSecNeutron.size(); ++i )
   {
      if ( i < 2 ) fHistoSecNeutron[i]->Scale(1./norm,"width");
      fHistoSecNeutron[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecNeutronTot.size(); ++i )
   {
      if ( i < 2 ) fHistoSecNeutronTot[i]->Scale(1./norm,"width");
      fHistoSecNeutronTot[i]->Write();
   }
   
   for ( NuERange r=kR_0_2; r<=kR_20_50; r=NuERange(r+1) )
   {
      fHistoPiMinusFrom1stInt[r][0]->Scale(1./norm,"width");
      fHistoPiMinusFrom1stInt[r][0]->Write();
      // no normalization for pT vs xF (profile histo)
      fHistoPiMinusFrom1stInt[r][1]->Write();

      fHistoPiPlusFrom1stInt[r][0]->Scale(1./norm,"width");
      fHistoPiPlusFrom1stInt[r][0]->Write();
      // no normalization for pT vs xF (profile histo)
      fHistoPiPlusFrom1stInt[r][1]->Write();
      
// NO normalization (yet)...
//
      fHistoPiPlusFromReint[kFromP][r][0]->Write();
      fHistoPiPlusFromReint[kFromPim][r][0]->Write();
      fHistoPiPlusFromReint[kFromPip][r][0]->Write();

   }

   return;

}

double Tst23NA49Histo::CalculateBinWeight( const G4LorentzVector& labp, double Ekin, double mass, double pT, int xFbin, double sqrt_s )
{

   double wei = 1.;
   int pTbin = -1;
   for (  int ib=0; ib<fNPiBinsPT; ++ib )
   {
      if ( pT >= fPiBinsPT[ib] && pT < fPiBinsPT[ib+1] )
      {
         pTbin = ib;
	 break;
      }
   }
   if ( pTbin == -1 ) return wei;

// NOTE: This fragment of code draws inspiration in a similar application
//       originally implemented by Mike Kordosky (W&M/MINERvA/NuMI).
//       Credits go to Mike !!!
//
   double dpT2 = fPiBinsPT[pTbin+1]*fPiBinsPT[pTbin+1] - fPiBinsPT[pTbin]*fPiBinsPT[pTbin];
   double pLmin = fPiBinsXF[xFbin]*sqrt_s/2.;
   double pLmax = fPiBinsXF[xFbin+1]*sqrt_s/2.;
   double EPCM1 = sqrt( (pLmin*pLmin+pT*pT) + mass*mass );
   double EPCM2 = sqrt( (pLmax*pLmax+pT*pT) + mass*mass );
   double pZmin = labp.boostVector().gamma()*(labp.beta()*EPCM1 + pLmin );
   double pZmax = labp.boostVector().gamma()*(labp.beta()*EPCM2 + pLmax );
   double dP3   = CLHEP::pi * dpT2 *(pZmax-pZmin);
   wei = (Ekin+mass)/dP3;

   return wei;
   
}

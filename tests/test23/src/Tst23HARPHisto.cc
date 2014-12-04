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

#include "G4PhysicalConstants.hh"


#include "Tst23HARPHisto.hh"
#include "Tst23ParticleChange.hh"

#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include <sstream>

Tst23HARPHisto::Tst23HARPHisto( G4String ht )
   : TstHistoSet( ht ), 
     fNThetaBinsFW(4), fThetaMinFW(0.05), fDeltaThetaFW(0.05),
     fNThetaBinsLA(9), fThetaMinLA(0.35), fDeltaThetaLA(0.2)
{

   G4String htitle;
   G4String hname;
  
   fHistoNSec = new TH1D( "NSec", ht.c_str(), 100, 0., 100. );

//   std::ostringstream osTitle1(std::ios_base::out|std::ios_base::app);
//   std::ostringstream osTitle2(std::ios_base::out|std::ios_base::app);
//   std::ostringstream osTitle3(std::ios_base::out|std::ios_base::app);
       
   G4double thetaMin = 0.;
   G4double thetaMax = 0.;
   std::string theta_bin_fw;
   std::string theta_bin_la;

   for ( G4int i=0; i<fNThetaBinsFW; i++ )
   {
      thetaMin = fThetaMinFW + fDeltaThetaFW*i;
      thetaMax = thetaMin + fDeltaThetaFW;
      
      std::ostringstream osTitle1;
      std::ostringstream osTitle2;
      std::ostringstream osTitle3;
      
      osTitle1.clear();
      osTitle1 << thetaMin;
      theta_bin_fw = osTitle1.str() + " < theta < ";
      osTitle2.clear();
      osTitle2 << thetaMax;
      theta_bin_fw += osTitle2.str();
      theta_bin_fw += "(rad)";
   
      osTitle3.clear();
      osTitle3 << i;
      
      htitle = ht + " -> X + pi-, " + theta_bin_fw;  
      hname = "piminus_FW_" + osTitle3.str();         
      fHistoSecPiMinusFW.push_back( new TH1D( hname.c_str(), htitle.c_str(), 80, 0., 8.0 ) );

      htitle = ht + " -> X + pi+, " + theta_bin_fw;
      hname = "piplus_FW_" + osTitle3.str();
      fHistoSecPiPlusFW.push_back( new TH1D( hname.c_str(), htitle.c_str(), 80, 0., 8.0 ) );
   }
   
   for ( G4int i=0; i<fNThetaBinsLA; i++ )
   {
     
      thetaMin = fThetaMinLA + fDeltaThetaLA*i;
      thetaMax = thetaMin + fDeltaThetaLA; 
     
      std::ostringstream osTitle1;
      std::ostringstream osTitle2;
      std::ostringstream osTitle3;

      osTitle1.clear();
      osTitle1 << thetaMin;
      theta_bin_la = osTitle1.str() + " < theta < ";
      osTitle2.clear();
      osTitle2 << thetaMax;
      theta_bin_la += osTitle2.str();
      theta_bin_la += "(rad)";
      
      osTitle3.clear();
      osTitle3 << i;
      
      htitle = ht + " -> X + pi-, " + theta_bin_la;  
      hname = "piminus_LA_" + osTitle3.str();         
      fHistoSecPiMinusLA.push_back( new TH1D( hname.c_str(), htitle.c_str(), 50, 0., 2.5 ) );
//      hname = "pimimus_mom_total_" + osTitle3.str();
//      fHistoSecPiMinusTot.push_back( new TH1D( hname.c_str(), htitle.c_str(), 50, 0., 2.5) );	

      htitle = ht + " -> X + pi+, " + theta_bin_la;
      hname = "piplus_LA_" + osTitle3.str();
      fHistoSecPiPlusLA.push_back( new TH1D( hname.c_str(), htitle.c_str(), 50, 0., 2.5 ) );
//      hname = "piplus_mom_total_" + osTitle3.str();
//      fHistoSecPiPlusTot.push_back( new TH1D( hname.c_str(), htitle.c_str(), 50, 0., 2.5 ) );	

   }

}

Tst23HARPHisto::~Tst23HARPHisto()
{

   if ( fHistoNSec ) delete fHistoNSec;
   
   for (size_t i=0; i<fHistoSecPiMinusLA.size(); i++ )
   {
      delete fHistoSecPiMinusLA[i];
   }
   for (size_t i=0; i<fHistoSecPiPlusLA.size(); i++ )
   {
      delete fHistoSecPiPlusLA[i];
   }

/*
   for (size_t i=0; i<fHistoSecPiMinusTot.size(); i++ )
   {
      delete fHistoSecPiMinusTot[i];
   }
   for (size_t i=0; i<fHistoSecPiPlusTot.size(); i++ )
   {
      delete fHistoSecPiPlusTot[i];
   }
*/

}

void Tst23HARPHisto::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector&  ) 
{

   
   G4bool isFirst = (dynamic_cast<Tst23ParticleChange*>(aChange))->IsFisrtInteraction();
      
   G4int NSec = aChange->GetNumberOfSecondaries();

   if ( isFirst && NSec > 0 ) fHistoNSec->Fill( (double)NSec );
     
   const G4DynamicParticle* sec = 0;
          
   for (G4int i=0; i<NSec; i++) 
   {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();			
	const G4String& pname = sec->GetDefinition()->GetParticleName();
		
	// G4ThreeVector mom = sec->GetMomentum() / GeV ;
	G4double pmom = sec->GetTotalMomentum() / GeV ;
	//double mass = sec->GetDefinition()->GetPDGMass() / GeV;
	//double ekin = sec->GetKineticEnergy() / GeV ;
	G4double theta = (sec->GetMomentum()).theta();
	
	if ( theta < fThetaMinFW ) continue;
	if ( theta < fThetaMinFW+fDeltaThetaFW*fNThetaBinsFW )
	{
	   G4int ith = ( theta - fThetaMinFW ) / fDeltaThetaFW;
	   if ( pname == "pi-" )
	   {
	      if ( isFirst ) fHistoSecPiMinusFW[ith]->Fill( pmom );
	   }
	   else if ( pname == "pi+" )
	   {
	      if ( isFirst ) fHistoSecPiPlusFW[ith]->Fill( pmom );
	   }
	}
		
	if ( theta < fThetaMinLA ) continue;
	if ( theta > fThetaMinLA+fDeltaThetaLA*fNThetaBinsLA ) continue;
	G4int    itheta = ( theta - fThetaMinLA ) / fDeltaThetaLA;
	if ( itheta < 0 || itheta >= fNThetaBinsLA ) continue;
	        
	if ( pname == "pi-" )
	{
	   if ( isFirst )
	   {
	      fHistoSecPiMinusLA[itheta]->Fill( pmom );
	   }
	   // fHistoSecPiMinusTot[itheta]->Fill( pmom );
	}
	else if ( pname == "pi+" )
	{
	   if ( isFirst )
	   {
	      fHistoSecPiPlusLA[itheta]->Fill( pmom );
	   }
	   // fHistoSecPiPlusTot[itheta]->Fill( pmom );
	}
   }
      
   return;
   
}

void Tst23HARPHisto::Write( G4int, G4double xsec )
{

   // fHistoFile->cd();

   G4double xbin = 1.;
   G4double norm = 1.;
   G4double scale = 1.;
   
   norm = fHistoNSec->Integral();
   
   xbin = (G4double)fHistoNSec->GetBinWidth(1);
   scale = 1. / (xbin*norm);
   fHistoNSec->Scale( scale ) ;
   fHistoNSec->Write();
   
   // secondary pi-
   //
   for ( size_t i=0; i<fHistoSecPiMinusFW.size(); ++i )
   {
      xbin = fHistoSecPiMinusFW[i]->GetBinWidth(1);
      scale = xsec / (xbin*norm*fDeltaThetaFW);
      fHistoSecPiMinusFW[i]->Scale(scale);
      fHistoSecPiMinusFW[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecPiMinusLA.size(); ++i )
   {
      xbin = fHistoSecPiMinusLA[i]->GetBinWidth(1);
      scale = xsec / (xbin*norm*fDeltaThetaLA);
      fHistoSecPiMinusLA[i]->Scale(scale);
      fHistoSecPiMinusLA[i]->Write();
   }
/*
   for ( size_t i=0; i<fHistoSecPiMinusTot.size(); ++i )
   {
      xbin = fHistoSecPiMinusTot[i]->GetBinWidth(1);
      scale = xsec / (xbin*norm*fDeltaTheta);
      fHistoSecPiMinusTot[i]->Scale(scale);
      fHistoSecPiMinusTot[i]->Write();
   }
*/
   // secondary pi+
   //
   for ( size_t i=0; i<fHistoSecPiPlusFW.size(); ++i )
   {
      xbin = fHistoSecPiPlusFW[i]->GetBinWidth(1);
      scale = xsec / (xbin*norm*fDeltaThetaFW);
      fHistoSecPiPlusFW[i]->Scale(scale);
      fHistoSecPiPlusFW[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecPiPlusLA.size(); ++i )
   {
      xbin = fHistoSecPiPlusLA[i]->GetBinWidth(1);
      scale = xsec / (xbin*norm*fDeltaThetaLA);
      fHistoSecPiPlusLA[i]->Scale(scale);
      fHistoSecPiPlusLA[i]->Write();
   }
/*
   for ( size_t i=0; i<fHistoSecPiPlusTot.size(); ++i )
   {
      xbin = fHistoSecPiPlusTot[i]->GetBinWidth(1);
      scale = xsec / (xbin*norm*fDeltaTheta);
      fHistoSecPiPlusTot[i]->Scale(scale);
      fHistoSecPiPlusTot[i]->Write();
   }
*/
   return;

}

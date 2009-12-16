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
// $Id: G4CrossSectionTest.cc,v 1.4 2009-12-16 17:50:31 gunter Exp $
// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4CrossSectionTest.cc
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it), 
// 
//      Creation date: 15 April 1999
//
//      Modifications:
//
//  Dec- 2008 GF : migrate to root 
//
// -------------------------------------------------------------------

#include "globals.hh"

#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <assert.h>

#include "Randomize.hh"

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

#include "G4AntiProton.hh"
#include "G4Proton.hh"
#include "G4AntiNeutron.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4Lambda.hh"

#include "G4ParticleDefinition.hh"

#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"

#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4VCrossSectionSource.hh"
#include "G4XAqmTotal.hh"
#include "G4XAqmElastic.hh"
#include "G4XNNTotalLowE.hh"
#include "G4XPDGTotal.hh"
#include "G4XNNTotal.hh"
#include "G4XNNElastic.hh"
#include "G4XPDGElastic.hh"
#include "G4XNNElasticLowE.hh"
#include "G4XnpElastic.hh"
#include "G4XResonance.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TTree.h"

#define G4FPE_DEBUG=1

#ifdef G4FPE_DEBUG
  #include "G4FPEDetection.hh"
#endif

G4ParticleDefinition* SelectTrackDef(int i);

int main()
{

#ifdef G4FPE_DEBUG
  InvalidOperationDetection();
#endif
    std::string hFile="xsec.root"; 
    TFile * rootFile = new TFile(hFile.c_str(),"CREATE");
    if ( ! rootFile ) 
    {
       std::cout << " Fail to create root file " << std::endl;
       exit(1);
    }   

        TTree * Txs = new TTree("SEC","Secondary info");

        struct xs_info {
             float cmsE;
             float sigma;
        };
        xs_info xsec;     
        Txs->Branch("secondaries",&xsec.cmsE,"cmsE/F:sigma/F");

 
  TH1 * hSigma = new TH1F("1","CrossSection", 100,0.,10.);


/*  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();

  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiProton = G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();   
  G4ParticleDefinition* antiNeutron = G4AntiNeutron::AntiNeutronDefinition();   

  G4ParticleDefinition* pionPlus = G4PionPlus::PionPlusDefinition();
  G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();
  G4ParticleDefinition* pionZero = G4PionZero::PionZeroDefinition();

  G4ParticleDefinition* kaonPlus = G4KaonPlus::KaonPlusDefinition();
  G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();
   
  G4ParticleDefinition* lambda = G4Lambda::LambdaDefinition();

  G4ParticleDefinition* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4ParticleDefinition* theMuonMinus = G4MuonMinus::MuonMinusDefinition();

  G4ParticleDefinition* theNeutrinoMu = G4NeutrinoMu::NeutrinoMuDefinition();
  G4ParticleDefinition* theAntiNeutrinoMu = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
*/
  G4BaryonConstructor Baryons;
  Baryons.ConstructParticle();
  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();
  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();

  G4ParticleDefinition* definition1= SelectTrackDef(1);
  G4double mass1 = definition1->GetPDGMass();

  // Formation time
  G4double time1 = 0.;
  //  G4cout << "Formation time (ns)" << G4endl;
  //  G4cin >> time1;
  //  if (time1 < 0.) time1 = 0.;
  //   time1 = time1 * ns;

  
  // Initial position
  G4double x1 = 0.;
  G4double y1 = 0.;
  G4double z1 = 0.;
  //  G4cout <<"Initial position (cm) " << G4endl;
  //  G4cin >> x1 >> y1 >> z1;
  G4ThreeVector position1(x1*cm, y1*cm, z1*cm);

  // Momentum
  G4double pzMin1, pzMax1;
  G4cout << "Pz min and max (GeV)" << G4endl;
  G4cin >> pzMin1 >> pzMax1;
  pzMin1 = pzMin1 * GeV;
  pzMax1 = pzMax1 * GeV;

  // Select track 2
 
  G4cout << "---- Set KineticTrack 2 ----" << G4endl;
  G4int id2 = 0;
  G4cout << "1 p   2 ap   3 n   4 an   5 pi+   6 pi-   7 K+   8 K- " << G4endl;
  G4cin >> id2;

  G4ParticleDefinition* definition2= SelectTrackDef(2);
  G4double mass2 = definition2->GetPDGMass();

  // Formation time
  G4double time2 = 0.;
  //  G4cout << "Formation time (ns)" << G4endl;
  //  G4cin >> time2;
  //  if (time2 < 0.) time2 = 0.;
  //  time2 = time2 * ns;
  
  // Initial position
  G4double x2 = 0.;
  G4double y2 = 0.;
  G4double z2 = 0.;
  //  G4cout <<"Initial position (cm) " << G4endl;
  //  G4cin >> x2 >> y2 >> z2;
  G4ThreeVector position2(x2*cm, y2*cm, z2*cm);

  // Momentum
  G4double pzMin2;
  G4double pzMax2;
  G4cout << "Pz min and max (GeV)" << G4endl;
  G4cin >> pzMin2 >> pzMax2;
  pzMin2 = pzMin2 * GeV;
  pzMax2 = pzMax2 * GeV;

  // Select CrossSectionSource

  G4cout << "CrossSectionSource: " << G4endl
	 <<  "1) AqmElastic " 
	 <<  "2) AqmTotal "
	 <<  "3) NNElastic "
	 <<  "4) NNElasticLowE "
	 <<  "5) NNTotal "         << G4endl
	 <<  "6) NNTotalLowE " 
	 <<  "7) PDGElastic " 
	 <<  "8) PDGTotal "
	 <<  "9) Resonance "
	 <<  "10) XnpElastic "
	 << G4endl;

  G4int ids;
  G4cin >> ids;

  G4VCrossSectionSource* source;

  switch (ids)
    {
    case 1:
      source = new G4XAqmElastic;
      break;
    case 2:
      source = new G4XAqmTotal;
      break;
    case 3:
      source = new G4XNNElastic;
      break;
    case 4:
      source = new G4XNNElasticLowE;
      break;
    case 5:
      source = new G4XNNTotal;
      break;
    case 6:
      source = new G4XNNTotalLowE;
      break;
    case 7:
      source = new G4XPDGElastic();
      break;
    case 8:
      source = new G4XPDGTotal();
      break;
    case 9:
      //      source = new G4XResonance;
      break;
    case 10:
            source = new G4XnpElastic;
      break;
    default:
      source = new G4XAqmTotal;
      break;
    }

  G4cout << "***** " << source->Name() <<  "***** created " << G4endl;
  source->Print();

  // Number of iterations
  G4int n = 0;
  G4cout << "---- Number of iterations ----" << G4endl;
  G4cin >> n;
  if (n < 1) n = 1;

  G4double step1 = (pzMax1 - pzMin1) / n;
  G4double step2 = (pzMax2 - pzMin2) / n;

  G4cout << "pz min max step " 
	 << pzMin1 / GeV << ", " 
	 << pzMax1 / GeV << ", " 
	 << step1  / GeV << ", " << G4endl 
	 << "pz min max step " 
	 << pzMin2 / GeV << ", " 
	 << pzMax2 / GeV << ", " 
	 << step2  / GeV << G4endl; 

  G4int i;
  for (i=0; i<n; i++)
    {
      G4double pStep1 = pzMin1 + step1 * i;
      G4ThreeVector p1(0.,0.,pStep1);
      G4double e1 = std::sqrt(mass1*mass1 + p1.mag()*p1.mag());
      G4LorentzVector p41(p1,e1);
      G4KineticTrack trk1(definition1,time1,position1,p41);

      G4double pStep2 = -(pzMin2 + step2 * i);
      G4ThreeVector p2(0.,0.,pStep2);
      G4double e2 = std::sqrt(mass2*mass2 + p2.mag()*p2.mag());
      G4LorentzVector p42(p2,e2);
      G4KineticTrack trk2(definition2,time2,position2,p42);

      G4double sqrts = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
      
      //      G4cout << "*-*-*-*-* Step " << i 
      //	     << " - pz1, pz2 = " << pStep1 << ", " << pStep2 
      //	     << " - e1, e2 = " << e1 << ", " << e2 
      //	     << " - p41.m(), p42.m() = " << trk1.Get4Momentum().m() << ", " << trk2.Get4Momentum().m()
      //	     << " ----- sqrtS = " << sqrts << G4endl;


      source->PrintAll(trk1,trk2);
     
      G4double sigma = source->CrossSection(trk1,trk2);
      
      //      G4cout << "Energy: " << sqrts / GeV 
      //	     << " ---- CrossSection = " << sigma / millibarn 
      //	     << G4endl;

      xsec.cmsE=sqrts/GeV;
      xsec.sigma=sigma/millibarn;
      Txs->Fill();
    }

  std::cout << "Committing..." << std::endl;
  rootFile->Write();


  delete source;

  return EXIT_SUCCESS;
}

G4ParticleDefinition* SelectTrackDef(int i)
{  
  G4cout << "---- Set KineticTrack " << i <<  " ----" << G4endl;
  G4int id = 0;
  G4cout << "1 p   2 ap   3 n   4 an   5 pi+   6 pi-   7 K+   8 K- " << G4endl;
  G4cin >> id;

  G4ParticleDefinition* pDef;
  switch (id)
    {
    case 1:
      pDef = G4Proton::Proton();
      break;
    case 2:
      pDef = G4AntiProton::AntiProton();
      break;
    case 3:
      pDef = G4Neutron::Neutron();
      break;
    case 4:
      pDef = G4AntiNeutron::AntiNeutron();
      break;
    case 5:
      pDef = G4PionPlus::PionPlus();
      break;
    case 6:
      pDef = G4PionMinus::PionMinus();
      break;
    case 7:
      pDef = G4KaonPlus::KaonPlus();
      break;
    case 8:
      pDef = G4KaonMinus::KaonMinus();
      break;
    default:
      pDef = G4Proton::Proton();
      break;
    }
   return pDef;
}

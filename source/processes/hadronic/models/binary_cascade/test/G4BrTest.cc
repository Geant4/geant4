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
// $Id: G4BrTest.cc,v 1.3 2009-12-16 17:50:25 gunter Exp $
// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4DecayTest.cc
//
//      Author:        Lorenzo Bellagamba (bellagamba@bo.infn.it), 
// 
//      Creation date: 30 October 1999
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"

#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <assert.h>

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "Randomize.hh"

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

#include "G4AntiProton.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4Lambda.hh"
#include "G4Eta.hh"
#include "G4LambdacPlus.hh"


#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"

#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"

#include "G4ParticleTable.hh"
#include "G4VCrossSectionSource.hh"
#include "G4XAqmTotal.hh"
#include "G4XAqmElastic.hh"
#include "G4XNNTotalLowE.hh"
#include "G4XPDGTotal.hh"
#include "G4XNNTotal.hh"
#include "G4XPDGElastic.hh"


int main()
{
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  hbookManager = new HBookFile("br.hbook", 58);
  assert (hbookManager != 0);

  // MGP ---- Book histograms
 
  // MGP ---- Book a ntuple
  HepTuple* ntuple;
  ntuple = hbookManager->ntuple("Decay ntuple");
  assert (ntuple != 0); 

  G4ShortLivedConstructor ShortLived;
  G4ShortLivedTable ShortLivedTab;
  ShortLived.ConstructParticle();

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();

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
  G4ParticleDefinition* eta = G4Eta::EtaDefinition();

  G4ParticleDefinition* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4ParticleDefinition* theMuonMinus = G4MuonMinus::MuonMinusDefinition();

  G4ParticleDefinition* theNeutrinoMu = G4NeutrinoMu::NeutrinoMuDefinition();
  G4ParticleDefinition* theAntiNeutrinoMu = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  //  G4ParticleDefinition* lambdacPlus = G4LambdacPlus::LambdacPlusDefinition();

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* N1440_Plus = particleTable->FindParticle("N(1440)+");
  G4ParticleDefinition* N1440_0 = particleTable->FindParticle("N(1440)0");
  G4ParticleDefinition* N1535_0 = particleTable->FindParticle("N(1535)0");
  // G4ParticleDefinition* deltaPluslus = particleTable->FindParticle("delta++");
  G4ParticleDefinition* deltaPlus = particleTable->FindParticle("delta+");
  G4ParticleDefinition* deltaMinus = particleTable->FindParticle("delta-");


  G4int DEBUG;
  cout << "insert Debug value 0/1" << G4endl;
  G4cin >> DEBUG;

  // Select track 1
  
  G4cout << "---- Set KineticTrack 1 ----" << G4endl;
  G4int id1 = 0;
  G4cout << "1 p   2 ap   3 n   4 an   5 pi+   6 pi-   7 K+   8 K-  9 N(1440)0  10 N(1440)+  11 delta+  12 N(1535)0   13 delta++" << G4endl;
  G4cin >> id1;

  G4ParticleDefinition* definition1;
  switch (id1)
    {
    case 1:
      definition1 = proton;
      break;
    case 2:
      definition1 = antiProton;
      break;
    case 3:
      definition1 = neutron;
      break;
    case 4:
      definition1 = antiNeutron;
      break;
    case 5:
      definition1 = pionPlus;
      break;
    case 6:
      definition1 = pionMinus;
      break;
    case 7:
      definition1 = kaonPlus;
      break;
    case 8:
      definition1 = kaonMinus;
      break;
    case 9:
      definition1 = N1440_0;
      break;
    case 10:
      definition1 = N1440_Plus;
      break;
    case 11:
      definition1 = deltaPlus;
      break;
    case 12:
      definition1 = N1535_0;
      break;
      //    case 13:
      //      definition1 = deltaPlusPlus;
      //      break;
    default:
      definition1 = proton;
      break;
    }
  if(DEBUG==1)cout << definition1 << G4endl;
  G4double mass0 = definition1->GetPDGMass();
  G4double width0 = definition1->GetPDGWidth();
  G4String name0;
  if(DEBUG==1)name0 = definition1->GetParticleName();
  if(DEBUG==1)cout << "Particle Type: " << name0 << G4endl;
  if(DEBUG==1)cout << "Particle Mass: " << mass0 << G4endl;
  if(DEBUG==1)cout << "Particle Width: " << width0 << G4endl;

  G4KineticTrack wig;
  cout << "Enter number of width to be consdered for BW mass distribtion" << G4endl;
  cout << " If Nwidth = 0  Mass = PDGMass " << G4endl;
  double nwid;
  G4cin >> nwid;

 // Formation time
  G4double time1 = 0.;

  // Initial position
  G4double x1 = 0.;
  G4double y1 = 0.;
  G4double z1 = 0.;
  G4ThreeVector pos1(x1*cm, y1*cm, z1*cm);

  // Momentum
  G4double px1 = 0. * GeV;
  G4double py1 = 0. * GeV;
  G4double pz1 = 0. * GeV; 
  cout << " enter Pz (Px=Py=0.)" << G4endl;
  G4cin >> pz1;

  G4int iter=0, iter_max;
  G4int nchan;
  cout << " Enter number of iteration" << G4endl;
  G4cin >> iter_max;
  G4double Dch[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  G4double max=wig.BrWig(width0, mass0, mass0);
  G4double Dmas=(2*nwid*width0)/iter_max;
  while(iter < iter_max){
    cout << " =================>   Processing event:......" << (iter+1) << G4endl;
  G4double mass1=(mass0-nwid*width0)+iter*Dmas;
  cout << mass1 << G4endl;
  G4double width1=width0;

  ntuple->column("mass",mass1);

  G4ThreeVector p1(px1,py1,pz1);
  G4double e1 = std::sqrt(mass1*mass1 + p1.mag()*p1.mag());
  G4LorentzVector p41(p1,e1);

  G4KineticTrack mytrack(definition1,time1,pos1,p41);

  ntuple->column("bwmas",mytrack.BrWig(width0,mass0,mass1));

  if(DEBUG==1)cout << mytrack.GetDefinition()->GetParticleName() << G4endl;
  G4DecayTable* TheDec=definition1->GetDecayTable();
  nchan=TheDec->entries(); 
  G4double* br=mytrack.GetActualWidth();
  G4int nn;
  G4double TotWd=0.;
  // Total Width
  for (nn=0; nn < nchan; nn++)TotWd = TotWd + br[nn];
  ntuple->column("totwd",TotWd);
  for (nn=0; nn < nchan; nn++){
  // fill decay channel histos
  switch (nn){
  case 0:
    ntuple->column("wd0",br[nn]);
    ntuple->column("br0",br[nn]/TotWd);
    break;
  case 1:
    ntuple->column("wd1",br[nn]);
    ntuple->column("br1",br[nn]/TotWd);
    break;
  case 2:
    ntuple->column("wd2",br[nn]);
    ntuple->column("br2",br[nn]/TotWd);
    break;
  case 3:
    ntuple->column("wd3",br[nn]);
    ntuple->column("br3",br[nn]/TotWd);
    break;
  case 4:
    ntuple->column("wd4",br[nn]);
    ntuple->column("br4",br[nn]/TotWd);
    break;
  case 5:
    ntuple->column("wd5",br[nn]);
    ntuple->column("br5",br[nn]/TotWd);
    break;
  case 6:
    ntuple->column("wd6",br[nn]);
    ntuple->column("br6",br[nn]/TotWd);
    break;
  case 7:
    ntuple->column("wd7",br[nn]);
    ntuple->column("br7",br[nn]/TotWd);
    break;
  case 8:
    ntuple->column("wd8",br[nn]);
    ntuple->column("br8",br[nn]/TotWd);
    break;
  case 9:
    ntuple->column("wd9",br[nn]);
    ntuple->column("br9",br[nn]/TotWd);
    break;
  case 10:
    ntuple->column("wd10",br[nn]);
    ntuple->column("br10",br[nn]/TotWd);
    break;
  }
  }
  //  ntuple->column("id",i);
  ntuple->dumpData();

  iter++;
  }  
  hbookManager->write();

  //  delete source;

  return EXIT_SUCCESS;
}
















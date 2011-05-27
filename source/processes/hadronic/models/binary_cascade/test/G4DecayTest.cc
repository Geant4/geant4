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
// $Id: G4DecayTest.cc,v 1.3 2009-12-16 17:50:33 gunter Exp $
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
  hbookManager = new HBookFile("decay.hbook", 58);
  assert (hbookManager != 0);

  // MGP ---- Book histograms
 
  HepHistogram* hSigma;

  hSigma = hbookManager->histogram("CrossSection", 100,0.,10.);
  assert (hSigma != 0);  

  // MGP ---- Book a ntuple
  HepTuple* ntuple;
  ntuple = hbookManager->ntuple("Decay ntuple");
  assert (ntuple != 0); 

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

  G4ParticleDefinition* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4ParticleDefinition* theMuonMinus = G4MuonMinus::MuonMinusDefinition();

  G4ParticleDefinition* theNeutrinoMu = G4NeutrinoMu::NeutrinoMuDefinition();
  G4ParticleDefinition* theAntiNeutrinoMu = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  G4ParticleDefinition* lambdacPlus = G4LambdacPlus::LambdacPlusDefinition();

  //  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //  G4ParticleDefinition* deltaPlus = particleTable->FindParticle("delta+");
  //  G4ParticleDefinition* deltaPlus = particleTable->FindParticle("N(1440)0");

  G4ShortLivedConstructor ShortLived;
  G4ShortLivedTable ShortLivedTab;
  ShortLived.ConstructParticle();

  G4int DEBUG;
  cout << "insert Debug value 0/1" << G4endl;
  G4cin >> DEBUG;

  // Select track 1
  
  G4cout << "---- Set KineticTrack 1 ----" << G4endl;
  G4int id1 = 0;
  G4cout << "1 p   2 ap   3 n   4 an   5 pi+   6 pi-   7 K+   8 K-  " << G4endl;
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
      //    case 9:
      //      definition1 = deltaPlus;
      //      break;
    default:
      definition1 = proton;
      break;
    }
  if(DEBUG==1)cout << definition1 << G4endl;
  G4double mass1 = definition1->GetPDGMass();
  G4String name1;
  if(DEBUG==1)name1 = definition1->GetParticleName();
  if(DEBUG==1)cout << "Particle Type: " << name1 << G4endl;
  if(DEBUG==1)cout << "Particle Mass: " << mass1 << G4endl;


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
  G4ThreeVector p1(px1,py1,pz1);
  G4double e1 = std::sqrt(mass1*mass1 + p1.mag()*p1.mag());
  G4LorentzVector p41(p1,e1);

  G4KineticTrack mytrack;

  G4int iter=0, iter_max;
  cout << " Enter number of iteration" << G4endl;
  G4cin >> iter_max;
  while(iter < iter_max){
    cout << " =================>   Processing event:......" << (iter+1) << G4endl;
  mytrack=G4KineticTrack(definition1,time1,pos1,p41);
  if(DEBUG==1)cout << mytrack.GetDefinition()->GetParticleName() << G4endl;
  //  mytrack.SetFormationTime(time1);
  //  mytrack.SetPosition(pos1);
  //  mytrack.Set4Momentum(p41);

  //  G4DecayTable Dectab=G4DecayTable();
  //  definition1->SetDecayTable(&Dectab);
  if(DEBUG==1){
   cout << "\n " << G4endl;
   definition1->DumpTable();
   cout << "\n " << G4endl;}
  G4DecayTable* TheDec=definition1->GetDecayTable();
  //  cout << " n decay channels: " << TheDec << G4endl;
  //  if(TheDec == NULL)cout << " **************** " << G4endl;

  //  cout << " Table not yet dumped " << G4endl;
  // TheDec->DumpInfo();
  //  cout << " Table dumped " << G4endl;

  G4KineticTrackVector* Dec = mytrack.Decay();
  if(Dec != NULL){
//    if(DEBUG==1)TheDec->DumpInfo();
    if(DEBUG==1)cout << " Pointer to Decay Products: "<< Dec << G4endl;  
    G4int nchan = TheDec->entries();
    /*    G4double TheBR[10];
    cout << nchan << G4endl;
    for(int index=nchan-1; index >=0; index--){
      TheBR[index]=TheDec->GetDecayChannel(index)->GetBR();
      cout << " BR: " << TheBR[index] << G4endl;} */
  } 

  // get kinetictracks from kinetictracksvector

  G4int nComponents = 0;
  if (Dec != 0)
    nComponents = Dec->entries();
  else
    continue;
  G4int i;
  G4LorentzVector Sumv=G4LorentzVector();
      G4double pxx1=-999.;
      G4double pyy1=-999.;
      G4double pzz1=-999.;
      G4double pxx2=-999.;
      G4double pyy2=-999.;
      G4double pzz2=-999.;
      G4double pxx3=-999.;
      G4double pyy3=-999.;
      G4double pzz3=-999.;
  for (i=0; i<nComponents; i++)
    {
      if(DEBUG==1)G4cout << "* Component " << i << ": " <<
	Dec->at(i)->GetDefinition()->GetParticleName() << "  " <<
	Dec->at(i)->GetPosition() << "  " <<
	Dec->at(i)->Get4Momentum() << "  " <<
	G4endl;
      
      switch (i) {
      case 0:
        pxx1=Dec->at(i)->Get4Momentum().px();
        pyy1=Dec->at(i)->Get4Momentum().py();
        pzz1=Dec->at(i)->Get4Momentum().pz();
        break;
      case 1:
        pxx2=Dec->at(i)->Get4Momentum().px();
        pyy2=Dec->at(i)->Get4Momentum().py();
        pzz2=Dec->at(i)->Get4Momentum().pz();
        break;
      case 2:
        pxx3=Dec->at(i)->Get4Momentum().px();
        pyy3=Dec->at(i)->Get4Momentum().py();
        pzz3=Dec->at(i)->Get4Momentum().pz();
        break;}      
        Sumv=Sumv+Dec->at(i)->Get4Momentum();
    }
  G4double mass=Sumv.m()/1000.;
  if(DEBUG==1)cout << "Mass:" << mass << G4endl;

  ntuple->column("mass",mass);
  ntuple->column("px1",pxx1);
  ntuple->column("px2",pxx2);
  ntuple->column("px3",pxx3);
  ntuple->column("py1",pyy1);
  ntuple->column("py2",pyy2);
  ntuple->column("py3",pyy3);
  ntuple->column("pz1",pzz1);
  ntuple->column("pz2",pzz2);
  ntuple->column("pz3",pzz3);
  //  ntuple->column("id",i);
  ntuple->dumpData();

  iter++;
  int k;
  }  
  hbookManager->write();

  //  delete source;

  return EXIT_SUCCESS;
}
















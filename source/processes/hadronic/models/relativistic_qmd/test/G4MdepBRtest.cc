// $Id: G4MdepBRtest.cc,v 1.1 2000-03-23 16:45:50 hweber Exp $
// -------------------------------------------------------------------
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
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
#include <fstream.h>
#include <iomanip.h>
#include <iostream.h>
#include <assert.h>
#include <vector>

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
#include "G4KaonZero.hh"
#include "G4KaonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Eta.hh"
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

#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

int main()
{
  // MGP ---- HBOOK initialization
  //  HepTupleManager* hbookManager;
  //  hbookManager = new HBookFile("decay.hbook", 58);
  //  assert (hbookManager != 0);

  // MGP ---- Book histograms
 
  //  HepHistogram* hSigma;

  //  hSigma = hbookManager->histogram("CrossSection", 100,0.,10.);
  //  assert (hSigma != 0);  

  // MGP ---- Book a ntuple
  //  HepTuple* ntuple;
  //  ntuple = hbookManager->ntuple("Decay ntuple");
  //  assert (ntuple != 0); 

  G4BosonConstructor Bosons;
  Bosons.ConstructParticle();

  G4MesonConstructor Mesons;
  Mesons.ConstructParticle();

  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();

  G4BaryonConstructor Barions;
  Barions.ConstructParticle();

  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  /*
  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();

  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiProton = G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();   
  G4ParticleDefinition* antiNeutron = G4AntiNeutron::AntiNeutronDefinition();   

  G4ParticleDefinition* pionPlus = G4PionPlus::PionPlusDefinition();
  G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();
  G4ParticleDefinition* pionZero = G4PionZero::PionZeroDefinition();

  G4ParticleDefinition* eta = G4Eta::EtaDefinition();
  
  */
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* N1440_0 = particleTable->FindParticle("N(1440)0");
  G4ParticleDefinition* N1535_m = particleTable->FindParticle("N(1535)-");
  G4ParticleDefinition* N1535_0 = particleTable->FindParticle("N(1535)0");
  G4ParticleDefinition* N1535_1 = particleTable->FindParticle("N(1535)+");
  G4ParticleDefinition* N1535_2 = particleTable->FindParticle("N(1535)++");
  G4ParticleDefinition* N1710_0 = particleTable->FindParticle("N(1710)0");
  /*
  G4ParticleDefinition* kaonPlus = G4KaonPlus::KaonPlusDefinition();
  G4ParticleDefinition* kaonZero = G4KaonZero::KaonZeroDefinition();
  G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();

  G4ParticleDefinition* rhoPlus = G4RhoPlus::RhoPlusDefinition();
  G4RhoZero::RhoZeroDefinition();
  G4ParticleDefinition* rhoMinus = G4RhoMinus::RhoMinusDefinition();
   
  G4ParticleDefinition* lambda = G4Lambda::LambdaDefinition();

  G4ParticleDefinition* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4ParticleDefinition* theMuonMinus = G4MuonMinus::MuonMinusDefinition();

  G4ParticleDefinition* theNeutrinoMu = G4NeutrinoMu::NeutrinoMuDefinition();
  G4ParticleDefinition* theAntiNeutrinoMu = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  G4ParticleDefinition* lambdacPlus = G4LambdacPlus::LambdacPlusDefinition();
  */

  G4int DEBUG;
  cout << "Insert DEBUG value (0|1)" << endl;
  cin >> DEBUG;

  // Select track 1
  
  G4cout << "---- Select KineticTrack 1 ----" << endl;
  G4int id1 = 0;
  G4cout << "1 N1535_0      2 N1440_0    3 N1535_1  4 N1710_0" << endl;
  cin >> id1;

  G4ParticleDefinition* definition1;
  switch (id1) 
    { 
    case 1:
      definition1 = N1535_0; break;
    case 2:
      definition1 = N1440_0; break;
    case 3:
      definition1 = N1535_1; break;
    case 4:
      definition1 = N1710_0; break;
    default:
      definition1 = N1535_0; break;
    }
  if (DEBUG) cout << definition1 << endl;
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
  cout << " Enter Pz (Px=Py=0.)" << endl;
  cin >> pz1;
  G4ThreeVector p1(px1,py1,pz1);

  G4double startMass = definition1->GetPDGMass()-definition1->GetPDGWidth()*6;
  G4double endMass = definition1->GetPDGMass()+definition1->GetPDGWidth()*6;
  G4String name1 = definition1->GetParticleName();
  const G4int nSteps=100;
  G4double stepMass = (endMass-startMass)/double(nSteps);
  typedef vector<double> graph;
  G4int nChannels = definition1->GetDecayTable()->entries();
  cout << "nChannels: " << nChannels << endl;
  vector<graph> plot(nChannels);
  if (DEBUG) cout << "Particle Type: " << name1 << endl;
  cout << "! title \"Mass Dependent Decay Width (" << definition1->GetParticleName() << ")\"" << endl;
  cout << "! xaxis \"mass [GeV]\" " << startMass << " " << endMass << endl;
  cout << "! yaxis \"Ratio\" 0  0.8" << endl;
  for (G4int i=0; i<nChannels; i++) {
    cout << "! col " << i << " \"";
    G4int nDaughters=definition1->GetDecayTable()->GetDecayChannel(i)->GetNumberOfDaughters();
    cout << definition1->GetDecayTable()->GetDecayChannel(i)->GetDaughterName(0);
    for (G4int d=1; d<nDaughters; d++) {
      cout << " + " << definition1->GetDecayTable()->GetDecayChannel(i)->GetDaughterName(d);
    } 
    cout << "\"" << endl;
  }
  cout << "! col " << i << " \"Breit-Wigner\" fill" << endl;

  for (G4int step=0; step<nSteps; step++) 
    {    
      G4double mass1 = startMass + step*stepMass;
      G4double e1    = sqrt(mass1*mass1 + p1.mag()*p1.mag());
      G4LorentzVector p41(p1,e1);

      //      if (DEBUG) cout << "Step " << step << "  Particle Mass: " << mass1 << endl;
      
      G4KineticTrack* mytrack = new G4KineticTrack(definition1,time1,pos1,p41);
      G4double* awidth = mytrack->GetActualWidth();

      //      cout << "ActualMass: " << endl;

      cout << mass1 << " ";
      G4double TotWidth=0.0;
      for (G4int i=0; i<nChannels; i++) {
	TotWidth+=awidth[i];
	//	cout << TotWidth << endl;
      }
      for (i=0; i<nChannels; i++) {
	cout << awidth[i]/TotWidth << " ";
      }
      //            cout << endl;
      G4double BreitWigner = 0.7*mytrack->BrWig(definition1->GetPDGWidth(),mass1,definition1->GetPDGMass()) / mytrack->BrWig(definition1->GetPDGWidth(),definition1->GetPDGMass(),definition1->GetPDGMass());
      cout << BreitWigner << endl;
      delete mytrack;
    }
  if (DEBUG) { 
    //    cout << endl;
    definition1->DumpTable();
    cout << endl;
  }
  //  for (G4int i=0; i<nChannels; i++) {
  // cout << i << plot[i] << endl;
  // }
  //  G4DecayTable* TheDec=definition1->GetDecayTable();
  //  cout << " n decay channels: " << TheDec << endl;
  //  if (TheDec == NULL) cout << " **************** " << endl;
  


  //  G4KineticTrackVector* Dec = mytrack->Decay();
  //  if(Dec != NULL){
    //    if(DEBUG==1)TheDec->DumpInfo();
  //   if(DEBUG==1)cout << " Pointer to Decay Products: "<< Dec << endl;  
  //G4int nchan = TheDec->entries();
    /*    G4double TheBR[10];
	  cout << nchan << endl;
	  for(int index=nchan-1; index >=0; index--){
	  TheBR[index]=TheDec->GetDecayChannel(index)->GetBR();
	  cout << " BR: " << TheBR[index] << endl;} */
  //} 
  
  

  // get kinetictracks from kinetictracksvector
  /*
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
		      endl;
	
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
    if(DEBUG==1)cout << "Mass:" << mass << endl;
    
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
    cout << iter << endl;
    iter++;
    cout << iter << endl;
    int k;
  }  
  hbookManager->write();
  
  //  delete source;
  */
  return EXIT_SUCCESS;
}
















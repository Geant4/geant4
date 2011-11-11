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
#include "G4AtomicTransitionManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <vector>
#include "G4DynamicParticle.hh"

#ifdef G4ANALYSIS_BUILD
#include "AIDA/AIDA.h"
#endif

#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4He3.hh"



using namespace CLHEP;

int main(int argc, char* argv[]){

  if (!argc) argc=0;

  time_t seconds = time(NULL);
  G4int seed = seconds;
  
  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);

  G4int Z;
  G4int a;
  G4int b;
  G4int startId;
  G4int vacancyIndex;
  G4int numberOfRun;
  G4String PIXEmodel;
  G4int batch=0;
  G4int element = 0;
  if (argv[1]) {batch = atoi(argv[1]);}
  G4String fileName;
  if (argv[3]) {element = atoi(argv[3]);}
  if (argv[4]) {fileName = argv[4];}
  else {fileName = "transitions.xml";}

#ifdef G4ANALYSIS_BUILD
  AIDA::ITree* tree;
  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITupleFactory* tupleFactory;
  AIDA::ITuple* tupleFluo = 0;
#endif

  if (batch != 1) {
    G4cout << "Enter Z " << G4endl;
    G4cin >> a;
    G4cout << "Enter the index of the vacancy" << G4endl;
    G4cin >> startId;
    G4cout<<"Enter model for PIXE XS ( Empirical / Analytical / ECPSSR_FormFactor )"<<G4endl;
    G4cin>> PIXEmodel;
    G4cout<<"Enter the number of runs "<<G4endl;
    G4cin>> numberOfRun;
  }
  else {
    
    a = 0;
    startId = -1;
    numberOfRun = atoi(argv[2]);
  }

#ifdef G4ANALYSIS_BUILD
  analysisFactory = AIDA_createAnalysisFactory();
  AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
  tree = treeFactory->create(fileName,"xml",false,true);
  tupleFactory = analysisFactory->createTupleFactory(*tree);
  // Book tuple column names
  std::vector<std::string> columnNames;
  // Book tuple column types
  std::vector<std::string> columnTypes;
#endif
  
  //if Z=0 a number of runs numberOfRun is generated for all the elements 
  if (a==0)
    {
      if (element == 0) { 
	a = 6;
	b = 98;
      }
      else {
	a = element;
	b = a;}
#ifdef G4ANALYSIS_BUILD
      columnNames.push_back("AtomicNumber");
      columnNames.push_back("Particle");
      columnNames.push_back("Energies");
      
      columnTypes.push_back("int");
      columnTypes.push_back("int");
      columnTypes.push_back("double");
      tupleFluo = tupleFactory->create("10", "Total Tuple", columnNames, columnTypes, "");
      assert(tupleFluo);
#endif
    }
  else { b = a;} 
  
  G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
 
  G4UAtomicDeexcitation* deexcitation = new G4UAtomicDeexcitation;
 

 std::map<G4int,G4int> shellNumberTable;
  
 deexcitation->SetPIXECrossSectionModel(PIXEmodel);
 deexcitation->SetAuger(true);
 deexcitation->SetPIXE(true);
 deexcitation->InitialiseForNewRun();    
  
  for (Z = a; Z<=b; Z++) {    
  G4cout << "******** Z = "<< Z << "*********" << G4endl;

  G4int numberOfPossibleShell = transitionManager->NumberOfShells(Z);

  shellNumberTable[Z] = numberOfPossibleShell;
  G4int min = 0;
  G4int max = 0;
  std::vector<G4DynamicParticle*>* vectorOfParticles = new std::vector<G4DynamicParticle*>;


  G4ParticleDefinition* particle;

  assert(vectorOfParticles);

    for(G4int i = 0; i<numberOfRun;i++){ 
      G4cout<<"**************"<<G4endl;
      G4cout<<"begin of run "<< i <<G4endl;
      G4cout<<"**************"<<G4endl;
      vectorOfParticles->clear();
      // if startId = -1 the test runs on every shell of the atom
      if (startId == -1){
	min = 0;
	max = shellNumberTable[Z];
      }
      else {
	min = startId;
	max = min;
      }

      for (vacancyIndex = min; vacancyIndex <= max; vacancyIndex++) { 

	G4AtomicShell* shell = transitionManager->Shell(Z, vacancyIndex);
	G4AtomicShellEnumerator as;

	if (shell->ShellId() == 1) {as = fKShell;}
	else if (shell->ShellId() == 3 ) {as = fL1Shell;}
	else if (shell->ShellId() == 5 ) {as = fL2Shell;}
	else if (shell->ShellId() == 6 ) {as = fL3Shell;}
	else if (shell->ShellId() == 8 ) {as = fM1Shell;}
	else if (shell->ShellId() == 10) {as = fM2Shell;}
	else if (shell->ShellId() == 11) {as = fM3Shell;}
	else if (shell->ShellId() == 13) {as = fM4Shell;}
	else if (shell->ShellId() == 14) {as = fM5Shell;}

	// testing how to get shell from deexcitation
	const G4AtomicShell* shell2 = deexcitation->GetAtomicShell(Z, as);

	if (shell2 == shell) {

	  G4cout << "**********************************" << G4endl;
	  G4cout << "deexcitation->GetAtomicShell WORKS" << G4endl;
	  G4cout << "**********************************" << G4endl;
	}

	//loop over the energy? no, let's try the "standard" ones
      
	particle = G4Proton::Proton();

      	G4double crossSecProton = deexcitation->GetShellIonisationCrossSectionPerAtom
	  (particle,Z,as,3.0 * MeV) /* * barn*/ ;

	particle = G4Alpha::Alpha();

      	G4double crossSecAlpha = deexcitation->GetShellIonisationCrossSectionPerAtom
	  (particle, Z, as, 5.8 * MeV) /* * barn */;

	particle = G4He3::He3();

      	G4double crossSecHe3 = deexcitation->GetShellIonisationCrossSectionPerAtom
	  (particle, Z, as, 5.8 * MeV) /* * barn */;
      
	G4cout << "Shell ID: " << shell->ShellId() << G4endl; 
	G4cout<<  vectorOfParticles->size()<<" particles in the vector "<<G4endl;
	G4cout<< as <<" XS for p @ 3 MeV: "<< crossSecProton/barn  << "barns" << G4endl;
	G4cout<< as <<" XS for alpha @ 5.8 MeV: "<< crossSecAlpha/barn << "barns" << G4endl; 
	G4cout<< as <<" XS for He3 @ 5.8 MeV: "<< crossSecHe3/barn << "barns" << G4endl; 

	deexcitation->GenerateParticles(vectorOfParticles, shell, Z, 0, 0);

	for (G4int k=0; k< vectorOfParticles->size();k++)
	  {

	    G4DynamicParticle* newParticle = (*vectorOfParticles)[k];

	    if ( newParticle->GetDefinition()->GetParticleName() == "e-")
	      {

		G4DynamicParticle* newElectron = (*vectorOfParticles)[k];
		G4ThreeVector augerDirection =newElectron ->GetMomentum();
		G4double  augerEnergy =newElectron ->GetKineticEnergy();


		if (startId==-1){
#ifdef G4ANALYSIS_BUILD		  
		  tupleFluo->fill(0,Z);
		  tupleFluo->fill(1,0);
		  tupleFluo->fill(2,augerEnergy);
		  tupleFluo->addRow();
#endif		  
		}
		else{	      
		  
		  G4cout <<" An auger has been generated"<<G4endl;
		  G4cout<<" vectorOfParticles ["<< k <<"]:"<<G4endl;
		  G4cout<<"Non zero particle. Index: "<< k <<G4endl;
		  G4cout<< "The Auger electron has a kinetic energy = "<<augerEnergy
			<<" MeV " <<G4endl;
		}
	      }
	    else{
	      //G4cout << "pippo" << G4endl; 

	      G4ThreeVector photonDirection = newParticle ->GetMomentum();
	      G4double  photonEnergy =newParticle ->GetKineticEnergy();

	      if (startId==-1){
#ifdef G4ANALYSIS_BUILD
		//G4cout << "pippo2" << G4endl;
		tupleFluo->fill(0,Z);

		tupleFluo->fill(1,1);
		tupleFluo->fill(2,photonEnergy);
		tupleFluo->addRow();		
#endif
	      }
	      else{
		
		G4cout<<" vectorOfParticles ["<<k<<"]:"<<G4endl;
		G4cout<<"Non zero particle. Index: "<<k<<G4endl;
		G4cout<< "The photon has a kinetic energy = "<<photonEnergy
		      <<" MeV " <<G4endl;
	      }
	    }
	  }
      }
      if (batch == 1){
#ifdef G4ANALYSIS_BUILD
	tree->commit(); // Write histos in file. 
	tree->close();
#endif
      }
      if (!vectorOfParticles) delete vectorOfParticles;
    }
  }  
  delete deexcitation;
  G4cout<<"END OF THE MAIN PROGRAM"<<G4endl;
}

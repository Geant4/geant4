// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyPhotoElectric.cc,v 1.12 1999-06-07 09:59:15 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4LowEnergyPhotoElectric physics process --------
//                   by Michel Maire, April 1996
// **************************************************************
// 12-06-96, Added SelectRandomAtom() method, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, PartialSumSigma(i)
//           split of ComputeBindingEnergy, M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 04-06-98, in DoIt, secondary production condition: range>min(threshold,safety)
// --------------------------------------------------------------

// This Class Header
#include "G4LowEnergyPhotoElectric.hh"

// Collaborating Class Headers
#include "G4EnergyLossTables.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Epdl97File.hh"
#include "G4Epdl89File.hh"
#include "G4EpdlTables.hh"
#include "G4PhysicsFreeVector.hh" 
#include "G4Step.hh"
#include "Randomize.hh" 

#include "CLHEP/String/Strings.h"

#include <rw/tpordvec.h>
typedef RWTPtrOrderedVector<G4DynamicParticle> G4ParticleVector;
// constructor
 
G4LowEnergyPhotoElectric::G4LowEnergyPhotoElectric(const G4String& processName)
  : G4VDiscreteProcess(processName),             // initialization
    LowestEnergyLimit (100*eV),
    HighestEnergyLimit(100*GeV),
    theCrossSectionTable(0),
    theBindingEnergyTable(0),
    theMeanFreePathTable(0),
    theFluorTransitionTable(0),
    theAugerTransitionTable(0),
    NumbBinTable(200)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
     G4cout << "LowestEnergy: " << LowestEnergyLimit/keV << "keV ";
     G4cout << "HighestEnergy: " << HighestEnergyLimit/MeV << "MeV " << endl;
   }
}
 
// destructor
 
G4LowEnergyPhotoElectric::~G4LowEnergyPhotoElectric()
{
   if (theCrossSectionTable) {
      theCrossSectionTable->clearAndDestroy();
      delete theCrossSectionTable;
   }

   if (theBindingEnergyTable) {
      theBindingEnergyTable->clearAndDestroy();
      delete theBindingEnergyTable;
   }

   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }

   // ClearAndDestroy of this tables is called in their destructors
   if (theFluorTransitionTable) {

      delete theFluorTransitionTable;
   }

   if (theAugerTransitionTable) {

      delete theAugerTransitionTable;
   }

}
 
 
// methods.............................................................................
 
void G4LowEnergyPhotoElectric::BuildPhysicsTable(const G4ParticleDefinition& PhotonType)

// Build microscopic cross section table and mean free path table
{

  BuildCrossSectionTable();

  BuildMeanFreePathTable();
  
  BuildBindingEnergyTable();

  BuildFluorTransitionTable();

   
}

G4VParticleChange* G4LowEnergyPhotoElectric::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){

  // incoming particle initialization
  if(getenv("GENERAL")) aParticleChange.Initialize(aTrack);

  G4Material* aMaterial = aTrack.GetMaterial();

  const G4DynamicParticle* aDynamicPhoton = aTrack.GetDynamicParticle();
  const G4double PhotonEnergy = aDynamicPhoton->GetKineticEnergy();
  const G4ParticleMomentum PhotonDirection = aDynamicPhoton->GetMomentumDirection();
   
  // select randomly one element constituing the material.
  G4Element* anElement = SelectRandomAtom(aDynamicPhoton, aMaterial);

  // PAY ATTENTION TO THE MEANING OF THIS NUMBER
  G4int AtomNum = (G4int) anElement->GetZ();

  // Select the subshell ionized the first one 
  // with binding energy less than the photon energy

  G4PhysicsVector* theBindEnVec = (*theBindingEnergyTable)(AtomNum-1);
  
  G4int g = 0;
  while(g < theBindEnVec->GetVectorLength() && PhotonEnergy < (*theBindEnVec)(g)){
    g++;
  }

  G4int thePrimaryShell = (G4int) theBindEnVec->GetLowEdgeEnergy(g);

  // Create lists of pointers to DynamicParticles (photons and electrons)
  G4ParticleVector photvec;
  G4int photInd = 0; 
  G4ParticleVector elecvec;
  G4int elecInd = 0; 

  // primary outcoming electron
  G4double ElecKineEnergy = (PhotonEnergy - (*theBindEnVec)(g))*MeV;
  G4double theEnergyDeposit = (PhotonEnergy - ElecKineEnergy)*MeV;

  if (G4EnergyLossTables::GetRange(G4Electron::Electron(),ElecKineEnergy,aMaterial)
      >= min(G4Electron::GetCuts(), aStep.GetPostStepPoint()->GetSafety()) ){

    // the electron is created in the direction of the incident photon ...  
    G4DynamicParticle* aElectron = new G4DynamicParticle (G4Electron::Electron(), 
							  PhotonDirection, ElecKineEnergy) ;
    elecvec.append(aElectron);

    // load the transition probability table for the element
    // theTable[i][j][k] 
    // i = subshell, j = type of information (second shell, transition energy , 
    // transition probability), k = previous vectors.

    if(AtomNum > 5){
      
      G4bool ThereAreShells = TRUE;
      oneAtomTable* oneAtomFluorTrans = (*theFluorTransitionTable)[AtomNum-6]; 
      
      while(ThereAreShells == TRUE){
	
	// Select the second transition from another subshell
	// fluorPar[0] = SubShell 
	// fluorPar[1] = Sec SubShell (if there is), 
	// fluorPar[2] = Transition Probability
	// fluorPar[3] = Transition Energy
	// the same for augerPar
	
	G4double fluorPar[4] = {0};
	ThereAreShells = SelectRandomTransition(thePrimaryShell, fluorPar, oneAtomFluorTrans);
	
	// Daugther dynamic particle
	G4DynamicParticle* newPart;
	
	// Direction of the outcoming particle isotropic selection
	G4double newcosTh = 1-2*G4UniformRand();
	G4double newsinTh = sqrt(1-newcosTh*newcosTh);
	G4double newPhi = twopi*G4UniformRand();
	
	G4double dirx, diry, dirz;
	dirz = newcosTh;
	diry = newsinTh*cos(newPhi);
	dirx = newsinTh*sin(newPhi);
	G4ThreeVector newPartDirection(dirx, diry, dirz);
	newPartDirection.rotateUz(PhotonDirection);

	if(ThereAreShells != FALSE){
	  
	  thePrimaryShell = (G4int) fluorPar[0];
	  newPart = new G4DynamicParticle (G4Gamma::Gamma(), newPartDirection, fluorPar[3]) ;
	  photvec.append(newPart);
	  theEnergyDeposit -= fluorPar[3]*MeV;
	  
	}
	else{
	  
	  
	  G4int k = 0;
	  while(thePrimaryShell != theBindEnVec->GetLowEdgeEnergy(k)) k++;
	  
	  G4double lastTransEnergy = (*theBindEnVec)(k);
	  
	  newPart = new G4DynamicParticle (G4Gamma::Gamma(), newPartDirection, lastTransEnergy) ;
	  photvec.append(newPart);
	  thePrimaryShell = (G4int) fluorPar[0];
	  theEnergyDeposit -= lastTransEnergy*MeV;
	  
	}
      }
     
    } //END OF THE CHECK ON ATOMIC NUMBER
    
    //controllare se il setnumberofsecondaries  si puo' cambiare
    G4int numOfElec = elecvec.entries(), numOfPhot = photvec.entries();
    G4int numOfDau = numOfElec + numOfPhot;

    aParticleChange.SetNumberOfSecondaries(numOfDau);
    G4int l = 0;
    for( l = 0; l<numOfElec; l++ ){
      
      aParticleChange.AddSecondary(elecvec[l]);
    }
    
    for(l = 0; l < numOfPhot; l++) {

      aParticleChange.AddSecondary(photvec[l]); 
    }
    
    photvec.clear();
    elecvec.clear();
  } // END OF CUTS
  
  else{
    
    ElecKineEnergy = 0. ;
    aParticleChange.SetNumberOfSecondaries(0) ;
  }

  // Kill the incident photon 
  aParticleChange.SetMomentumChange( 0., 0., 0. );
  aParticleChange.SetEnergyChange( 0. );

  // the energy deposit with fluorescence made like this is ZERO.
  if(theEnergyDeposit < 0){
    theEnergyDeposit = 0;
  }

  aParticleChange.SetLocalEnergyDeposit(theEnergyDeposit) ;  
  aParticleChange.SetStatusChange( fStopAndKill ) ; 
  
  // Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

void G4LowEnergyPhotoElectric::BuildCrossSectionTable(){
 
  if (theCrossSectionTable) {
    
    theCrossSectionTable->clearAndDestroy(); 
    delete theCrossSectionTable; 
  }

  G4int par[4] = {73, 0, 0, 0}; G4String name("epdl97");
  G4Epdl97File File(name,par);
  G4EpdlTables table(File);
  table.FillDataTable();
  theCrossSectionTable = table.GetFstDataTable();

}

void G4LowEnergyPhotoElectric::BuildBindingEnergyTable(){
 
  if (theBindingEnergyTable) {
    
    theBindingEnergyTable->clearAndDestroy(); 
    delete theBindingEnergyTable; 
  }

  G4int par[4] = {91, 913, 0, 0}; G4String name("eadl.asc");
  G4Epdl89File File(name,par);
  G4EpdlTables table(File);
  table.FillDataTable();
  theBindingEnergyTable = table.GetFstDataTable();

}

void G4LowEnergyPhotoElectric::BuildFluorTransitionTable(){

   if (theFluorTransitionTable) {

    delete theFluorTransitionTable;
  }

   theFluorTransitionTable = new allAtomTable();
   G4int dataNum = 3;
 
   for(G4int TableInd = 5; TableInd < 100; TableInd++){

     oneAtomTable* oneAtomShellFL = BuildTables(TableInd, dataNum, "fl-tr-pr-");
     
     theFluorTransitionTable->insert(oneAtomShellFL);

   }//end for on atoms
}

oneAtomTable* G4LowEnergyPhotoElectric::BuildTables(const G4int TableInd, 
						    const G4int ParNum, 
						    const char* prename){

  HepString Znum(TableInd+1);
  HepString name = prename + Znum + ".dat"; 
  
  char* path = getenv("G4LEDATA");
  if(!path){ 
    
    G4Exception("G4LEDATA environment variable not set");
  }
  
  HepString path_string(path);
  HepString dir_file = path_string + "/" + name;
  ifstream file(dir_file);
  filebuf* lsdp = file.rdbuf();
  
  if(!lsdp->is_open()){
    
      HepString excep = "Error!!!! data file: " + dir_file + " NOT found";
      G4Exception(excep);
  }
  
  oneAtomTable* oneAtomPar = new oneAtomTable();
  oneShellTable* oneShellPar = new oneShellTable();
  
  for(G4int j = 0; j < ParNum; j++){ 
    
    oneShellPar->insertAt(j,new G4Data());
  }
  
  G4double a = 0;
  G4int k = 1, s = 0;
  
  do{
    
    file>>a;
    
    if(a == -1){
      
      if(s == 0){
	
	oneAtomPar->insert(oneShellPar);
	oneShellPar = new oneShellTable();
	
	for(G4int j = 0; j < ParNum; j++){ 
	  
	  oneShellPar->insertAt(j,new G4Data());
	}
      }
      
      s++;
      
      if(s == ParNum){
	
	s = 0;
      }
    }

    else if(a == -2){
      
      delete oneShellPar;
    }

    else{
      
      if(k%ParNum != 0){	
	
	(*oneShellPar)[k-1]->insert(a);
	k++;
      }
      else if(k%ParNum == 0){
	
	(*oneShellPar)[k-1]->insert(a);
	k = 1;
      }
    }

  }while(a != -2); //end for on file
  
  file.close();
  return oneAtomPar;
}

void G4LowEnergyPhotoElectric::BuildAugerTransitionTable(G4int  atomNum){
 
  if (theAugerTransitionTable) {
    
    //  theAugerTransitionTable->clearAndDestroy(); 
    delete theAugerTransitionTable; 
  }

  G4int par[4] = {92, 932, 91, 0}; G4String name("eadl.asc");
  G4Epdl89File File(name,par);
  G4EpdlTables table(File);
  theAugerTransitionTable = table.FillTheTable(atomNum);

}

void G4LowEnergyPhotoElectric::BuildMeanFreePathTable(){

  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy(); 
    delete theMeanFreePathTable; }

  // material
  G4double NumbOfMaterials = G4Material::GetNumberOfMaterials();
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;
  G4Material* material;

  // MeanFreePath 
  G4double LowEdgeEnergy, Value;
  theMeanFreePathTable = new G4PhysicsTable(NumbOfMaterials);
  G4PhysicsLogVector* ptrVector;

  for ( G4int J = 0 ; J < NumbOfMaterials; J++ ) { // For each material 
  
    //create physics vector then fill it ....
    // WARNING!!! below 50 ev cross section lower limit depend on the element  
    ptrVector = new  G4PhysicsLogVector(LowestEnergyLimit, HighestEnergyLimit, NumbBinTable);
    
    material = (*theMaterialTable)(J);
    const G4ElementVector* theElementVector = material->GetElementVector();
    const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();   
    
    for ( G4int i = 0 ; i < NumbBinTable ; i++ ){ 
      //For each energy
      
      LowEdgeEnergy = ptrVector->GetLowEdgeEnergy(i);
      
      const G4double BigPath= DBL_MAX;
      G4double SIGMA = 0 ;
      
      for ( G4int k=0 ; k < material->GetNumberOfElements() ; k++ ){ 
	// For each element            
	G4int tableIndex = (G4int) (*theElementVector)(k)->GetZ()-1;
	G4double interCrsSec = DataLogInterpolation(LowEdgeEnergy, tableIndex, theCrossSectionTable)*barn;
	SIGMA += theAtomNumDensityVector[k]*interCrsSec;
      }       
      
      Value = SIGMA<=0.0 ? BigPath : 1./SIGMA ;
      
      ptrVector->PutValue( i , Value ) ;
    }
    
    theMeanFreePathTable->insertAt( J , ptrVector ) ;
  }
}

G4Element*
G4LowEnergyPhotoElectric::SelectRandomAtom(const G4DynamicParticle* aDynamicPhoton, G4Material* aMaterial){

  // select randomly 1 element within the material
  G4double GammaEnergy = aDynamicPhoton->GetKineticEnergy();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  if (NumberOfElements == 1) return (*theElementVector)(0);

  const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();

  G4double PartialSumSigma = 0. ;

  // G4int materialIndex = aMaterial->GetIndex();
  //  MeanFreePath = DataLogInterpolation(GammaEnergy, materialIndex, theMeanFreePathTable);

  G4double rval = G4UniformRand()/MeanFreePath;
 
  for ( G4int i=0 ; i < NumberOfElements ; i++ ){ 

    G4double crossSection;
    if (GammaEnergy <  LowestEnergyLimit)

      crossSection = 0. ;

    else {

      if (GammaEnergy > HighestEnergyLimit) GammaEnergy = 0.99*HighestEnergyLimit ;

      G4int tableIndex = (G4int) (*theElementVector)(i)->GetZ()-1;
      crossSection = DataLogInterpolation(GammaEnergy, tableIndex, theCrossSectionTable)*barn;
    }

    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;

    if (rval <= PartialSumSigma) return ((*theElementVector)(i));

  }

  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
	 << "' has no elements" << endl;
  return (*theElementVector)(0);
}


G4bool G4LowEnergyPhotoElectric::SelectRandomTransition(G4int thePrimShell, 
							G4double* TransParam,
							const oneAtomTable* TransitionTable){
  
  G4int SubShellCol, SecShellCol, ProbCol, EnergyCol;

  // too check when the subshell are finished
  G4bool ColIsFull = TRUE;
  
  if(TransParam[0] == 0){  
    
    SubShellCol = 0; ProbCol = 1; EnergyCol = 2;
  }
  
  if(TransParam[0] == 1){  
    
    G4int SubShellCol = 0; SecShellCol = 1; ProbCol = 2; EnergyCol = 3;
  }
  
  G4int ShellNum = 0;
  G4double TotalSum = 0; 
  while(thePrimShell != (*(*(*TransitionTable)[ShellNum])[0])[0]){
  
    if(ShellNum == TransitionTable->entries()-1){
      break;
    }
    ShellNum++;
  }

  if(ShellNum != TransitionTable->entries()-1) {
    
    //TransProb start from 1 because the first element of the list is the primary shall id number
    G4int TransProb = 1;
    for(TransProb = 1; TransProb < (*(*TransitionTable)[ShellNum])[ProbCol]->length(); TransProb++){ 
      
      TotalSum += (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
    }
    
    G4double PartialProb = G4UniformRand()*TotalSum;
    G4double PartSum = 0;

    TransProb = 1; 
    while(TransProb < (*(*TransitionTable)[ShellNum])[ProbCol]->length()){
      
      PartSum += (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];

      if(PartialProb <= PartSum){
	
	TransParam[0] = (*(*(*TransitionTable)[ShellNum])[SubShellCol])[TransProb];

	// This if will be needed when Auger Effect will be added TransPar[i] initialized 
	// to 0 for fluorescence and to 1 for auger effect at the moment this distinction is 
	// tricky

	if(TransParam[1] == 0){
	  
	  TransParam[1] = 0;
	}
	else{
	  
	  TransParam[1] = (*(*(*TransitionTable)[ShellNum])[SecShellCol])[TransProb];
	}
	
	TransParam[2] = (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
	TransParam[3] = (*(*(*TransitionTable)[ShellNum])[EnergyCol])[TransProb];
	break;
      }
      
      TransProb++;
    }

    if(TransProb == (*(*TransitionTable)[ShellNum])[ProbCol]->length()-1) {

      ColIsFull = FALSE;
    }
  }
  else{
   
    ColIsFull = FALSE;
  }
  
  return ColIsFull;
}







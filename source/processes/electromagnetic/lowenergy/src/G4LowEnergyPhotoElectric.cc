// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyPhotoElectric.cc,v 1.16 1999-07-06 15:03:03 aforti Exp $
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
#include "G4Electron.hh"

typedef RWTPtrOrderedVector<G4DynamicParticle> G4ParticleVector;
// constructor
 
G4LowEnergyPhotoElectric::G4LowEnergyPhotoElectric(const G4String& processName)
  : G4VDiscreteProcess(processName),             // initialization
    LowestEnergyLimit (250*eV),
    HighestEnergyLimit(100*GeV),
    theCrossSectionTable(0),
    theBindingEnergyTable(0),
    theMeanFreePathTable(0),
    theFluorTransitionTable(0),
    allAtomShellCrossSec(0),
    CutForLowEnergySecondaryPhotons(0.),
    ZNumVec(0),
    ZNumVecFluor(0),
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
     //   theCrossSectionTable->clearAndDestroy();
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

   if (allAtomShellCrossSec) {

     delete allAtomShellCrossSec;
   }

   if(ZNumVec){
     
     ZNumVec->clear();
     delete ZNumVec;
   }
   
   if(ZNumVecFluor){

     ZNumVecFluor->clear();
     delete ZNumVecFluor;
   }

}
 
 
// methods.............................................................................
void G4LowEnergyPhotoElectric::SetCutForLowEnSecPhotons(G4double cut){

  CutForLowEnergySecondaryPhotons = cut;
}

void G4LowEnergyPhotoElectric::BuildPhysicsTable(const G4ParticleDefinition& PhotonType)

// Build microscopic cross section table and mean free path table
{

  BuildZVec();

  BuildCrossSectionTable();

  BuildShellCrossSectionTable();

  BuildMeanFreePathTable();
  
  BuildBindingEnergyTable();

  BuildFluorTransitionTable();

   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyPhotoElectric::BuildCrossSectionTable(){
 
  if (theCrossSectionTable) {
    
    //theCrossSectionTable->clearAndDestroy(); 
    delete theCrossSectionTable; 
  }

  theCrossSectionTable = new G4SecondLevel();
  G4int dataNum = 2;
 
  for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){

    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];

    G4FirstLevel* oneAtomCS = util.BuildFirstLevelTables(AtomInd, dataNum, "phot/pe-cs-");
     
     theCrossSectionTable->insert(oneAtomCS);
   
  }//end for on atoms
}

void G4LowEnergyPhotoElectric::BuildShellCrossSectionTable(){

   if (allAtomShellCrossSec) {

    delete allAtomShellCrossSec;
   }

   allAtomShellCrossSec = new allAtomTable();
   G4int dataNum = 2;
 
   for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){

     G4int AtomInd = (G4int) (*ZNumVec)[TableInd];

     oneAtomTable* oneAtomShellCS = util.BuildSecondLevelTables(AtomInd, dataNum, "phot/pe-ss-cs-");
     
     allAtomShellCrossSec->insert(oneAtomShellCS);
   
   }//end for on atoms
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyPhotoElectric::BuildBindingEnergyTable(){

  if (theBindingEnergyTable) {

    delete theBindingEnergyTable;
  }

  G4int dataNum = 2;
  theBindingEnergyTable = util.BuildSecondLevelTables(0,dataNum,"fluor/binding");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyPhotoElectric::BuildFluorTransitionTable(){

  if (theFluorTransitionTable) {
    
    delete theFluorTransitionTable;
  }
  
  theFluorTransitionTable = new allAtomTable();
  ZNumVecFluor = new G4Data(*ZNumVec);
  G4int dataNum = 3;
  
  for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    if(AtomInd > 5){
      
      oneAtomTable* oneAtomShellFL = util.BuildSecondLevelTables(AtomInd, dataNum, "fluor/fl-tr-pr-");
      theFluorTransitionTable->insert(oneAtomShellFL);
    }
    else{
      ZNumVecFluor->remove(AtomInd);
    }
  }//end for on atoms
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyPhotoElectric::BuildZVec(){

  const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

  if(ZNumVec){

    ZNumVec->clear();
    delete ZNumVec;
  }

  ZNumVec = new G4Data(); 
  for (G4int J=0 ; J < numOfMaterials; J++){ 
 
    const G4Material* material= (*theMaterialTable)[J];        
    const G4ElementVector* theElementVector = material->GetElementVector();
    const G4int NumberOfElements = material->GetNumberOfElements() ;

    for (G4int iel=0; iel<NumberOfElements; iel++ ){

      G4double Zel = (*theElementVector)(iel)->GetZ();

      if(ZNumVec->contains(Zel) == FALSE){

	ZNumVec->insert(Zel);
      }
      else{
	
	continue;
      }
    }
  }
}

G4double G4LowEnergyPhotoElectric::ComputeCrossSection(const G4double AtomIndex,
						    const G4double IncEnergy){
  // calculates the microscopic cross section from subshell cross sections
  //(it is called for elements , AtomicNumber = Z )
 
  G4double TotalCrossSection(0.);

  const oneAtomTable* oneAtomCS
    = (*allAtomShellCrossSec)[ZNumVec->index(AtomIndex)];

  for(G4int ind = 0; ind < oneAtomCS->entries(); ind++){

    G4double crossSec = 0;
    G4Data* EnergyVector = (*(*oneAtomCS)[ind])[0];
    G4Data* CrossSecVector = (*(*oneAtomCS)[ind])[1];

    if(IncEnergy < (*EnergyVector)[1]){ // First element is the shell number

      crossSec = 0;
    }

    else{

      crossSec = util.DataLogInterpolation(IncEnergy, (*EnergyVector), (*CrossSecVector));

    }

    TotalCrossSection += crossSec;
  }

  return TotalCrossSection ;
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
	G4int AtomIndex = (G4int) (*theElementVector)(k)->GetZ();
	const G4FirstLevel* oneAtomCS
	  = (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];

	G4double interCrsSec = util.DataLogInterpolation(LowEdgeEnergy, (*(*oneAtomCS)[0]), (*(*oneAtomCS)[1]))*barn;
	  //DataLogInterpolation(LowEdgeEnergy, tableIndex, theCrossSectionTable)*barn;
	SIGMA += theAtomNumDensityVector[k]*interCrsSec;
      }       
      
      Value = SIGMA<=0.0 ? BigPath : 1./SIGMA ;
      
      ptrVector->PutValue( i , Value ) ;
    }
    
    theMeanFreePathTable->insertAt( J , ptrVector ) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4VParticleChange* G4LowEnergyPhotoElectric::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){

  // incoming particle initialization
  aParticleChange.Initialize(aTrack);

  G4Material* aMaterial = aTrack.GetMaterial();

  const G4DynamicParticle* aDynamicPhoton = aTrack.GetDynamicParticle();
  const G4double PhotonEnergy = aDynamicPhoton->GetKineticEnergy();
  if(PhotonEnergy <= LowestEnergyLimit){
    
    aParticleChange.SetStatusChange(fStopAndKill);
    aParticleChange.SetEnergyChange(0.);
    aParticleChange.SetLocalEnergyDeposit(PhotonEnergy);
    
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  const G4ParticleMomentum PhotonDirection = aDynamicPhoton->GetMomentumDirection();
   
  // select randomly one element constituing the material.
  G4Element* anElement = SelectRandomAtom(aDynamicPhoton, aMaterial);

  // PAY ATTENTION TO THE MEANING OF THIS NUMBER
  G4int AtomNum = (G4int) anElement->GetZ();

  // Select the subshell WARNING!!!!
  G4int subShellIndex = SelectRandomShell(AtomNum, PhotonEnergy);

  G4FirstLevel* theBindEnVec = (*theBindingEnergyTable)[AtomNum-1];
  G4int thePrimaryShell = (G4int) (*(*theBindEnVec)[0])[subShellIndex];
  G4double BindingEn = (*(*theBindEnVec)[1])[subShellIndex];
  
  if(thePrimShVec.length() != 0){
    
    thePrimShVec.clear();
  }

  thePrimShVec.insert(thePrimaryShell);

  // Create lists of pointers to DynamicParticles (photons and electrons)
  G4ParticleVector photvec;
  G4int photInd = 0; 
  G4ParticleVector elecvec;
  G4int elecInd = 0; 

  // primary outcoming electron
  G4double ElecKineEnergy = (PhotonEnergy - BindingEn)*MeV;

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
      G4int AtomInd = ZNumVecFluor->index(AtomNum);
      oneAtomTable* oneAtomFluorTrans = (*theFluorTransitionTable)[AtomInd];
      
      while(ThereAreShells == TRUE){
	
	// Select the second transition from another subshell
	// fluorPar[0] = SubShell 
	// fluorPar[1] = Sec SubShell (if there is), 
	// fluorPar[2] = Transition Probability
	// the same for augerPar
	
	G4double fluorPar[3] = {0};
	ThereAreShells = SelectRandomTransition(thePrimaryShell, 
						fluorPar, 
						oneAtomFluorTrans);

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

	  if(fluorPar[2] >= CutForLowEnergySecondaryPhotons){

	    theEnergyDeposit -= fluorPar[2]*MeV;
	    
	    newPart = new G4DynamicParticle (G4Gamma::Gamma(), 
					     newPartDirection, 
					     fluorPar[2]) ;
	    photvec.append(newPart);
	  }
	}
	else{
	  
	  
	  G4int k = 0;
	  while(thePrimaryShell != (*(*theBindEnVec)[0])[k]) k++;
	  
	  G4double lastTransEnergy = (*(*theBindEnVec)[1])[k];
	  thePrimaryShell = (G4int) fluorPar[0];

	  if(fluorPar[2] >= CutForLowEnergySecondaryPhotons){

	    theEnergyDeposit -= lastTransEnergy*MeV;
	
	    newPart = new G4DynamicParticle (G4Gamma::Gamma(), 
					     newPartDirection, 
					     lastTransEnergy) ;
	    photvec.append(newPart);

	  }
	}

	  thePrimShVec.insert(thePrimaryShell);
      }
    } //END OF THE CHECK ON ATOMIC NUMBER
    
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

    if(theEnergyDeposit < 0){
      
      theEnergyDeposit = 0;
    }
    
  } // END OF CUTS
  
  else{
    
    ElecKineEnergy = 0. ;
    aParticleChange.SetNumberOfSecondaries(0) ;
  }

  // Kill the incident photon 
  aParticleChange.SetMomentumChange( 0., 0., 0. );
  aParticleChange.SetEnergyChange( 0. );

  if(theEnergyDeposit < 0){
    theEnergyDeposit = 0;
  }

  aParticleChange.SetLocalEnergyDeposit(theEnergyDeposit) ;  
  aParticleChange.SetStatusChange( fStopAndKill ) ; 
#ifdef G4VERBOSE
  if(verboseLevel > 15){
    G4cout<<"LE PhotoElectric PostStepDoIt"<<endl;
  }
#endif
  // Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );

}

G4int G4LowEnergyPhotoElectric::SelectRandomShell(const G4int AtomIndex, const G4double IncEnergy){
  
  G4double partialSum = 0;
  G4double totalSum = ComputeCrossSection(AtomIndex,IncEnergy);

  G4double rval = totalSum*G4UniformRand();
  const oneAtomTable* oneAtomCS 
    = (*allAtomShellCrossSec)[ZNumVec->index(AtomIndex)];

  for(G4int ind = 0; ind < oneAtomCS->entries(); ind++){

    G4double crossSec;
    G4Data* EnergyVector = (*(*oneAtomCS)[ind])[0];
    G4Data* CrossSecVector = (*(*oneAtomCS)[ind])[1];
    if(IncEnergy < (*EnergyVector)[0]){ //First element is the shell number

      crossSec = 0;
    }

    else{

      crossSec = util.DataLogInterpolation(IncEnergy, (*EnergyVector), (*CrossSecVector));

    }
    
    partialSum += crossSec;

    if(rval <= partialSum) return ind;
  }
  
  G4Exception("LEPhotoElectric: Cannot select a shell");
  return 0;
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

  G4double rval = G4UniformRand()/MeanFreePath;
 
  for ( G4int i=0 ; i < NumberOfElements ; i++ ){ 

    G4double crossSection;
    if (GammaEnergy <  LowestEnergyLimit)

      crossSection = 0. ;

    else {

      if (GammaEnergy > HighestEnergyLimit) GammaEnergy = 0.99*HighestEnergyLimit ;

      G4int AtomIndex = (G4int) (*theElementVector)(i)->GetZ();
      const G4FirstLevel* oneAtomCS
	= (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];

      crossSection =  util.DataLogInterpolation(GammaEnergy, (*(*oneAtomCS)[0]), (*(*oneAtomCS)[1]))*barn;
	//DataLogInterpolation(GammaEnergy, tableIndex, theCrossSectionTable)*barn;
    }

    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;

    if (rval <= PartialSumSigma) return ((*theElementVector)(i));

  }

  //  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
  // << "' has no elements" << endl;
  return (*theElementVector)(0);
}

G4bool G4LowEnergyPhotoElectric::SelectRandomTransition(G4int thePrimShell, 
							G4double* TransParam,
							const oneAtomTable* TransitionTable){
  
  G4int SubShellCol = 0, ProbCol = 1, EnergyCol = 2;
  //transitionTable means  for one atom not for one shell

  // too check when the subshell are finished
  G4bool ColIsFull = TRUE;
  G4int ShellNum = 0;
  G4double TotalSum = 0; 
  G4int maxNumOfShells = TransitionTable->entries()-1;

  if(thePrimShell <= (*(*(*TransitionTable)[maxNumOfShells])[0])[0]){

    while(thePrimShell != (*(*(*TransitionTable)[ShellNum])[0])[0]){
  
      if(ShellNum == maxNumOfShells){
	break;
      }
      
      ShellNum++;
    }
    
    //    if(ShellNum <= maxNumOfShells) {
      
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
	  TransParam[1] = (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
	  TransParam[2] = (*(*(*TransitionTable)[ShellNum])[EnergyCol])[TransProb];
	  break;
	}
	
	TransProb++;
      }

    //if(TransProb == (*(*TransitionTable)[ShellNum])[ProbCol]->length()-1) {

    //ColIsFull = FALSE;
    //}
      //}
      //else{
      
      // ColIsFull = FALSE;
      //}
  }
  else{

    ColIsFull = FALSE;
  }

  return ColIsFull;
}







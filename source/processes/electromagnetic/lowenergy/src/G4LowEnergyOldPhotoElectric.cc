//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4LowEnergyOldPhotoElectric.cc,v 1.2 2001-09-23 20:07:07 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      ------------ G4LowEnergyPhotoelctric: low energy modifications --------
//                   by Alessandra Forti, October 1998
// **************************************************************
//   10.04.2000 VL
// - Correcting Fluorescence transition probabilities in order to take into account 
//   non-radiative transitions. No Auger electron simulated yet: energy is locally deposited.
// 17.02.2000 Veronique Lefebure
// - bugs corrected in fluorescence simulation: 
//   . when final use of binding energy: no photon was ever created
//   . no Fluorescence was simulated when the photo-electron energy
//     was below production threshold.
//
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Added EnergySampling method A. Forti
// Modified PostStepDoIt to insert sampling with EPDL97 data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
// 07-09-99, if no e- emitted: edep=photon energy, mma
// 24.04.01 V.Ivanchenko remove RogueWave 
//                                  
// --------------------------------------------------------------

// This Class Header
#include "G4LowEnergyOldPhotoElectric.hh"

// Collaborating Class Headers
#include "G4EnergyLossTables.hh"
#include "G4Electron.hh"

typedef G4std::vector<G4DynamicParticle*> G4ParticleVector;

//    ..

// constructor
 
G4LowEnergyOldPhotoElectric::G4LowEnergyOldPhotoElectric(const G4String& processName)
  : G4VDiscreteProcess(processName),             // initialization
  lowestEnergyLimit (250*eV),
  highestEnergyLimit(100*GeV),
  NumbBinTable(200),
  CutForLowEnergySecondaryPhotons(0.),
  theCrossSectionTable(0),
  theMeanFreePathTable(0),
  allAtomShellCrossSec(0),
  theFluorTransitionTable(0),
  theBindingEnergyTable(0),
  ZNumVec(0),
  ZNumVecFluor(0),
  MeanFreePath(0.)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< G4endl;
     G4cout << "lowestEnergy: " << lowestEnergyLimit/keV << "keV ";
     G4cout << "highestEnergy: " << highestEnergyLimit/MeV << "MeV " << G4endl;
   }
}

//    ..
 
// destructor
 
G4LowEnergyOldPhotoElectric::~G4LowEnergyOldPhotoElectric()
{
   if (theCrossSectionTable) {
      delete theCrossSectionTable;
   }

   if (theBindingEnergyTable) {
     //      theBindingEnergyTable->clearAndDestroy();
      theBindingEnergyTable->clear();
      delete theBindingEnergyTable;
   }

   if (theMeanFreePathTable) {
     //      theMeanFreePathTable->clearAndDestroy();
      theMeanFreePathTable->clear();
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
     ZNumVecFluor->erase(ZNumVecFluor->begin(),ZNumVecFluor->end());
     delete ZNumVecFluor;
   }

}
 
//    ..
 
void G4LowEnergyOldPhotoElectric::SetCutForLowEnSecPhotons(G4double cut){

  CutForLowEnergySecondaryPhotons = cut;
}

//    ..

void G4LowEnergyOldPhotoElectric::BuildPhysicsTable(const G4ParticleDefinition& PhotonType)

// Build microscopic cross section table and mean free path table
{

  BuildZVec();

  BuildCrossSectionTable();

  BuildShellCrossSectionTable();

  BuildMeanFreePathTable();
  
  BuildBindingEnergyTable();

  BuildFluorTransitionTable();
   
}

//    ..

// CONSTRUCT THE CROSS SECTION TABLE FOR THE ELEMENTS MAPPED IN ZNUMVEC USING EPDL DATA
void G4LowEnergyOldPhotoElectric::BuildCrossSectionTable(){

  if (theCrossSectionTable) {
    
    delete theCrossSectionTable; 
  }

  theCrossSectionTable = new G4SecondLevel();
  G4int dataNum = 2;
 
  for(size_t TableInd = 0; TableInd < ZNumVec->size(); TableInd++){

    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];

    G4FirstLevel* oneAtomCS = util.BuildFirstLevelTables(AtomInd, dataNum, "phot/pe-cs-");
     
    //     theCrossSectionTable->insert(oneAtomCS);
     theCrossSectionTable->push_back(oneAtomCS);
   
  }//end for on atoms
}

//    ..

// CONSTRUCT THE SUBSHELL CS TABLE FOR THE ELEMENTS MAPPED IN ZNUMVEC USING EPDL DATA
void G4LowEnergyOldPhotoElectric::BuildShellCrossSectionTable(){

   if (allAtomShellCrossSec) {

    delete allAtomShellCrossSec;
   }

   allAtomShellCrossSec = new allAtomTable();
   G4int dataNum = 2;
 
   for(size_t TableInd = 0; TableInd < ZNumVec->size(); TableInd++){

     G4int AtomInd = (G4int) (*ZNumVec)[TableInd];

     oneAtomTable* oneAtomShellCS = util.BuildSecondLevelTables(AtomInd, dataNum, "phot/pe-ss-cs-");
     
     //     allAtomShellCrossSec->insert(oneAtomShellCS);
     allAtomShellCrossSec->push_back(oneAtomShellCS);
   
   }//end for on atoms
}

//    ..

// CONSTRUCT THE BE TABLE FOR THE ELEMENTS MAPPED IN ZNUMVEC USING EADL DATA
void G4LowEnergyOldPhotoElectric::BuildBindingEnergyTable(){

  if (theBindingEnergyTable) {

    delete theBindingEnergyTable;
  }

  G4int dataNum = 2;
  theBindingEnergyTable = util.BuildSecondLevelTables(0,dataNum,"fluor/binding");
}

//    ..

// CONSTRUCT THE FTP TABLE FOR THE ELEMENTS MAPPED IN ZNUMVEC USING EADL DATA
void G4LowEnergyOldPhotoElectric::BuildFluorTransitionTable(){

  if (theFluorTransitionTable) {
    
    delete theFluorTransitionTable;
  }
  
  theFluorTransitionTable = new allAtomTable();
  ZNumVecFluor = new G4DataVector(*ZNumVec);
  G4int dataNum = 3;
  
  for(size_t TableInd = 0; TableInd < ZNumVec->size(); TableInd++){
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    if(AtomInd > 5){
      
      oneAtomTable* oneAtomShellFL = util.BuildSecondLevelTables(AtomInd, dataNum, "fluor/fl-tr-pr-");
      //      theFluorTransitionTable->insert(oneAtomShellFL);
      theFluorTransitionTable->push_back(oneAtomShellFL);
    }
    else{
      ZNumVecFluor->remove(AtomInd);
    }
  }//end for on atoms
}

//    ..

//
// vector mapping the elements of the material table 
// 
void G4LowEnergyOldPhotoElectric::BuildZVec(){

  const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

  if(ZNumVec){
    ZNumVec->clear();
    delete ZNumVec;
  }

  ZNumVec = new G4DataVector(); 
  for (G4int J=0 ; J < numOfMaterials; J++){ 
 
    const G4Material* material= (*theMaterialTable)[J];        
    const G4ElementVector* theElementVector = material->GetElementVector();
    const G4int NumberOfElements = material->GetNumberOfElements() ;

    for (G4int iel=0; iel<NumberOfElements; iel++ ){

      G4double Zel = (*theElementVector)(iel)->GetZ();

      if(ZNumVec->contains(Zel) == FALSE){
	ZNumVec->push_back(Zel);
      } else{
	continue;
      }
    }
  }
}

//    ..

// Compute total cross section from subshell integrated cross section: needed for 
// selection of the first subshell ionized.

G4double G4LowEnergyOldPhotoElectric::ComputeCrossSection(const G4double AtomIndex,
						    const G4double IncEnergy){
  // calculates the microscopic cross section from subshell cross sections
  //(it is called for elements , AtomicNumber = Z )
 
  G4double TotalCrossSection(0.);

  const oneAtomTable* oneAtomCS
    = (*allAtomShellCrossSec)[ZNumVec->index(AtomIndex)];

  for(size_t ind = 0; ind < oneAtomCS->size(); ind++){

    G4double crossSec = 0;
    G4DataVector* EnergyVector = (*(*oneAtomCS)[ind])[0];
    G4DataVector* CrossSecVector = (*(*oneAtomCS)[ind])[1];

    if(IncEnergy < (*EnergyVector)[1]){ // First element is the shell number

      crossSec = 0;
    }

    else{

      crossSec = util.DataLogInterpolation(IncEnergy, (*EnergyVector), (*CrossSecVector))*barn;

    }

    TotalCrossSection += crossSec;
  }

  return TotalCrossSection ;
}

//    ..

void G4LowEnergyOldPhotoElectric::BuildMeanFreePathTable(){

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
    // WARNING: Lower limit of total cross sections in the data is the binding energy 
    // of the relative subshell. MeanFreePath table require a common lowest limit.
    // This lowestEnergyLimit is at the moment fixed at 250 ev.
 
    ptrVector = new  G4PhysicsLogVector(lowestEnergyLimit, highestEnergyLimit, NumbBinTable);
    
    material = (*theMaterialTable)(J);
    const G4ElementVector* theElementVector = material->GetElementVector();
    const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();   
    
    for ( G4int i = 0 ; i < NumbBinTable ; i++ ){ 
      //For each energy
      
      LowEdgeEnergy = ptrVector->GetLowEdgeEnergy(i);
      
      G4double SIGMA = 0;
      
      for ( size_t k=0 ; k < material->GetNumberOfElements() ; k++ ){ 
	// For each element            
	G4int AtomIndex = (G4int) (*theElementVector)(k)->GetZ();
	const G4FirstLevel* oneAtomCS
	  = (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];

	G4double interCrsSec = util.DataLogInterpolation(LowEdgeEnergy, (*(*oneAtomCS)[0]), (*(*oneAtomCS)[1]))*barn;

	SIGMA += theAtomNumDensityVector[k]*interCrsSec;
      }       
      
      Value = SIGMA > DBL_MIN ? 1./SIGMA : DBL_MAX ;
      
      ptrVector->PutValue( i , Value ) ;
    }
    
    theMeanFreePathTable->insertAt( J , ptrVector ) ;
  }
}

//    .. 

G4VParticleChange* G4LowEnergyOldPhotoElectric::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){

  // Fluorescence (as reported by stepanek):
  // J. Stepanek " A program to determine the radiation spectra due to a single atomic 
  // subshell ionisation by a particle or due to deexcitation or decay of radionuclides", 
  // Comp. Phys. Comm. 1206 pp 1-1-9 (1997)
  // 
  // incoming particle initialization
  aParticleChange.Initialize(aTrack);

  G4Material* aMaterial = aTrack.GetMaterial();

  const G4DynamicParticle* aDynamicPhoton = aTrack.GetDynamicParticle();
  const G4double PhotonEnergy = aDynamicPhoton->GetKineticEnergy();
  if(PhotonEnergy <= lowestEnergyLimit){
    
    aParticleChange.SetStatusChange(fStopAndKill);
    aParticleChange.SetEnergyChange(0.);
    aParticleChange.SetLocalEnergyDeposit(PhotonEnergy);
    
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  const G4ParticleMomentum PhotonDirection = aDynamicPhoton->GetMomentumDirection();
   
  // select randomly one element constituing the material.
  G4Element* anElement = SelectRandomAtom(aDynamicPhoton, aMaterial);

  // PAY ATTENTION TO THE MEANING OF THIS NUMBER!!! SelectRandomAtom requires to use AtomNum
  // the BindingEnergyTable requires AtomNum-1
  G4int AtomNum = (G4int) anElement->GetZ();

  // First Ionised subshell is chosen basing on subshell integrated cross section EPDL97
  // using the partial sum method.
  // Select the subshell WARNING!!!!: it returns the subshell index in the table.

  G4int subShellIndex = SelectRandomShell(AtomNum, PhotonEnergy);

  G4FirstLevel* theBindEnVec = (*theBindingEnergyTable)[AtomNum-1];
  G4int thePrimaryShell = (G4int) (*(*theBindEnVec)[0])[subShellIndex];
  G4double BindingEn = ((*(*theBindEnVec)[1])[subShellIndex])*MeV;
  
  if(thePrimShVec.size() != 0){
    
    thePrimShVec.clear();
  }

  thePrimShVec.push_back(thePrimaryShell);

  // Create lists of pointers to DynamicParticles (photons and electrons)
  G4ParticleVector photvec;
  // G4int photInd = 0; 
  G4ParticleVector elecvec;
  // G4int elecInd = 0; 

  // primary outcoming electron
  G4double ElecKineEnergy = (PhotonEnergy - BindingEn);

  G4double theEnergyDeposit = BindingEn;

  if (G4EnergyLossTables::GetRange(G4Electron::Electron(),ElecKineEnergy,aMaterial)
      >= G4std::min(G4Electron::GetCuts(), aStep.GetPostStepPoint()->GetSafety())){

    // the electron is created in the direction of the incident photon ...  
    
    G4DynamicParticle* aElectron = new G4DynamicParticle (G4Electron::Electron(), 
							  PhotonDirection, ElecKineEnergy) ;
    elecvec.push_back(aElectron);
  } // END OF CUTS
  
  else{
    theEnergyDeposit += ElecKineEnergy;    
  }

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
	/////newPartDirection.rotateUz(PhotonDirection);
	
	if(ThereAreShells != FALSE){
	  
	  thePrimaryShell = (G4int) fluorPar[0];

	  if(fluorPar[2]*MeV >= CutForLowEnergySecondaryPhotons){

	    theEnergyDeposit -= fluorPar[2]*MeV;
	    
	    newPart = new G4DynamicParticle (G4Gamma::Gamma(), 
					     newPartDirection, 
					     fluorPar[2]*MeV);
	    //	    photvec.append(newPart);
	    photvec.push_back(newPart);
	  }
	}
	else{
	  	  
	  /////Energy deposition vl
	  ////=================NEW================vl
	  
	  /*
	  G4int k = 0;
	  while(thePrimaryShell != (*(*theBindEnVec)[0])[k]) k++;
	  
	  G4double lastTransEnergy = ((*(*theBindEnVec)[1])[k])*MeV;
	  thePrimaryShell = (G4int) fluorPar[0];

	  if(lastTransEnergy >= CutForLowEnergySecondaryPhotons){

	    theEnergyDeposit -= lastTransEnergy;
	
	    newPart = new G4DynamicParticle (G4Gamma::Gamma(), 
					     newPartDirection, 
					     lastTransEnergy) ;
	    photvec.push_back(newPart);
           
	  }
	  thePrimShVec.insert(thePrimaryShell);
	  */
	}

      }
    } //END OF THE CHECK ON ATOMIC NUMBER
    
    G4int numOfElec = elecvec.size();
    G4int numOfPhot = photvec.size();
    G4int numOfDau  = numOfElec + numOfPhot;

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
    

  // Kill the incident photon 
  aParticleChange.SetMomentumChange( 0., 0., 0. );
  aParticleChange.SetEnergyChange( 0. );

  if (theEnergyDeposit < 0) theEnergyDeposit = 0;
  aParticleChange.SetLocalEnergyDeposit(theEnergyDeposit);  
  aParticleChange.SetStatusChange( fStopAndKill ); 

  // Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );

}

//    ..

G4int G4LowEnergyOldPhotoElectric::SelectRandomShell(const G4int AtomIndex,
                                                  const G4double IncEnergy)
{  
  G4double partialSum = 0;
  G4double totalSum = ComputeCrossSection(AtomIndex,IncEnergy);

  G4double rval = totalSum*G4UniformRand();
  const oneAtomTable* oneAtomCS 
    = (*allAtomShellCrossSec)[ZNumVec->index(AtomIndex)];

  for(size_t ind = 0; ind < oneAtomCS->size(); ind++){

    G4double crossSec;
    G4DataVector* EnergyVector = (*(*oneAtomCS)[ind])[0];
    G4DataVector* CrossSecVector = (*(*oneAtomCS)[ind])[1];
    if(IncEnergy < (*EnergyVector)[0]){ //First element is the shell number

      crossSec = 0;
    }

    else{

      crossSec = util.DataLogInterpolation(IncEnergy, (*EnergyVector), (*CrossSecVector))*barn;

    }
    
    partialSum += crossSec;

    if(rval <= partialSum) return ind;
  }
  
  G4Exception("LEPhotoElectric: Cannot select a shell");
  return 0;
}

//    ..

G4Element*
G4LowEnergyOldPhotoElectric::SelectRandomAtom(const G4DynamicParticle* aDynamicPhoton,
                                           G4Material* aMaterial)
{
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
    if (GammaEnergy <  lowestEnergyLimit)

      crossSection = 0. ;

    else {

      if (GammaEnergy > highestEnergyLimit) GammaEnergy = 0.99*highestEnergyLimit ;

      G4int AtomIndex = (G4int) (*theElementVector)(i)->GetZ();
      const G4FirstLevel* oneAtomCS
	= (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];

      crossSection =  util.DataLogInterpolation(GammaEnergy, (*(*oneAtomCS)[0]), (*(*oneAtomCS)[1]))*barn;
    }

    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;

    if (rval <= PartialSumSigma) return ((*theElementVector)(i));

  }
  return (*theElementVector)(0);
}
//    ..

//
// Select a random transition with the transition probabilities and the partial sum 
// method using EADL data (A. Forti)
//

G4bool G4LowEnergyOldPhotoElectric::SelectRandomTransition(G4int thePrimShell, 
							G4double* TransParam,
							const oneAtomTable* TransitionTable){
  
  G4int SubShellCol = 0, ProbCol = 1, EnergyCol = 2;
  // transitionTable contains all the transition probabilities of one atom:
  // loop on subshell is inside the method.

  // when the last subshell is reached CollIsFull becomes FALSE.
  G4bool ColIsFull = FALSE;
  G4int ShellNum = 0;
  //  G4double TotalSum = 0; 
  G4int maxNumOfShells = TransitionTable->size()-1;

  if(thePrimShell <= 0) {
     G4cerr<<"*** Unvalid Primary shell: "<<thePrimShell<<G4endl;
     return FALSE;
  }   
  if(thePrimShell <= (*(*(*TransitionTable)[maxNumOfShells])[0])[0]){

    while(thePrimShell != (*(*(*TransitionTable)[ShellNum])[0])[0]){
  
      if(ShellNum == maxNumOfShells){
	break;
      }
      
      ShellNum++;
    }
    
    // TransProb is the index of the loop and of the table of transition. it starts from 1 
    // because the first element of the data table is the primary shell id number and not a 
    // transition probability: it must not be added to TotalSum. 

      G4int TransProb = 1;
      
     // Include non-radiative transitions (vl):
     //// for(TransProb = 1; TransProb < (*(*TransitionTable)[ShellNum])[ProbCol]->length(); TransProb++){ 
     ////	TotalSum += (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
     //// }
     ////G4double PartialProb = G4UniformRand()*TotalSum;
     ////
     G4double PartialProb = G4UniformRand();
    
    
    //vl.


      G4double PartSum = 0;
      
      TransProb = 1; 
      G4int trSize = (*(*TransitionTable)[ShellNum])[ProbCol]->size();
      while(TransProb < trSize){
	
	PartSum += (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
	
	if(PartialProb <= PartSum){
	  
	  TransParam[0] = (*(*(*TransitionTable)[ShellNum])[SubShellCol])[TransProb];
	  TransParam[1] = (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
	  TransParam[2] = (*(*(*TransitionTable)[ShellNum])[EnergyCol])[TransProb];
 	  ColIsFull = TRUE;
	  break;
	}
	
	TransProb++;
      }
  }
  else{

    ColIsFull = FALSE;
  }

  return ColIsFull;
}

//    ..






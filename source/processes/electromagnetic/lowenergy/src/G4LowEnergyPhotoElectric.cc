// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyPhotoElectric.cc,v 1.1 1999-03-02 17:17:54 aforti Exp $
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
#include <rw/tvordvec.h>

// constructor
 
G4LowEnergyPhotoElectric::G4LowEnergyPhotoElectric(const G4String& processName)
  : G4VDiscreteProcess(processName),             // initialization
    LowestEnergyLimit (100*eV),
    HighestEnergyLimit(100*GeV),
    NumbBinTable(100)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
     G4cout << "LowestEnergy: " << LowestEnergyLimit/keV << "keV ";
     G4cout << "HighestEnergy: " << HighestEnergyLimit/MeV << "MeV " << endl;
   }
   theCrossSectionTable = 0;
   theBindingEnergyTable = 0;
   theMeanFreePathTable = 0;
   theFluorTransitionTable = 0;
   theAugerTransitionTable = 0;
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

   if (theFluorTransitionTable) {
      theFluorTransitionTable->clearAndDestroy();
      delete theFluorTransitionTable;
   }

   if (theAugerTransitionTable) {
      theAugerTransitionTable->clearAndDestroy();
      delete theAugerTransitionTable;
   }
}
 
 
// methods.............................................................................
 
void G4LowEnergyPhotoElectric::BuildPhysicsTable(const G4ParticleDefinition& PhotonType)

// Build microscopic cross section table and mean free path table
{

   // Build microscopic cross section tables for the Photo Electric Effect
   BuildCrossSectionTable();

   // Build Binding Energy Table
   BuildBindingEnergyTable();

   // Build mean free path table for the Compton Scattering process
   BuildMeanFreePathTable();


}

G4VParticleChange* G4LowEnergyPhotoElectric::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){

   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();

   const G4DynamicParticle* aDynamicPhoton = aTrack.GetDynamicParticle();
   G4double PhotonEnergy = aDynamicPhoton->GetKineticEnergy();
   G4ParticleMomentum PhotonDirection = aDynamicPhoton->GetMomentumDirection();
   
   // select randomly one element constituing the material.
   G4Element* anElement = SelectRandomAtom(aDynamicPhoton, aMaterial);

   G4int AtomNum = anElement->GetZ();

   // load the transition probability table for the element
   BuildFluorTransitionTable(AtomNum);

   BuildAugerTransitionTable(AtomNum);

   // Photo electron
   G4PhysicsVector* theBindEnVec = (*theBindingEnergyTable)(AtomNum);
    
   G4int g = 0;
   while(g < theBindEnVec->GetVectorLength() && PhotonEnergy < (*theBindEnVec)(g)) g++;

   G4int thePrimaryShell = theBindEnVec->GetLowEdgeEnergy(g);

   G4double ElecKineEnergy = PhotonEnergy - (*theBindEnVec)(g);

   if (G4EnergyLossTables::GetRange(G4Electron::Electron(),ElecKineEnergy,aMaterial)
        >= min(G4Electron::GetCuts(), aStep.GetPostStepPoint()->GetSafety()) ){
     
     // the electron is created in the direction of the incident photon ...  
     G4DynamicParticle* aElectron = new G4DynamicParticle (G4Electron::Electron(), 
							  PhotonDirection, ElecKineEnergy) ;
     aParticleChange.SetNumberOfSecondaries(1) ;
     aParticleChange.AddSecondary( aElectron ) ; 
   }
   
   else{

     ElecKineEnergy = 0. ;
     aParticleChange.SetNumberOfSecondaries(0) ;
   }

   G4int fluorPar[2] = {0};
   SelectRandomTransition(thePrimaryShell, fluorPar, theFluorTransitionTable);

   G4int augerPar[2] = {0};
   SelectRandomTransition(thePrimaryShell, augerPar, theAugerTransitionTable);

   G4double theFluorProb = (*(*(*theFluorTransitionTable)[fluorPar[0]])[1])[fluorPar[1]];
   G4double theAugerProb = (*(*(*theAugerTransitionTable)[augerPar[0]])[1])[augerPar[1]];
   G4double theSh, theSecSh, theEn;

   G4double newcosTh = 1-2*G4UniformRand();
   G4double newsinTh = sqrt(1-newcosTh*newcosTh);
   G4double newPhi = twopi*G4UniformRand();

   G4double dirx, diry, dirz;
   dirz = newcosTh;
   diry = newsinTh*cos(newPhi);
   dirx = newsinTh*sin(newPhi);
   G4ThreeVector newPartDirection(dirx, diry,dirz);
   G4DynamicParticle* newPart;

   if((theFluorProb + theAugerProb)*G4UniformRand() < theFluorProb){

     theSh = (*(*(*theFluorTransitionTable)[fluorPar[0]])[0])[fluorPar[1]];
     theSecSh = 0;
     theEn = (*(*(*theFluorTransitionTable)[fluorPar[0]])[2])[fluorPar[1]];
     newPart = new G4DynamicParticle (G4Gamma::Gamma(), newPartDirection, theEn) ;
   }
   else{
     
     theSh = (*(*(*theAugerTransitionTable)[augerPar[0]])[0])[augerPar[1]];
     theSecSh = (*(*(*theAugerTransitionTable)[augerPar[0]])[2])[augerPar[1]];
     theEn = (*(*(*theAugerTransitionTable)[augerPar[0]])[3])[augerPar[1]];
     newPart = new G4DynamicParticle (G4Electron::Electron(), newPartDirection, theEn) ;
   }


   //controllare se il setnumberofsecondaries  si puo' cambiare
   aParticleChange.SetNumberOfSecondaries(1) ;
   aParticleChange.AddSecondary( newPart ) ; 

   // Kill the incident photon 
   aParticleChange.SetMomentumChange( 0., 0., 0. ) ;
   aParticleChange.SetEnergyChange( 0. ) ;
   aParticleChange.SetLocalEnergyDeposit( PhotonEnergy - ElecKineEnergy ) ;  
   aParticleChange.SetStatusChange( fStopAndKill ) ; 

   // WARNING!!!! FLUORESCENCE EFFECT MUST BE ADDED

   //  Reset NbOfInteractionLengthLeft and return aParticleChange
   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

void G4LowEnergyPhotoElectric::BuildCrossSectionTable(){
 
  if (theCrossSectionTable) {
    
    theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable; 
  }

  G4int par[4] = {73, 0, 0, 0}; G4String name("epdl97");
  G4Epdl97File File(name,par);
  G4EpdlTables table(File);
  table.FillDataTable();
  theCrossSectionTable = new G4PhysicsTable(*(table.GetFstDataTable())) ;

}

void G4LowEnergyPhotoElectric::BuildBindingEnergyTable(){
 
  if (theBindingEnergyTable) {
    
    theBindingEnergyTable->clearAndDestroy(); delete theBindingEnergyTable; 
  }

  G4int par[4] = {91, 913, 0, 0}; G4String name("eadl.asc");
  G4Epdl97File File(name,par);
  G4EpdlTables table(File);
  table.FillDataTable();
  theBindingEnergyTable = new G4PhysicsTable(*(table.GetFstDataTable())) ;

}

void G4LowEnergyPhotoElectric::BuildFluorTransitionTable(G4int atomNum){
 
  if (theFluorTransitionTable) {
  
    theFluorTransitionTable->clearAndDestroy(); delete theFluorTransitionTable;
  }

  G4String name("eadl.asc");
  G4int par[4] = {92, 931, 91, 0}; 
  G4Epdl89File File(name,par);

  G4EpdlTables table(File);
  table.FillTheTable(atomNum);
  theFluorTransitionTable = table.GetGlobalList() ;
}

void G4LowEnergyPhotoElectric::BuildAugerTransitionTable(G4int  atomNum){
 
  if (theAugerTransitionTable) {
    
    theAugerTransitionTable->clearAndDestroy(); delete theAugerTransitionTable; 
  }

  G4int par[4] = {92, 932, 91, 0}; G4String name("eadl.asc");
  G4Epdl97File File(name,par);
  G4EpdlTables table(File);
  table.FillTheTable(atomNum);
  theAugerTransitionTable = table.GetGlobalList() ;
}

void G4LowEnergyPhotoElectric::BuildMeanFreePathTable(){

  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable; }

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
	
	G4double interCrsSec = DataLogInterpolation(LowEdgeEnergy, (*theElementVector)(k)->GetZ(), theCrossSectionTable)*barn;
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
  G4double rval = G4UniformRand()/MeanFreePath;
 
  for ( G4int i=0 ; i < NumberOfElements ; i++ ){ 

    G4double crossSection;
    if (GammaEnergy <  LowestEnergyLimit)
      crossSection = 0. ;
    else {
      if (GammaEnergy > HighestEnergyLimit) GammaEnergy = 0.99*HighestEnergyLimit ;
      crossSection = DataLogInterpolation(GammaEnergy, (*theElementVector)(i)->GetZ(), theCrossSectionTable)*barn;
    }

    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;
    if (rval <= PartialSumSigma) return ((*theElementVector)(i));
  }
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << endl;
  return 0;
}


void G4LowEnergyPhotoElectric::SelectRandomTransition(G4int thePrimShell, 
						      G4int* TransParam,
						      RWTPtrSlist< RWTPtrSlist<G4DataVector> >* TransitionTable){

  G4int SubShellCol = 0, ProbCol = 1, EnergyCol = 2;
  G4int ShellNum = 0;
  G4double TotalSum = 0; 
  
  while(thePrimShell != (*(*(*TransitionTable)[ShellNum])[0])[0]) ShellNum++;

  //TransProb start from 1 because the first element of the list is the primary shall id number
  for(G4int TransProb = 1; TransProb < (*(*TransitionTable)[ShellNum])[SubShellCol]->length(); TransProb++){ 
    
    TotalSum += (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
  }
  
  G4double PartialProb = G4UniformRand()*TotalSum;

  G4double PartSum = 0;
  
  TransProb = 1; 
  while(TransProb < (*(*TransitionTable)[ShellNum])[SubShellCol]->length()){
    
    PartSum += (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];

    if(PartialProb < PartSum){
      
      TransParam[0] = ShellNum;
      TransParam[1] = TransProb;
      break;
    }
    TransProb++;
  }
}







// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyRayleigh.cc,v 1.11 1999-07-30 15:18:46 aforti Exp $
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
//      ------------ G4LowEnergyRayleigh physics process --------
//                   by Michel Maire, April 1996
// **************************************************************
// 28-05-96, DoIt() small change in ElecDirection, by M.Maire
// 10-06-96, simplification in ComputeMicroscopicCrossSection(), by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 13-09-96, small changes in DoIt for better efficiency. Thanks to P.Urban
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 05-03-97, new Physics scheme, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 07-04-98, remove 'tracking cut' of the scattered gamma, MMa
// 04-06-98, in DoIt, secondary production condition: range>min(threshold,safety)
// --------------------------------------------------------------

// This Class Header
#include "G4LowEnergyRayleigh.hh"

// Collaborating Class Headers
#include "G4EnergyLossTables.hh"
#include "G4Electron.hh"

// constructor
 
G4LowEnergyRayleigh::G4LowEnergyRayleigh(const G4String& processName)
  : G4VDiscreteProcess(processName),
    theCrossSectionTable(0),
    theMeanFreePathTable(0),
    theFormFactorTable(0),
    ZNumVec(0),
    LowestEnergyLimit (250*eV),              // initialization
    HighestEnergyLimit(100*GeV),
    NumbBinTable(200)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
     G4cout << "LowestEnergy: " << LowestEnergyLimit/keV << "keV ";
     G4cout << "HighestEnergy: " << HighestEnergyLimit/TeV << "TeV " << endl;
   }
}
 
// destructor
 
G4LowEnergyRayleigh::~G4LowEnergyRayleigh()
{
   if (theCrossSectionTable) {

      delete theCrossSectionTable;
   }

   if(theFormFactorTable){

     delete theFormFactorTable;
   }

   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }

   if(ZNumVec){
     
     ZNumVec->clear();
     delete ZNumVec;
   }
}
 
 
// methods.............................................................................
 
void G4LowEnergyRayleigh::BuildPhysicsTable(const G4ParticleDefinition& GammaType){

  BuildZVec();

  // Build microscopic cross section tables for the Rayleigh process
  BuildCrossSectionTable();
  
  // Build mean free path table for the Rayleigh Scattering process
  BuildMeanFreePathTable();
  
  // build the scattering function table
  BuildFormFactorTable();
}

void G4LowEnergyRayleigh::BuildCrossSectionTable(){
 
  if (theCrossSectionTable) {
    
    delete theCrossSectionTable; 
  }

  theCrossSectionTable = new G4SecondLevel();
  G4int dataNum = 2;
  
  for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){
    
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    
    G4FirstLevel* oneAtomCS = util.BuildFirstLevelTables(AtomInd, dataNum, "rayl/re-cs-");
    
    theCrossSectionTable->insert(oneAtomCS);
    
  }//end for on atoms
}

void G4LowEnergyRayleigh::BuildFormFactorTable(){
 
  if (theFormFactorTable) {
    
    delete theFormFactorTable; 
  }

  theFormFactorTable = new G4SecondLevel();
  G4int dataNum = 2;
  
  for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){
    
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    
    G4FirstLevel* oneAtomFF = util.BuildFirstLevelTables(AtomInd, dataNum, "rayl/re-ff-");
    
    theFormFactorTable->insert(oneAtomFF);
    
  }//end for on atoms
}

void G4LowEnergyRayleigh::BuildZVec(){

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

G4VParticleChange* G4LowEnergyRayleigh::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){

//
// The scattered gamma energy is sampled according to Form Factors 
// multiplied by the Rayleigh distribution with a pure rejection technique.  
// EGS4 W.R. Nelson et al. The EGS4 Code System. SLAC-Report-265 , December 1985 
// Expression of the angular distribution as Rayleigh distribution and Form factors 
// is taken from D. E. Cullen "A simple model of photon transport" Nucl. Instr. Meth. 
// Phys. Res. B 101 (1995). Method of sampling with form factors is different.
// Reference to the article is from J. Stepanek New Photon, Positron
// and Electron Interaction Data for GEANT in Energy Range from 1 eV to 10
// TeV (draft). 


  aParticleChange.Initialize(aTrack);
  // Dynamic particle quantities  
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  G4double GammaEnergy0 = aDynamicGamma->GetKineticEnergy();
  
  if(GammaEnergy0 <= LowestEnergyLimit){
    
    aParticleChange.SetStatusChange(fStopAndKill);
    aParticleChange.SetEnergyChange(0.);
    aParticleChange.SetLocalEnergyDeposit(GammaEnergy0);
    
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  G4double E0_m = GammaEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum GammaDirection0 = aDynamicGamma->GetMomentumDirection();
  
  // Select randomly one element
  G4Material* aMaterial = aTrack.GetMaterial();
  const G4int numOfElem = aMaterial->GetNumberOfElements();
  G4Element* theElement = SelectRandomAtom(aDynamicGamma, aMaterial);
  
  // sample the energy rate of the scattered gamma 

  G4double wlGamma = h_Planck*c_light/GammaEnergy0;
  G4int elementZ = (G4int) theElement->GetZ();
  G4double tableIndex = elementZ - 1;

  G4double Theta, DataFormFactor;
  G4double cosTheta, greject;
  G4double Theta_Half, x, SinThHalf, RandomFormFactor;
  G4double sinTheta;
  do{
    
    Theta_Half = G4UniformRand()*pi/2;
    SinThHalf = sin(Theta_Half);
    x = SinThHalf/wlGamma;
    
    const G4FirstLevel* oneAtomFF
      = (*theFormFactorTable)[ZNumVec->index(elementZ)];

    DataFormFactor = util.DataLogInterpolation(x, (*(*oneAtomFF)[0]), 
					       (*(*oneAtomFF)[1]))/cm;
    RandomFormFactor = G4UniformRand()*elementZ;

    Theta = Theta_Half*2;
    cosTheta = cos(Theta);
    sinTheta = sin(Theta);
    
    greject = cosTheta*cosTheta*DataFormFactor;

  }while( greject < RandomFormFactor);

  
  // scattered gamma angles. ( Z - axis along the parent gamma)
  G4double Phi = twopi * G4UniformRand() ;
  G4double dirx = sinTheta*cos(Phi) , diry = sinTheta*sin(Phi) , dirz = cosTheta ;
  
  // update G4VParticleChange for the scattered gamma 
  G4ThreeVector GammaDirection1(dirx, diry, dirz);

  GammaDirection1.rotateUz(GammaDirection0);
  aParticleChange.SetEnergyChange(GammaEnergy0);
  aParticleChange.SetMomentumChange(GammaDirection1);
  
  aParticleChange.SetNumberOfSecondaries(0);
#ifdef G4VERBOSE
  if(verboseLevel > 15){
    G4cout<<"LE Rayleigh PostStepDoIt"<<endl;
  }
#endif

  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

void G4LowEnergyRayleigh::BuildMeanFreePathTable(){

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

	G4double AtomIndex = (*theElementVector)(k)->GetZ();

	const G4FirstLevel* oneAtomCS
	  = (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];
	
	G4double interCrsSec = util.DataLogInterpolation(LowEdgeEnergy, 
							 (*(*oneAtomCS)[0]), 
							 (*(*oneAtomCS)[1]))*barn;

	SIGMA += theAtomNumDensityVector[k]*interCrsSec;

      }
      
      Value = SIGMA<=0.0 ? BigPath : 1./SIGMA ;
      ptrVector->PutValue( i , Value ) ;
    }
    
    theMeanFreePathTable->insertAt( J , ptrVector ) ;
  }
}


G4Element* G4LowEnergyRayleigh::SelectRandomAtom(const G4DynamicParticle* aDynamicGamma,
                                               G4Material* aMaterial) {

  // select randomly 1 element within the material 
  G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  if (NumberOfElements == 1) return (*theElementVector)(0);

  const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();

  G4double PartialSumSigma = 0.;
  G4double rval = G4UniformRand()/MeanFreePath;

  for ( G4int i=0 ; i < NumberOfElements ; i++ ){ 

    G4double crossSection;
    if (GammaEnergy <  LowestEnergyLimit)
      crossSection = 0. ;
    else {
      if (GammaEnergy > HighestEnergyLimit) GammaEnergy = 0.99*HighestEnergyLimit ;

      G4double AtomIndex = (*theElementVector)(i)->GetZ();

      const G4FirstLevel* oneAtomCS
	= (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];
      
      crossSection =  util.DataLogInterpolation(GammaEnergy, 
						(*(*oneAtomCS)[0]), 
						(*(*oneAtomCS)[1]))*barn;
    }
    
    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;
    if(rval <= PartialSumSigma) return ((*theElementVector)(i));
  }
  //  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
  //   << "' has no elements" << endl;
  return (*theElementVector)(0);
}









// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyCompton.cc,v 1.7 1999-06-06 16:26:16 aforti Exp $
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
//      ------------ G4LowEnergyCompton physics process --------
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
#include "G4LowEnergyCompton.hh"

// C Headers

// C++ Headers

// Collaborating Class Headers
#include "G4EnergyLossTables.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Epdl97File.hh"
#include "G4EpdlTables.hh"
#include "G4PhysicsFreeVector.hh" 
#include "G4Step.hh"
#include "Randomize.hh" 

// RW Headers
#include <rw/tvordvec.h>

// constructor
 
G4LowEnergyCompton::G4LowEnergyCompton(const G4String& processName)
  : G4VDiscreteProcess(processName),
    LowestEnergyLimit (100*eV),              // initialization
    HighestEnergyLimit(100*GeV),
    NumbBinTable(100)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
     G4cout << "LowestEnergy: " << LowestEnergyLimit/keV << "keV ";
     G4cout << "HighestEnergy: " << HighestEnergyLimit/TeV << "TeV " << endl;
   }

   theCrossSectionTable = 0;
   theMeanFreePathTable = 0; 
   theScatteringFunctionTable = 0;
}
 
// destructor
 
G4LowEnergyCompton::~G4LowEnergyCompton()
{
   if (theCrossSectionTable) {
      theCrossSectionTable->clearAndDestroy();
      delete theCrossSectionTable;
   }

   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }

   if (theScatteringFunctionTable) {
      theScatteringFunctionTable->clearAndDestroy();
      delete theScatteringFunctionTable;
   }
}
 
 
// methods.............................................................................

// to change with other functions like in G4eIonization
 
void G4LowEnergyCompton::BuildPhysicsTable(const G4ParticleDefinition& GammaType){

  // Build microscopic cross section table and mean free path table
  BuildCrossSectionTable();

  // Build mean free path table for the Compton Scattering process
  BuildMeanFreePathTable();

  // build the scattering function table
  BuildScatteringFunctionTable();
}

G4VParticleChange* G4LowEnergyCompton::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){

//
// The scattered gamma energy is sampled according to Klein - Nishina formula.
// The random number techniques of Butcher & Messel are used 
// (Nuc Phys 20(1960),15).
// GEANT4 internal units
//
  if(getenv("GENERAL")) aParticleChange.Initialize(aTrack);
 
  // Dynamic particle quantities  
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  G4double GammaEnergy0 = aDynamicGamma->GetKineticEnergy();
  G4double E0_m = GammaEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum GammaDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element 
  G4Material* aMaterial = aTrack.GetMaterial();
  const G4int numOfElem = aMaterial->GetNumberOfElements();
  
  G4Element* theElement = SelectRandomAtom(aDynamicGamma, aMaterial);
  G4int elementZ = (G4int) theElement->GetZ();
  G4double epsilon, epsilonsq, onecost, sint2, greject ;

  G4double epsilon0 = 1./(1. + 2*E0_m) , epsilon0sq = epsilon0*epsilon0;
  G4double alpha1   = - log(epsilon0)  , alpha2 = 0.5*(1.- epsilon0sq);
  G4double ScatteringFunction, x;
  G4double wlGamma = h_Planck*c_light/GammaEnergy0;
  
  // sample the energy rate of the scattered gamma 
  do{
    
    if ( alpha1/(alpha1+alpha2) > G4UniformRand()){
      
      epsilon   = exp(-alpha1*G4UniformRand());  // pow(epsilon0,G4UniformRand())
      epsilonsq = epsilon*epsilon; 
    }
    else{
      
      epsilonsq = epsilon0sq + (1.- epsilon0sq)*G4UniformRand();
      epsilon   = sqrt(epsilonsq);
    }
    
    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    
    x = sqrt(onecost/2)/wlGamma;
    
    ScatteringFunction = DataLogInterpolation(x, elementZ - 1, theScatteringFunctionTable)/cm;
    greject = (1. - epsilon*sint2/(1.+ epsilonsq))*ScatteringFunction;
    
  }  while(2*greject < elementZ*G4UniformRand());
  
  G4double cosTeta = 1. - onecost , sinTeta = sqrt (sint2);
  G4double Phi     = twopi * G4UniformRand() ;
  G4double dirx = sinTeta*cos(Phi) , diry = sinTeta*sin(Phi) , dirz = cosTeta ;

  //
  // update G4VParticleChange for the scattered gamma 
  //
  
  G4ThreeVector GammaDirection1 ( dirx,diry,dirz );
  GammaDirection1.rotateUz(GammaDirection0);
  aParticleChange.SetMomentumChange( GammaDirection1 ) ;
  G4double GammaEnergy1 = epsilon*GammaEnergy0;
  if (GammaEnergy1 > 0.)
    {
      aParticleChange.SetEnergyChange( GammaEnergy1 ) ;
    }
  else
    {    
      aParticleChange.SetEnergyChange(0.) ;
      aParticleChange.SetStatusChange(fStopAndKill);
      
    }
  
  //
  // kinematic of the scattered electron
  //
  
  G4double ElecKineEnergy = GammaEnergy0 - GammaEnergy1 ;
  
  if (G4EnergyLossTables::GetRange(G4Electron::Electron(), ElecKineEnergy, aMaterial)
      >= min(G4Electron::GetCuts(), aStep.GetPostStepPoint()->GetSafety())){

    G4double ElecMomentum = sqrt(ElecKineEnergy*(ElecKineEnergy+2.*electron_mass_c2));
    G4ThreeVector ElecDirection((GammaEnergy0*GammaDirection0 - 
				 GammaEnergy1*GammaDirection1)*(1./ElecMomentum));
    
    // create G4DynamicParticle object for the electron.  
    G4DynamicParticle* aElectron= new G4DynamicParticle (G4Electron::Electron(),
							  ElecDirection, ElecKineEnergy) ;
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary( aElectron );
    aParticleChange.SetLocalEnergyDeposit (0.); 
  }
  else{
    
    aParticleChange.SetNumberOfSecondaries(0);
    aParticleChange.SetLocalEnergyDeposit (ElecKineEnergy);
   }
  
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

void G4LowEnergyCompton::BuildCrossSectionTable(){
 
  if (theCrossSectionTable) {
    
    theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable; 
  }

  G4int par[4] = {72, 0, 0, 0}; G4String name("epdl97");
  G4Epdl97File File(name,par);
  G4EpdlTables table(File);
  table.FillDataTable();
  theCrossSectionTable = table.GetFstDataTable();
}

void G4LowEnergyCompton::BuildScatteringFunctionTable(){
 
  if (theScatteringFunctionTable) {
    
    theScatteringFunctionTable->clearAndDestroy(); delete theScatteringFunctionTable; 
  }

  G4int par[4] = {93, 942, 0, 0}; G4String name("epdl97");
  G4Epdl97File File(name,par);
  G4EpdlTables table(File);
  table.FillDataTable();
  theScatteringFunctionTable = table.GetFstDataTable();
}

void G4LowEnergyCompton::BuildMeanFreePathTable(){

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
	
	G4double interCrsSec = DataLogInterpolation(LowEdgeEnergy, (*theElementVector)(k)->GetZ(), theCrossSectionTable)*barn;
	SIGMA += theAtomNumDensityVector[k]*interCrsSec;
      }       
      
      Value = SIGMA<=0.0 ? BigPath : 1./SIGMA ;

      ptrVector->PutValue( i , Value ) ;

    }
    
    theMeanFreePathTable->insertAt( J , ptrVector );
  }
}


G4Element* G4LowEnergyCompton::SelectRandomAtom(const G4DynamicParticle* aDynamicGamma,
                                               G4Material* aMaterial){
  // select randomly 1 element within the material 
  G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  if (NumberOfElements == 1) return (*theElementVector)(0);

  const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  //GetMeanFreePath
  G4double PartialSumSigma = 0.;

  G4double rval = 0;
  rval = G4UniformRand()/MeanFreePath;

  for ( G4int i=0 ; i < NumberOfElements ; i++ ){ 

    G4double crossSection;
    if (GammaEnergy <  LowestEnergyLimit)
      crossSection = 0. ;
    else {
      if (GammaEnergy > HighestEnergyLimit) GammaEnergy = 0.99*HighestEnergyLimit ;
      crossSection = DataLogInterpolation(GammaEnergy, (*theElementVector)(i)->GetZ(), theCrossSectionTable)*barn;
    }

    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;
    if(rval <= PartialSumSigma) return ((*theElementVector)(i));
  }

  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
	 << "' has no elements" << endl;
  return (*theElementVector)(0);
}













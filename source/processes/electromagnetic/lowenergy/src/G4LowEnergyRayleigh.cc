// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyRayleigh.cc,v 1.6 1999-06-05 14:04:07 aforti Exp $
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

#include "G4LowEnergyRayleigh.hh"
#include "G4EnergyLossTables.hh"
#include "G4Epdl97File.hh"
#include "G4EpdlTables.hh"
#include "CLHEP/String/Strings.h"
#include <rw/tvordvec.h>

// constructor
 
G4LowEnergyRayleigh::G4LowEnergyRayleigh(const G4String& processName)
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
   theFormFactorTable = 0;
}
 
// destructor
 
G4LowEnergyRayleigh::~G4LowEnergyRayleigh()
{
   if (theCrossSectionTable) {
      theCrossSectionTable->clearAndDestroy();
      delete theCrossSectionTable;
   }

   if(theFormFactorTable){
     theFormFactorTable->clearAndDestroy();
     delete theFormFactorTable;
   }

   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }
}
 
 
// methods.............................................................................
 
void G4LowEnergyRayleigh::BuildPhysicsTable(const G4ParticleDefinition& GammaType){

   // Build microscopic cross section tables for the Rayleigh process
   BuildCrossSectionTable();

   // Build mean free path table for the Rayleigh Scattering process
   BuildMeanFreePathTable();

   // build the scattering function table
   BuildFormFactorTable();
   
}

G4VParticleChange* G4LowEnergyRayleigh::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){

//
// The scattered gamma energy is sampled according to Form Factors and 
// then accepted or rejected based on Rayleigh distribution.
// The random number techniques of Butcher & Messel are used 
// (Nuc Phys 20(1960),15). GEANT4 internal units
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
    
    DataFormFactor = DataLogInterpolation(x, tableIndex, theFormFactorTable)/cm;
    
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
  
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

void G4LowEnergyRayleigh::BuildCrossSectionTable(){
 
  if (theCrossSectionTable) {
    
    theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable; 
  }

  G4int par[4] = {71, 0, 0, 0}; G4String name("epdl97");
  G4Epdl97File File(name,par);
  G4EpdlTables table(File);
  table.FillDataTable();
  theCrossSectionTable = table.GetFstDataTable();
}

void G4LowEnergyRayleigh::BuildFormFactorTable(){
 
  if (theFormFactorTable) {
    
    theFormFactorTable->clearAndDestroy(); delete theFormFactorTable; 
  }

  G4int par[4] = {93, 941, 0, 0}; G4String name("epdl97");
  G4Epdl97File File(name,par);
  G4EpdlTables table(File);
  table.FillDataTable();
  theFormFactorTable = table.GetFstDataTable();
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

	G4double tableIndex = (*theElementVector)(k)->GetZ() - 1;
	G4double interCrsSec = DataLogInterpolation(LowEdgeEnergy,tableIndex, theCrossSectionTable)*barn;

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

      G4double tableIndex = (*theElementVector)(i)->GetZ() - 1;
      crossSection = DataLogInterpolation(GammaEnergy, tableIndex, theCrossSectionTable)*barn;
      
    }
    
    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;
    if(rval <= PartialSumSigma) return ((*theElementVector)(i));
  }
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements" << endl;
  return (*theElementVector)(0);
}
























// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyPhotoElectric.cc,v 1.1 1999-01-08 14:16:35 gunter Exp $
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

#include "G4LowEnergyPhotoElectric.hh"
#include "G4EnergyLossTables.hh"
#include "G4EpdlData.hh"
#include "CLHEP/String/Strings.h"
#include <rw/tvordvec.h>
// constructor
 
G4LowEnergyPhotoElectric::G4LowEnergyPhotoElectric(const G4String& processName)
  : G4PhotoElectricEffect(processName),             // initialization
    LowestEnergyLimit (1*keV),
    HighestEnergyLimit(50*MeV),
    MeanFreePath(0),
    NumbBinTable(100)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
     G4cout << "LowestEnergy: " << LowestEnergyLimit/keV << "keV ";
     G4cout << "HighestEnergy: " << HighestEnergyLimit/MeV << "MeV " << endl;
   }
   theCrossSectionTable = 0;
   theMeanFreePathTable = 0;
   theBindingEnergyTable = 0;
}
 
// destructor
 
G4LowEnergyPhotoElectric::~G4LowEnergyPhotoElectric()
{
   if (theCrossSectionTable) {
      theCrossSectionTable->clearAndDestroy();
      delete theCrossSectionTable;
   }

   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }

   if (theBindingEnergyTable) {
      theBindingEnergyTable->clearAndDestroy();
      delete theBindingEnergyTable;
   }
}
 
 
// methods.............................................................................
 
void G4LowEnergyPhotoElectric::BuildPhysicsTable(const G4ParticleDefinition& PhotonType)

// Build microscopic cross section table and mean free path table
{
   G4double Value;
   G4PhysicsFreeVector* ptrVector;

// Build microscopic cross section tables for the Photo Electric Effect

   if (theCrossSectionTable) {
           theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable; }

   G4double param[4] = {73, 0, 0, 0};
   HepString file("epdl.asc");
   G4EpdlData c(file, param);
   c.FillDataTable(1., 1.);
   theCrossSectionTable = new G4PhysicsTable(*(c.GetFstDataTable())) ;


   // Build the binding energy table
   if(theBindingEnergyTable){
     theBindingEnergyTable->clearAndDestroy(); delete theBindingEnergyTable;
   }
   
   G4double param1[4] = {91, 913, 0, 0};
   HepString file1("eadl.asc");
   G4EpdlData c1(file1, param1);
   c1.FillDataTable( 1, 1.);
   theBindingEnergyTable = new G4PhysicsTable(*(c1.GetFstDataTable())) ;

   // Build mean free path table for the Compton Scattering process
   BuildMeanFreePathTable();

}

G4VParticleChange* G4LowEnergyPhotoElectric::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep)
//
// Generate an electron resulting of a photo electric effect.
// The incident photon disappear.
// GEANT4 internal units
//
 
{
   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();

   const G4DynamicParticle* aDynamicPhoton = aTrack.GetDynamicParticle();
   G4double PhotonEnergy = aDynamicPhoton->GetKineticEnergy();
   G4ParticleMomentum PhotonDirection = aDynamicPhoton->GetMomentumDirection();
   
   // select randomly one element constituing the material.
   G4double Z = SelectRandomAtom(aDynamicPhoton, aMaterial)->GetZ();

   // Photo electron
   G4double ElecKineEnergy = PhotonEnergy - ComputeKBindingEnergy(Z) ;

   if (ElecKineEnergy < 0.) ElecKineEnergy = PhotonEnergy - ComputeL1BindingEnergy(Z); 
   if (ElecKineEnergy < 0.) ElecKineEnergy = PhotonEnergy - ComputeL2BindingEnergy(Z);
   if (ElecKineEnergy < 0.) ElecKineEnergy = PhotonEnergy ;

   if (G4EnergyLossTables::GetRange(G4Electron::Electron(),ElecKineEnergy,aMaterial)
        >= min(G4Electron::GetCuts(), aStep.GetPostStepPoint()->GetSafety()) ){
     
     // the electron is created in the direction of the incident photon ...  
     G4DynamicParticle* aElectron= new G4DynamicParticle (G4Electron::Electron(), PhotonDirection, ElecKineEnergy) ;
     aParticleChange.SetNumberOfSecondaries(1) ;
     aParticleChange.AddSecondary( aElectron ) ; 
   }
   
   else{
     ElecKineEnergy = 0. ;
     aParticleChange.SetNumberOfSecondaries(0) ;
   }

   // Kill the incident photon 
   aParticleChange.SetMomentumChange( 0., 0., 0. ) ;
   aParticleChange.SetEnergyChange( 0. ) ;
   aParticleChange.SetLocalEnergyDeposit( PhotonEnergy - ElecKineEnergy ) ;  
   aParticleChange.SetStatusChange( fStopAndKill ) ; 

   // WARNING!!!! FLUORESCENCE EFFECT MUST BE ADDED

   //  Reset NbOfInteractionLengthLeft and return aParticleChange
   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

G4double G4LowEnergyPhotoElectric::DataLogInterpolation(G4double Argument, G4double AtomicNumber, G4PhysicsTable* Table){

  G4PhysicsVector* theVec = (*Table)(AtomicNumber);
  G4int theLoc = FindBinLocation(Argument, theVec); 

  G4double val1 = (*theVec)(theLoc), val2 = (*theVec)(theLoc+1);
  G4double arg1 = theVec->GetLowEdgeEnergy(theLoc), arg2 = theVec->GetLowEdgeEnergy(theLoc+1);

  G4double theVal = (log10(val1)*log10(arg2/Argument)
		     +log10(val2)*log10(Argument/arg1))/log10(arg2/arg1);
  
  theVal = exp(theVal);

  return theVal;
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


G4int G4LowEnergyPhotoElectric::FindBinLocation(G4double arg, G4PhysicsVector* vec){

  G4int numberOfBin = vec->GetVectorLength();
  G4int lowerBound = 0;
  G4int upperBound = numberOfBin-1;

  do {
    G4int midBin = (lowerBound + upperBound)/2;
    if( arg < vec->GetLowEdgeEnergy(midBin) )
       upperBound = midBin-1;
    else

       lowerBound = midBin+1;
  } while (lowerBound <= upperBound); 

  return upperBound;
}









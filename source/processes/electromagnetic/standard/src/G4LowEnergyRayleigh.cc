// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyRayleigh.cc,v 1.1 1999-01-08 14:16:37 gunter Exp $
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
#include "G4EpdlData.hh"
#include "CLHEP/String/Strings.h"
//#include "globals.hh"
#include <rw/tvordvec.h>

// constructor
 
G4LowEnergyRayleigh::G4LowEnergyRayleigh(const G4String& processName)
  : G4ComptonScattering(processName),
    LowestEnergyLimit (1*keV),              // initialization
    HighestEnergyLimit(100*GeV),
    MeanFreePath(0),
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

// Build microscopic cross section table and mean free path table
   G4double LowEdgeEnergy, Value;
   G4PhysicsFreeVector* ptrVector;

   // Build microscopic cross section tables for the Rayleigh Scattering process
    if (theCrossSectionTable) {
            theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable;}
   
   G4double param[4] = {71, 0, 0, 0};
   HepString file("epdl.asc");
   G4EpdlData c(file, param);
   c.FillDataTable(1., 1.);
   theCrossSectionTable = new G4PhysicsTable(*(c.GetFstDataTable())) ;

   // build the scattering function table
   if (theFormFactorTable) {
              theFormFactorTable->clearAndDestroy(); delete theFormFactorTable; 
   }
   
   G4double param1[4] = {71, 941, 0, 0};
 
   HepString file1("epdl.asc");
   G4EpdlData c1(file1, param1);
   c1.FillDataTable(1., 1.);
   theFormFactorTable = new G4PhysicsTable(*(c1.GetFstDataTable())) ;
   cout<<" *************** Rayleigh Form Factor Table **************** "<<endl;
   for(G4int U =0; U<theFormFactorTable->length(); U++){
     cout<<"ELEMENT "<<U+1<<":"<<endl;
     
     for(G4int f =0; f<(*theCrossSectionTable)(U)->GetVectorLength() ; f++){
       
       cout<<"datvec["<<f<<"]: "<<(*(*theCrossSectionTable)(U))(f)<<"    binvec["<<f<<"]: "<<(*theCrossSectionTable)(U)->GetLowEdgeEnergy(f)<<endl;
     }
   }


   // Build mean free path table for the Rayleigh Scattering process
   BuildMeanFreePathTable();
}

G4VParticleChange* G4LowEnergyRayleigh::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){
//
// The scattered gamma energy is sampled according to Form Factors and then accepted or rejected
// based on Rayleigh distribution.
// The random number techniques of Butcher & Messel are used (Nuc Phys 20(1960),15).
// GEANT4 internal units
//
// Note : Effects due to binding of atomic electrons are negliged.

  aParticleChange.Initialize(aTrack);
  
  // Dynamic particle quantities  
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  G4double GammaEnergy0 = aDynamicGamma->GetKineticEnergy();
  G4double E0_m = GammaEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum GammaDirection0 = aDynamicGamma->GetMomentumDirection();
  
  //  const G4ElementVector* theElementVector = aMaterial->GetElementVector() ;

  // Select randomly one element
  G4Material* aMaterial = aTrack.GetMaterial();
  const G4int numOfElem = aMaterial->GetNumberOfElements();
  G4Element* theElement = 0;

  G4int count =0;
  while(!theElement){
    
    theElement = SelectRandomAtom(aDynamicGamma, aMaterial);
    count++;
    
    if(count>10){
      
      cout<<"Ci ho provato "<<count<<" volte e non funziona!"<<endl;
      exit(1);
    }
  }

  // sample the energy rate of the scattered gamma 
  G4double Val, ffPar;
  G4double cosTh2, sinTheta, greject;


  do {
      
    // For the Rayleigh effect dataVector (form factors) and binVector inside a PhysicsVector
    // are inverted to sample the cosTheta with the data distribution: Val IS the form factor

    G4double hc = h_Planck*c_light;
    G4double wlGamma = hc/GammaEnergy0;

    Val = G4UniformRand()*theElement->GetZ();
    G4double sinThOnLambda = DataLogInterpolation(Val, theElement->GetZ(), theFormFactorTable)/cm;
    sinTheta = wlGamma*sinThOnLambda;
    cosTh2 = 1-sinTheta*sinTheta;
    greject = (1+cosTh2)/2;

  }while( (greject < G4UniformRand()));


   // scattered gamma angles. ( Z - axis along the parent gamma)
   G4double cosTheta = sqrt (cosTh2);
   G4double Phi     = twopi * G4UniformRand() ;
   G4double dirx = sinTheta*cos(Phi) , diry = sinTheta*sin(Phi) , dirz = cosTheta ;

   // update G4VParticleChange for the scattered gamma 
   G4ThreeVector GammaDirection1 (dirx, diry, dirz);
   GammaDirection1.rotateUz(GammaDirection0);
   aParticleChange.SetMomentumChange(GammaDirection1) ;
   if (GammaEnergy0 > 0.)
     {
       aParticleChange.SetEnergyChange(GammaEnergy0) ;
     }
   else
     {    
       aParticleChange.SetEnergyChange(0.) ;
       aParticleChange.SetStatusChange(fStopAndKill);
     }
   /*
   // kinematic of the scattered electron that in Rayleigh effect remain at rest
   if (G4EnergyLossTables::GetRange(G4Electron::Electron(), 0, aMaterial)
   >= min(G4Electron::GetCuts(), aStep.GetPostStepPoint()->GetSafety()) ){
     
     G4double ElecMomentum = 0;
     G4ThreeVector ElecDirection(0., 0., 0.);
     
     // create G4DynamicParticle object for the electron.  
     G4DynamicParticle* aElectron= new G4DynamicParticle (G4Electron::Electron(), ElecDirection, 0) ;
     aParticleChange.SetNumberOfSecondaries(1) ;
     aParticleChange.AddSecondary( aElectron ) ;
     aParticleChange.SetLocalEnergyDeposit (0.) ; 
   }	
   else{
     
     aParticleChange.SetNumberOfSecondaries(0) ;
     aParticleChange.SetLocalEnergyDeposit (0) ;
   }
   */
   //  Reset NbOfInteractionLengthLeft and return aParticleChange
   
   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}


G4double G4LowEnergyRayleigh::DataLogInterpolation(G4double Argument, G4double AtomicNumber, G4PhysicsTable* Table){

  G4PhysicsVector* theVec = (*Table)(AtomicNumber);
  G4int theLoc = FindBinLocation(Argument, theVec); 

  G4double val1 = (*theVec)(theLoc), val2 = (*theVec)(theLoc+1);
  G4double arg1 = theVec->GetLowEdgeEnergy(theLoc), arg2 = theVec->GetLowEdgeEnergy(theLoc+1);

  G4double theVal = (log10(val1)*log10(arg2/Argument)
		     +log10(val2)*log10(Argument/arg1))/log10(arg2/arg1);
  
  theVal = exp(theVal);

  return theVal;
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
	G4double interCrsSec = DataLogInterpolation(LowEdgeEnergy, (*theElementVector)(k)->GetZ(), theCrossSectionTable)*barn;
	SIGMA += theAtomNumDensityVector[k]*interCrsSec;
      }       
      
      Value = SIGMA<=0.0 ? BigPath : 1./SIGMA ;
      
      ptrVector->PutValue( i , Value ) ;
    }
    
    theMeanFreePathTable->insertAt( J , ptrVector ) ;
  }
}


G4Element* G4LowEnergyRayleigh::SelectRandomAtom(const G4DynamicParticle* aDynamicGamma,
                                               G4Material* aMaterial) const {
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
      crossSection = DataLogInterpolation(GammaEnergy, (*theElementVector)(i)->GetZ(), theCrossSectionTable)*barn;
      
    }
    
    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;
    if(rval <= PartialSumSigma) return ((*theElementVector)(i));
  }
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << endl;
  return 0;
}

G4int G4LowEnergyRayleigh::FindBinLocation(G4double arg, G4PhysicsVector* vec){

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


























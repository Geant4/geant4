// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyCompton.cc,v 1.1 1999-01-08 14:16:32 gunter Exp $
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

#include "G4LowEnergyCompton.hh"
#include "G4EnergyLossTables.hh"
#include "G4EpdlData.hh"
#include "G4PhysicsLogVector.hh" 
//#include "CLHEP/String/Strings.h"
//#include "globals.hh"
#include <rw/tvordvec.h>

// constructor
 
G4LowEnergyCompton::G4LowEnergyCompton(const G4String& processName)
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

   G4double Value;
   G4PhysicsFreeVector* ptrVector;

   // Build microscopic cross section tables for the Compton Scattering process

   if (theCrossSectionTable) {
     theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable;}
   
   G4double param[4] = {72, 0, 0, 0};
   HepString file("epdl.asc");
   G4EpdlData c(file, param);
   c.FillDataTable(1., 1.);
   theCrossSectionTable = new G4PhysicsTable(*(c.GetFstDataTable())) ;
   /* cout<<"************ the CST ****************"<<endl;
   for(G4int U =0; U<100; U++){
     cout<<"ELEMENT "<<U+1<<":"<<endl;
     
     for(G4int f =0; f<(*theCrossSectionTable)(U)->GetVectorLength() ; f++){

     cout<<"datvec["<<f<<"]: "<<(*(*theCrossSectionTable)(U))(f)<<"    binvec["<<f<<"]: "<<(*theCrossSectionTable)(U)->GetLowEdgeEnergy(f)<<endl;
     }
   }*/
   // build the scattering function table

   if (theScatteringFunctionTable) {
     theScatteringFunctionTable->clearAndDestroy(); delete theScatteringFunctionTable; 
   }

   G4double param1[4] = {72, 942, 0, 0};
   HepString file1("epdl.asc");
   G4EpdlData c1(file1, param1);
   c1.FillDataTable(1., 1.);
   
   theScatteringFunctionTable = new G4PhysicsTable(*(c1.GetFstDataTable())) ;
   
   // Build mean free path table for the Compton Scattering process
   BuildMeanFreePathTable();
}

G4VParticleChange* G4LowEnergyCompton::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep){

//
// The scattered gamma energy is sampled according to Klein - Nishina formula.
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

  // Select randomly one element 
  G4Material* aMaterial = aTrack.GetMaterial();
  const G4int numOfElem = aMaterial->GetNumberOfElements();
  G4Element* theElement = 0;
  G4int count = 0;

  while(!theElement){

    theElement = SelectRandomAtom(aDynamicGamma, aMaterial);
    count++;
    
    if(count>10){
      
      cout<<"Ci ho provato "<<count<<" volte e non funziona!"<<endl;
      exit(1);
    }
  }

  // sample the energy rate of the scattered gamma 
  G4double Val, sfPar;
  G4double ERate, ERate2, onecost, sint2; 
  
  G4double ERateMin  = 1./(1.+2*E0_m) , ERateMin2 = ERateMin*ERateMin;
  G4double NormF1 = - log(ERateMin),  NormF2 = 0.5*(1.- ERateMin2);
  
  do {
    if ( NormF1/(NormF1+NormF2) > G4UniformRand() ){ 

      ERate  = exp(-NormF1*G4UniformRand());     // pow(ERateMin,G4UniformRand())
      ERate2 = ERate*ERate; 
    }

    else {

      ERate2 = ERateMin2 + (1.- ERateMin2)*G4UniformRand();
      ERate  = sqrt(ERate2);
    }

    onecost = (1.- ERate)/(ERate*E0_m);
    sint2 = onecost*(2.-onecost);
    //    greject = 1. - ERate*sint2/(1.+ ERate2);

    G4double hc = h_Planck*c_light;
    G4double wlGamma = hc/GammaEnergy0;
    sfPar = sqrt(onecost/2)*wlGamma;
      
    Val = DataLogInterpolation(sfPar, theElement->GetZ(), theScatteringFunctionTable)/cm;
    
  }  while(Val >=  theElement->GetZ()*G4UniformRand());
  
   G4double cosTeta = 1. - onecost , sinTeta = sqrt (sint2);
   G4double Phi     = twopi * G4UniformRand() ;
   G4double dirx = sinTeta*cos(Phi) , diry = sinTeta*sin(Phi) , dirz = cosTeta ;

   //
   // update G4VParticleChange for the scattered gamma 
   //
   
   G4ThreeVector GammaDirection1 ( dirx,diry,dirz );
   GammaDirection1.rotateUz(GammaDirection0);
   aParticleChange.SetMomentumChange( GammaDirection1 ) ;
   G4double GammaEnergy1 = ERate*GammaEnergy0;
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
     G4ThreeVector ElecDirection((GammaEnergy0*GammaDirection0 - GammaEnergy1*GammaDirection1)*(1./ElecMomentum));
     
     // create G4DynamicParticle object for the electron.  
     G4DynamicParticle* aElectron= new G4DynamicParticle (G4Electron::Electron(), ElecDirection, ElecKineEnergy) ;
     aParticleChange.SetNumberOfSecondaries(1) ;
     aParticleChange.AddSecondary( aElectron ) ;
     aParticleChange.SetLocalEnergyDeposit (0.) ; 
   }
   else{
     
     aParticleChange.SetNumberOfSecondaries(0) ;
     aParticleChange.SetLocalEnergyDeposit (ElecKineEnergy) ;
   }

   //  Reset NbOfInteractionLengthLeft and return aParticleChange

   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

G4double G4LowEnergyCompton::DataLogInterpolation(G4double Argument, G4double AtomicNumber, G4PhysicsTable* Table){

  G4PhysicsVector* theVec = (*Table)(AtomicNumber);
  G4int theLoc = FindBinLocation(Argument, theVec); 

  G4double val1 = (*theVec)(theLoc), val2 = (*theVec)(theLoc+1);
  G4double arg1 = theVec->GetLowEdgeEnergy(theLoc), arg2 = theVec->GetLowEdgeEnergy(theLoc+1);

  G4double theVal = (log10(val1)*log10(arg2/Argument)
		     +log10(val2)*log10(Argument/arg1))/log10(arg2/arg1);
  
  theVal = exp(theVal);
  //  cout<<"ATOMICN: "<<AtomicNumber<<"  theLoc: "<<theLoc<<endl;
  //  cout<<"    ******************************   "<<endl;
  //  cout<<"Argument: "<<Argument<<"  arg1: "<<arg1<<"  arg2: "<<arg2<<endl;
  //  cout<<"theVal: "<<theVal<<"  val1: "<<val1<<"  val2: "<<val2<<endl;
  return theVal;
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
	// For each element            
	//	cout<<"BMFP: dens: "<<theAtomNumDensityVector[k]<<" DLInerp: "
	//          <<DataLogInterpolation(LowEdgeEnergy, (*theElementVector)(k)->GetZ(), theCrossSectionTable)<<endl;
	//	cout<<"**************************************"<<endl;
	G4double interCrsSec = DataLogInterpolation(LowEdgeEnergy, (*theElementVector)(k)->GetZ(), theCrossSectionTable)*barn;
	SIGMA += theAtomNumDensityVector[k]*interCrsSec;
	  
	//	cout<<"BMFP: sigma: "<<SIGMA<<endl;

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
     cout<<"MFP: "<<MeanFreePath<<"   crossSection: "<<crossSection<<endl;
    }

    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;
    cout<<"rval: "<<rval<<"   PSS: "<<PartialSumSigma<<endl;
    if(rval <= PartialSumSigma) return ((*theElementVector)(i));
  }
  cout<<"         *************************     "<<endl;
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
	 << "' has no elements, NULL pointer returned." << endl;
  return 0;
}

G4int G4LowEnergyCompton::FindBinLocation(G4double arg, G4PhysicsVector* vec){

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












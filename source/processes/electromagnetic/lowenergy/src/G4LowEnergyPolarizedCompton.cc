// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyPolarizedCompton.cc,v 1.2 2001-05-24 18:19:35 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group

// --------- G4LowEnergyPolarizedCompton class -----
//
//           by G.Depaola & F.Longo (21 may 2001)
// 21 May 2001 - MGP      Modified to inherit from G4VDiscreteProcess
//                        Applies same algorithm as LowEnergyCompton
//                        if the incoming photon is not polarised
//                        Temporary protection to avoid crash in the case 
//                        of polarisation || incident photon direction
//
// ************************************************************
//
// Corrections by Rui Curado da Silva (2000)
// New Implementation by G.Depaola & F.Longo      
//
// - sampling of phi
// - polarization of scattered photon
//
// --------------------------------------------------------------

#include "G4LowEnergyPolarizedCompton.hh"
#include "G4Electron.hh"
#include "G4EnergyLossTables.hh"
#include "G4Gamma.hh"
#include "G4SecondLevel.hh"
#include "G4PhysicsTable.hh"
#include "G4DataVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"

// constructor
 
G4LowEnergyPolarizedCompton::G4LowEnergyPolarizedCompton(const G4String& processName)
  : G4VDiscreteProcess(processName),
    theCrossSectionTable(0),
    theScatteringFunctionTable(0),
    theMeanFreePathTable(0),
    ZNumVec(0),
    lowestEnergyLimit (250*eV),              // initialization
    highestEnergyLimit(100*GeV),
    numbBinTable(200)
{
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
    G4cout << "LowestEnergy: " << lowestEnergyLimit/keV << "keV ";
    G4cout << "HighestEnergy: " << highestEnergyLimit/TeV << "TeV " << G4endl;
  }
}

// destructor
 
G4LowEnergyPolarizedCompton::~G4LowEnergyPolarizedCompton()
{
   if (theCrossSectionTable) {

      delete theCrossSectionTable;
   }

   if (theScatteringFunctionTable) {

      delete theScatteringFunctionTable;
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
 
 void G4LowEnergyPolarizedCompton::BuildPhysicsTable(const G4ParticleDefinition& GammaType){

  BuildZVec();

  // Build microscopic cross section table and mean free path table
  BuildCrossSectionTable();

  // Build mean free path table for the Compton Scattering process
  BuildMeanFreePathTable();

  // build the scattering function table
  BuildScatteringFunctionTable();

}
void G4LowEnergyPolarizedCompton::BuildCrossSectionTable(){

  // BUILD THE CS TABLE FOR THE ELEMENTS MAPPED IN ZNUMVEC
 
  if (theCrossSectionTable) {
    
    delete theCrossSectionTable; 
  }
  
  theCrossSectionTable = new G4SecondLevel();
  G4int dataNum = 2;
  
  for(size_t TableInd = 0; TableInd < ZNumVec->size(); TableInd++){
    
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    
    G4FirstLevel* oneAtomCS = util.BuildFirstLevelTables(AtomInd, dataNum, "comp/ce-cs-");
    
    //    theCrossSectionTable->insert(oneAtomCS);
    theCrossSectionTable->push_back(oneAtomCS);
    
  }//end for on atoms
}

void G4LowEnergyPolarizedCompton::BuildScatteringFunctionTable(){

// BUILD THE SF TABLE FOR THE ELEMENTS MAPPED IN ZNUMVEC

  if (theScatteringFunctionTable) {
    
    delete theScatteringFunctionTable; 
  }

  theScatteringFunctionTable = new G4SecondLevel();
  G4int dataNum = 2;
 
  for(size_t TableInd = 0; TableInd < ZNumVec->size(); TableInd++){

    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];

    G4FirstLevel* oneAtomSF = util.BuildFirstLevelTables(AtomInd, dataNum, "comp/ce-sf-");
     
    //     theScatteringFunctionTable->insert(oneAtomSF);
     theScatteringFunctionTable->push_back(oneAtomSF);
   
  }//end for on atoms
}


void G4LowEnergyPolarizedCompton::BuildZVec(){

// vector mapping the elements in the material table

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
      }  else{
	continue;
      }
    }
  }
}


void G4LowEnergyPolarizedCompton::BuildMeanFreePathTable(){

// used log-log interpolation instead of linear interpolation to build the MFP 
// as reported in the stepanek paper 

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
    ptrVector = new  G4PhysicsLogVector(lowestEnergyLimit, highestEnergyLimit, numbBinTable);
    
    material = (*theMaterialTable)(J);
    const G4ElementVector* theElementVector = material->GetElementVector();
    const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();   
    
    for ( G4int i = 0 ; i < numbBinTable ; i++ ){ 
      //For each energy
      
      LowEdgeEnergy = ptrVector->GetLowEdgeEnergy(i);
      
      const G4double BigPath= DBL_MAX;
      G4double SIGMA = 0 ;
      for ( size_t k=0 ; k < material->GetNumberOfElements() ; k++ ){ 

	G4int AtomIndex = (G4int) (*theElementVector)(k)->GetZ();
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
    
    theMeanFreePathTable->insertAt( J , ptrVector );
  }
}


G4Element* G4LowEnergyPolarizedCompton::SelectRandomAtom(const G4DynamicParticle* aDynamicGamma,
                                               G4Material* aMaterial){

// METHOD BELOW  FROM STANDARD E_M PROCESSES CODE MODIFIED TO USE 
// LIVERMORE DATA (using log-log interpolation as reported in stepanek paper)


  // select randomly 1 element within the material 
  G4double gammaEnergy = aDynamicGamma->GetKineticEnergy();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  if (NumberOfElements == 1) return (*theElementVector)(0);

  const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();
 
  G4double PartialSumSigma = 0.;

  G4double rval = 0;
  rval = G4UniformRand()/meanFreePath;

  for ( G4int i=0 ; i < NumberOfElements ; i++ ){ 

    G4double crossSection;
    if (gammaEnergy <  lowestEnergyLimit)
      crossSection = 0. ;
    else {
      if (gammaEnergy > highestEnergyLimit) gammaEnergy = 0.99*highestEnergyLimit ;
      
      G4int AtomIndex = (G4int) (*theElementVector)(i)->GetZ();
      const G4FirstLevel* oneAtomCS
	= (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];

      crossSection =  util.DataLogInterpolation(gammaEnergy, 
						(*(*oneAtomCS)[0]), 
						(*(*oneAtomCS)[1]))*barn;
    }

    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;
    if(rval <= PartialSumSigma) return ((*theElementVector)(i));
  }

  return (*theElementVector)(0);
}


G4VParticleChange* G4LowEnergyPolarizedCompton::PostStepDoIt(const G4Track& aTrack,
							     const G4Step&  aStep)
{
  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used (Nuc Phys 20(1960),15).
  // GEANT4 internal units
  //
  // Note : Effects due to binding of atomic electrons are negliged.

  aParticleChange.Initialize(aTrack);

  // Dynamic particle quantities  
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  G4double gammaEnergy0 = aDynamicGamma->GetKineticEnergy();
  G4ThreeVector gammaPolarization0 = aDynamicGamma->GetPolarization(); 
  G4double polarisation = gammaPolarization0.mag();

  // Check magnitude of polarisation vector
  G4bool isPolarised = false;
  if (polarisation > 0. && polarisation <= 1.) 
    {
      isPolarised = true;
    }

  // Temporary protection: for some unclear reason a polarisation parallel to the 
  // direction causes problems; in that case apply the regular LowEnergyCompton algorithm
  G4ThreeVector gammaDirection = aDynamicGamma->GetMomentumDirection();
  G4double angle = gammaPolarization0.angle(gammaDirection);
  if (angle == 0.)
    {
      isPolarised = false;
    }
  // End of temporary protection

  // Within energy limit?

  if(gammaEnergy0 <= lowestEnergyLimit)
{
  aParticleChange.SetStatusChange(fStopAndKill);
  aParticleChange.SetEnergyChange(0.);
  aParticleChange.SetLocalEnergyDeposit(gammaEnergy0);
  
  return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
}
  
  // Select randomly one element 
  G4Material* aMaterial = aTrack.GetMaterial();
  G4Element* theElement = SelectRandomAtom(aDynamicGamma, aMaterial);
  G4int elementZ = (G4int) theElement->GetZ();

  G4double E0_m = gammaEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum gammaDirection0 = aDynamicGamma->GetMomentumDirection();

  G4double epsilon, epsilonsq, onecost, sint2, greject ;

  G4double epsilon0 = 1./(1. + 2*E0_m);
  G4double epsilon0sq = epsilon0*epsilon0;
  G4double alpha1   = - log(epsilon0);
  G4double alpha2 = 0.5*(1.- epsilon0sq);
  G4double ScatteringFunction;
  G4double x;
  G4double wlGamma = h_Planck*c_light/gammaEnergy0;
  G4double gammaEnergy1;
  G4ThreeVector gammaDirection1;

  if (isPolarised)
    {
      do {
	if ( alpha1/(alpha1+alpha2) > G4UniformRand() )
	  { 
	    epsilon   = exp(-alpha1*G4UniformRand());  
	    epsilonsq = epsilon*epsilon; 
	  }
	else 
	  {
	    epsilonsq = epsilon0sq + (1.- epsilon0sq)*G4UniformRand();
	    epsilon   = sqrt(epsilonsq);
	  }
	
	onecost = (1.- epsilon)/(epsilon*E0_m);
	sint2   = onecost*(2.-onecost);
	greject = 1. - epsilon*sint2/(1.+ epsilonsq);
      } while (greject < G4UniformRand());
      
      
      // ****************************************************
      //		Phi determination
      // ****************************************************
      
      G4double Phi = SetPhi(epsilon,sint2);
      
      //
      // scattered gamma angles. ( Z - axis along the parent gamma)
      //
      
      G4double cosTeta = 1. - onecost, sinTeta = sqrt (sint2);
      G4double dirx = sinTeta*cos(Phi), diry = sinTeta*sin(Phi), dirz = cosTeta ;
      
      //
      // update G4VParticleChange for the scattered gamma 
      //
      
      gammaEnergy1 = epsilon*gammaEnergy0;
      
      // New polarization
      
      G4ThreeVector gammaPolarization1 = SetNewPolarization(epsilon,sint2,Phi,
							    cosTeta);
      
      // Set new direction 
      G4ThreeVector tmpDirection1( dirx,diry,dirz );
      gammaDirection1 = tmpDirection1;
      
      // Change reference frame.
      SystemOfRefChange(gammaDirection0,gammaDirection1,
			gammaPolarization0,gammaPolarization1);   
      
      if (gammaEnergy1 > 0.)
	{
	  aParticleChange.SetEnergyChange( gammaEnergy1 ) ;
	}
      else
	{    
	  aParticleChange.SetEnergyChange(0.) ;
	  aParticleChange.SetStatusChange(fStopAndKill);
	}
      
    }
  else
    {
      // Temporary, same algorithm as G4LowEnergyCompton

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
	
	x = sqrt(onecost/2)/(wlGamma/cm);
	
	const G4FirstLevel* oneAtomSF
	  = (*theScatteringFunctionTable)[ZNumVec->index(elementZ)];
	
	ScatteringFunction = util.DataLogInterpolation(x, (*(*oneAtomSF)[0]), 
						       (*(*oneAtomSF)[1]));
	greject = (1. - epsilon*sint2/(1.+ epsilonsq))*ScatteringFunction;
	
      }  while(greject < G4UniformRand()*elementZ);
      
      G4double cosTeta = 1. - onecost , sinTeta = sqrt (sint2);
      G4double Phi     = twopi * G4UniformRand() ;
      G4double dirx = sinTeta*cos(Phi) , diry = sinTeta*sin(Phi) , dirz = cosTeta ;
      
      //
      // update G4VParticleChange for the scattered gamma 
      //
      
      G4ThreeVector tmpGammaDirection( dirx,diry,dirz );
      gammaDirection1 = tmpGammaDirection;
      gammaDirection1.rotateUz(gammaDirection0);
      aParticleChange.SetMomentumChange( gammaDirection1 ) ;
      gammaEnergy1 = epsilon*gammaEnergy0;
      if (gammaEnergy1 > 0.)
	{
	  aParticleChange.SetEnergyChange( gammaEnergy1 ) ;
	}
      else
	{    
	  aParticleChange.SetEnergyChange(0.) ;
	  aParticleChange.SetStatusChange(fStopAndKill);
	}
      
    }
  
  //
  // kinematic of the scattered electron
  //
  
  G4double ElecKineEnergy = gammaEnergy0 - gammaEnergy1 ;
  
  if((G4EnergyLossTables::GetRange(G4Electron::Electron(),
				   ElecKineEnergy,aMaterial)>aStep.GetPostStepPoint()->GetSafety())
     ||
     (ElecKineEnergy > 
      (G4Electron::Electron()->GetCutsInEnergy())[aMaterial->GetIndex()]))              
    {
      G4double ElecMomentum = sqrt(ElecKineEnergy*(ElecKineEnergy+2.*electron_mass_c2));
      G4ThreeVector ElecDirection (
				   (gammaEnergy0*gammaDirection0 - gammaEnergy1*gammaDirection1)*(1./ElecMomentum) );
      
      // create G4DynamicParticle object for the electron.  
      G4DynamicParticle* aElectron= new G4DynamicParticle (G4Electron::Electron(),
							   ElecDirection, ElecKineEnergy) ;
      
      aParticleChange.SetNumberOfSecondaries(1) ;
      aParticleChange.AddSecondary( aElectron ) ;
      aParticleChange.SetLocalEnergyDeposit (0.) ; 
    }
  else
    {
      aParticleChange.SetNumberOfSecondaries(0) ;
      aParticleChange.SetLocalEnergyDeposit (ElecKineEnergy) ;
    }
  
  
  
  // --- The end ---
  
  //  Reset NbOfInteractionLengthLeft and return aParticleChange
  
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
  
}


G4double G4LowEnergyPolarizedCompton::SetPhi(G4double EnergyRate,
					     G4double sinsqrth)
{
  G4double Rand1;
  G4double Rand2;
  G4double PhiProbability;
  G4double Phi;
  G4double a, b;
  
  do
    {
      Rand1 = G4UniformRand();
      Rand2 = G4UniformRand();
      PhiProbability=0.;
      Phi = twopi*Rand1;
      
      a = 2*sinsqrth;
      b = EnergyRate + 1/EnergyRate;
      
      PhiProbability = 1 - (a/b)*(cos(Phi)*cos(Phi)); 
    }
  while ( Rand2 > PhiProbability );
  
  return Phi;
}


G4ThreeVector G4LowEnergyPolarizedCompton::SetNewPolarization(G4double epsilon, G4double sinsqrth, 
							      G4double phi,G4double costheta) 
{
  G4double Rand1;
  G4double Rand2;
  G4double cosphi = cos(phi), sinphi = sin(phi);
  G4double sintheta = sqrt(sinsqrth);
  G4double cossqrphi = cosphi*cosphi;
  //  G4double cossqrth = 1.-sinsqrth;
  //  G4double sinsqrphi = sinphi*sinphi;
  G4double Normalisation = sqrt(1-cossqrphi*sinsqrth);
  
  // Determination of Theta 
  
  G4double ThetaProbability;
  G4double Theta;
  G4double a, b;

  do
    {
      Rand1 = G4UniformRand();
      Rand2 = G4UniformRand();
      ThetaProbability=0.;
      Theta = twopi*Rand1;
      
      a = 4;
      b = (epsilon + 1/epsilon) - 2;
      
      ThetaProbability = (b + a*cos(Theta)*cos(Theta))/(a+b); 
    }
  while ( Rand2 > ThetaProbability );
  
  G4double cosTheta = cos(Theta);

  G4double cosbeta = cosTheta/Normalisation;
  G4double sinbeta = sqrt(1-cosbeta*cosbeta);
  G4ThreeVector gammaPolarization1;

  G4double Xparallel = Normalisation*cosbeta;
  G4double Yparallel = -(sinsqrth*cosphi*sinphi)*cosbeta/Normalisation;
  G4double Zparallel = -(costheta*sintheta*cosphi)*cosbeta/Normalisation;

  G4double Xperpendicular = 0.;
  G4double Yperpendicular = (costheta)*sinbeta/Normalisation;
  G4double Zperpendicular = -(sintheta*sinphi)*sinbeta/Normalisation;

  G4double Xtotal = (Xparallel + Xperpendicular);
  G4double Ytotal = (Yparallel + Yperpendicular);
  G4double Ztotal = (Zparallel + Zperpendicular);

  gammaPolarization1.setX(Xtotal);
  gammaPolarization1.setY(Ytotal);
  gammaPolarization1.setZ(Ztotal);

  return gammaPolarization1;
}


void G4LowEnergyPolarizedCompton::SystemOfRefChange(G4ThreeVector& direction0,
						    G4ThreeVector& direction1,
						    G4ThreeVector& polarization0,
						    G4ThreeVector& polarization1)

{
  // Angles for go back to the original RS
  G4double cosTeta0 = direction0.cosTheta(), sinTeta0 = sin(direction0.theta());
  G4double cosPhi0  = cos(direction0.phi()), sinPhi0  = sin(direction0.phi());

  G4double cosPsi, sinPsi;

  if (sinTeta0 != 0. ) {
    cosPsi = -polarization0.z()/sinTeta0;
    if (cosPhi0 != 0. ) {
      sinPsi = (polarization0.y() - cosTeta0*sinPhi0*cosPsi)/cosPhi0;
    } else {
      sinPsi = -polarization0.x()/sinPhi0;
    }
  } else {
    cosPsi = polarization0.x()/cosTeta0;
    sinPsi = polarization0.y();
  }
  G4double Psi = atan(sinPsi/cosPsi); 
  

  // Rotation along Z axis
  direction1.rotateZ(Psi);
  // 
  direction1.rotateUz(direction0);
  aParticleChange.SetMomentumChange( direction1 ) ;  

  // 3 Euler angles rotation for scattered photon polarization

  polarization1.rotateZ(Psi);
  polarization1.rotateUz(direction0);
  aParticleChange.SetPolarizationChange( polarization1 );
}


G4bool G4LowEnergyPolarizedCompton::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}


G4double G4LowEnergyPolarizedCompton::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*)
{

  // returns the gamma mean free path in GEANT4 internal units
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  G4double gammaEnergy = aDynamicGamma->GetKineticEnergy();
  G4Material* aMaterial = aTrack.GetMaterial();
  
  //  G4bool isOutRange ;
  
  
  if (gammaEnergy > highestEnergyLimit)
    {
      meanFreePath = DBL_MAX;
    }	
  else if(gammaEnergy < lowestEnergyLimit)
    {
      meanFreePath = DBL_MIN;
    }
  else 
    {
      meanFreePath = util.DataLogInterpolation(gammaEnergy, 
					       aMaterial->GetIndex(), 
					       theMeanFreePathTable);	
    }                                     
  
  return meanFreePath;
}

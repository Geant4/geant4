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
// $Id: G4LowEnergyPolarizedCompton.cc,v 1.6 2001-07-11 10:02:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//

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
    numbBinTable(200),
    meanFreePath(0)
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
  
  for(size_t tableInd = 0; tableInd < ZNumVec->size(); tableInd++){
    
    G4int atomInd = (G4int) (*ZNumVec)[tableInd];
    
    G4FirstLevel* oneAtomCS = util.BuildFirstLevelTables(atomInd, dataNum, "comp/ce-cs-");
    
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
 
  for(size_t tableInd = 0; tableInd < ZNumVec->size(); tableInd++){

    G4int atomInd = (G4int) (*ZNumVec)[tableInd];

    G4FirstLevel* oneAtomSF = util.BuildFirstLevelTables(atomInd, dataNum, "comp/ce-sf-");
     
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
    const G4int numberOfElements = material->GetNumberOfElements() ;

    for (G4int iel=0; iel<numberOfElements; iel++ ){

      G4double Zel = (*theElementVector)(iel)->GetZ();

      if(!(ZNumVec->contains(Zel))){
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
  G4double lowEdgeEnergy, value;
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
      
      lowEdgeEnergy = ptrVector->GetLowEdgeEnergy(i);
      
      const G4double bigPath = DBL_MAX;
      G4double sigma = 0. ;
      for ( size_t k=0 ; k < material->GetNumberOfElements() ; k++ ){ 

	G4int atomIndex = (G4int) (*theElementVector)(k)->GetZ();
	const G4FirstLevel* oneAtomCS
	  = (*theCrossSectionTable)[ZNumVec->index(atomIndex)];
	
	G4double interCrsSec = util.DataLogInterpolation(lowEdgeEnergy, 
							 (*(*oneAtomCS)[0]), 
							 (*(*oneAtomCS)[1]))*barn;
	sigma += theAtomNumDensityVector[k]*interCrsSec;
      }       
      
      value = sigma<=0.0 ? bigPath : 1./sigma ;
      ptrVector->PutValue( i , value ) ;

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
 
  G4double partialSumSigma = 0.;

  G4double rval = 0;
  rval = G4UniformRand()/meanFreePath;

  for ( G4int i=0 ; i < NumberOfElements ; i++ ){ 

    G4double crossSection;
    if (gammaEnergy <  lowestEnergyLimit)
      crossSection = 0. ;
    else {
      if (gammaEnergy > highestEnergyLimit) gammaEnergy = 0.99*highestEnergyLimit ;
      
      G4int atomIndex = (G4int) (*theElementVector)(i)->GetZ();
      const G4FirstLevel* oneAtomCS
	= (*theCrossSectionTable)[ZNumVec->index(atomIndex)];

      crossSection =  util.DataLogInterpolation(gammaEnergy, 
						(*(*oneAtomCS)[0]), 
						(*(*oneAtomCS)[1]))*barn;
    }

    partialSumSigma += theAtomNumDensityVector[i] * crossSection;
    if(rval <= partialSumSigma) return ((*theElementVector)(i));
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

  // Temporary protection: a polarisation parallel to the 
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

  G4ThreeVector gammaDirection0 = aDynamicGamma->GetMomentumDirection();

  G4double epsilon, epsilonSq, onecost, sinThetaSqr, greject ;

  G4double epsilon0 = 1./(1. + 2*E0_m);
  G4double epsilon0Sq = epsilon0*epsilon0;
  G4double alpha1   = - log(epsilon0);
  G4double alpha2 = 0.5*(1.- epsilon0Sq);
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
	    epsilonSq = epsilon*epsilon; 
	  }
	else 
	  {
	    epsilonSq = epsilon0Sq + (1.- epsilon0Sq)*G4UniformRand();
	    epsilon   = sqrt(epsilonSq);
	  }
	
	onecost = (1.- epsilon)/(epsilon*E0_m);
	sinThetaSqr   = onecost*(2.-onecost);

	// Protection
	if (sinThetaSqr > 1.)
	  {
	    if (verboseLevel>0) G4cout 
	      << " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	      << "sin(theta)**2 = " 
	      << sinThetaSqr
	      << "; set to 1" 
	      << G4endl;
	    sinThetaSqr = 1.;
	  }
	if (sinThetaSqr < 0.)
	  {
	    if (verboseLevel>0) G4cout 
	      << " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	      << "sin(theta)**2 = " 
	      << sinThetaSqr
	      << "; set to 0" 
	      << G4endl;
	    sinThetaSqr = 0.;
	  }
        // End protection
	
	greject = 1. - epsilon*sinThetaSqr/(1.+ epsilonSq);
      } while (greject < G4UniformRand());
      
      // ****************************************************
      //		Phi determination
      // ****************************************************
      
      G4double phi = SetPhi(epsilon,sinThetaSqr);

      //
      // scattered gamma angles. ( Z - axis along the parent gamma)
      //
      
      G4double cosTheta = 1. - onecost;

      // Protection
      
      if (cosTheta > 1.)
	{
	  if (verboseLevel>0) G4cout 
	    << " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	    << "cosTheta = " 
	    << cosTheta
	    << "; set to 1" 
	    << G4endl;
	  cosTheta = 1.;
	}
      if (cosTheta < -1.)
	{
	  if (verboseLevel>0) G4cout 
	    << " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	    << "cosTheta = " 
	    << cosTheta
	    << "; set to -1" 
	    << G4endl;
	  cosTheta = -1.;
	}
      // End protection      
      
      
      G4double sinTheta = sqrt (sinThetaSqr);
      
      // Protection
      if (sinTheta > 1.)
	{
	  if (verboseLevel>0) G4cout 
	    << " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	    << "sinTheta = " 
	    << sinTheta
	    << "; set to 1" 
	    << G4endl;
	  sinTheta = 1.;
	}
      if (sinTheta < -1.)
	{
	  if (verboseLevel>0) G4cout 
	    << " -- Warning -- G4LowEnergyPolarizedCompton::PostStepDoIt "
	    << "sinTheta = " 
	    << sinTheta
	    << "; set to -1" 
	    << G4endl;
	  sinTheta = -1.;
	}
      // End protection
      
      
      G4double dirx = sinTheta*cos(phi);
      G4double diry = sinTheta*sin(phi);
      G4double dirz = cosTheta ;
      
      //
      // update G4VParticleChange for the scattered gamma 
      //
      
      gammaEnergy1 = epsilon*gammaEnergy0;
      
      // New polarization


      G4ThreeVector gammaPolarization1 = SetNewPolarization(epsilon,
							    sinThetaSqr,
							    phi,
							    cosTheta);
      
      // Set new direction 
      //G4ThreeVector tmpDirection1( dirx,diry,dirz );
      G4ParticleMomentum tmpDirection1( dirx,diry,dirz );
      
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
	  epsilonSq = epsilon*epsilon; 
	}
	else{
	  
	  epsilonSq = epsilon0Sq + (1.- epsilon0Sq)*G4UniformRand();
	  epsilon   = sqrt(epsilonSq);
	}
	
	onecost = (1.- epsilon)/(epsilon*E0_m);
	sinThetaSqr   = onecost*(2.-onecost);
	
	x = sqrt(onecost/2)/(wlGamma/cm);
	
	const G4FirstLevel* oneAtomSF
	  = (*theScatteringFunctionTable)[ZNumVec->index(elementZ)];
	
	ScatteringFunction = util.DataLogInterpolation(x, (*(*oneAtomSF)[0]), 
						       (*(*oneAtomSF)[1]));
	greject = (1. - epsilon*sinThetaSqr/(1.+ epsilonSq))*ScatteringFunction;
	
      }  while(greject < G4UniformRand()*elementZ);
      
      G4double cosTheta = 1. - onecost ;
      G4double sinTheta = sqrt (sinThetaSqr);
      G4double phi     = twopi * G4UniformRand() ;
      G4double dirx = sinTheta*cos(phi) , diry = sinTheta*sin(phi) , dirz = cosTheta ;
      
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
							   ElecDirection, 
							   ElecKineEnergy) ;
      
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


G4double G4LowEnergyPolarizedCompton::SetPhi(G4double energyRate,
					     G4double sinSqrTh)
{
  G4double rand1;
  G4double rand2;
  G4double phiProbability;
  G4double phi;
  G4double a, b;
  
  do
    {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      phiProbability=0.;
      phi = twopi*rand1;
      
      a = 2*sinSqrTh;
      b = energyRate + 1/energyRate;
      
      phiProbability = 1 - (a/b)*(cos(phi)*cos(phi));

      
 
    }
  while ( rand2 > phiProbability );
  return phi;
}


G4ThreeVector G4LowEnergyPolarizedCompton::SetNewPolarization(G4double epsilon, 
							      G4double sinSqrTh, 
							      G4double phi,
							      G4double costheta) 
{
  G4double rand1;
  G4double rand2;
  G4double cosPhi = cos(phi);
  G4double sinPhi = sin(phi);
  G4double sinTheta = sqrt(sinSqrTh);
  G4double cosSqrPhi = cosPhi*cosPhi;
  //  G4double cossqrth = 1.-sinSqrTh;
  //  G4double sinsqrphi = sinPhi*sinPhi;
  G4double normalisation = sqrt(1. - cosSqrPhi*sinSqrTh);
  

  // Determination of Theta 
  
  G4double thetaProbability;
  G4double theta;
  G4double a, b;
  G4double cosTheta;

  do
    {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      thetaProbability=0.;
      theta = twopi*rand1;
      a = 4;
      b = (epsilon + 1/epsilon) - 2;
      thetaProbability = (b + a*cos(theta)*cos(theta))/(a+b);
      cosTheta = cos(theta);       
    }
  while ( rand2 > thetaProbability || abs(cosTheta) > abs(normalisation) );
  
  G4double cosBeta = cosTheta/normalisation;
  G4double sinBeta = sqrt(1-cosBeta*cosBeta);

  G4ThreeVector gammaPolarization1;
  
  G4double xParallel = normalisation*cosBeta;
  G4double yParallel = -(sinSqrTh*cosPhi*sinPhi)*cosBeta/normalisation;
  G4double zParallel = -(cosTheta*sinTheta*cosPhi)*cosBeta/normalisation;
  G4double xPerpendicular = 0.;
  G4double yPerpendicular = (cosTheta)*sinBeta/normalisation;
  G4double zPerpendicular = -(sinTheta*sinPhi)*sinBeta/normalisation;
  
  G4double xTotal = (xParallel + xPerpendicular);
  G4double yTotal = (yParallel + yPerpendicular);
  G4double zTotal = (zParallel + zPerpendicular);
  
  gammaPolarization1.setX(xTotal);
  gammaPolarization1.setY(yTotal);
  gammaPolarization1.setZ(zTotal);
  
  return gammaPolarization1;

}


void G4LowEnergyPolarizedCompton::SystemOfRefChange
(G4ThreeVector& direction0,G4ThreeVector& direction1,
 G4ThreeVector& polarization0,G4ThreeVector& polarization1)
  
{
  // Angles for go back to the original RS


  
  G4double cosTheta0 = direction0.cosTheta();
  G4double sinTheta0 = sin(direction0.theta());
  G4double cosPhi0 = cos(direction0.phi());
  G4double sinPhi0 = sin(direction0.phi());
  G4double cosPsi, sinPsi;
  
  //  if (sinTheta0 >= 1.e-14 ) {
  if (sinTheta0 != 0. ) 
    {
      cosPsi = -polarization0.z()/sinTheta0;
      //  if (cosPhi0 >=1.e-14 ) {
      if (cosPhi0 != 0. ) 
	{
	  sinPsi = (polarization0.y() - cosTheta0*sinPhi0*cosPsi)/cosPhi0;
	} 
      else 
	{
	  sinPsi = -polarization0.x()/sinPhi0;
	}
    } 
  else 
    {
      cosPsi = polarization0.x()/cosTheta0;
      sinPsi = polarization0.y();
    }
  
  // Added protection
  G4double psi = 0;
  if (sinPsi < 0.) psi = -pi/2.;
  if (sinPsi > 0.) psi = pi/2.;

  if (cosPsi != 0.)
    {
      psi = atan(sinPsi/cosPsi); 
    }
  
  // Rotation along Z axis
  direction1.rotateZ(psi);
  // 
  direction1.rotateUz(direction0);

  aParticleChange.SetMomentumChange( direction1 ) ;  

  // 3 Euler angles rotation for scattered photon polarization
  
  polarization1.rotateZ(psi);
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
















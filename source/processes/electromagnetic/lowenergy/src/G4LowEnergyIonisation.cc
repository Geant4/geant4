// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyIonisation.cc,v 1.19 1999-07-05 14:27:34 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4LowEnergyIonisation physics process -----------
//                by Laszlo Urban, 20 March 1997 
// **************************************************************
// It is the first implementation of the NEW IONISATION PROCESS.
// It calculates the ionisation of e+/e-.
// **************************************************************
//
// 07-04-98: remove 'tracking cut' of the ionizing particle, MMa 
// 04-09-98: new methods SetBining() PrintInfo()
// 07-09-98: Cleanup
// --------------------------------------------------------------
 
// This Class Header
#include "G4LowEnergyIonisation.hh"

// Collaborating Class Headers
#include "G4EnergyLossTables.hh"
#include "G4Gamma.hh"
#include "G4UnitsTable.hh"

typedef RWTPtrOrderedVector<G4DynamicParticle> G4ParticleVector;

// constructor and destructor
G4LowEnergyIonisation::G4LowEnergyIonisation(const G4String& processName)
   : G4eEnergyLoss(processName),
     allAtomShellCrossSec(0),
     theBindingEnergyTable(0),
     theFluorTransitionTable(0),
     theSamplingCoeffTable(0),
     LowestKineticEnergy(250.*eV),
     HighestKineticEnergy(100.*GeV),
     CutForLowEnergySecondaryPhotons(0.),
     CutForLowEnergySecondaryElectrons(0.),
     ZNumVec(0),
     ZNumVecFluor(0),
     TotBin(200)
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LowEnergyIonisation::~G4LowEnergyIonisation() 
{

  if (allAtomShellCrossSec) {

    delete allAtomShellCrossSec;
  }

  if (theBindingEnergyTable) {

    delete theBindingEnergyTable;
  }

  if (theFluorTransitionTable) {

    delete theFluorTransitionTable;
  }

  if(theSamplingCoeffTable){

    delete theSamplingCoeffTable;
  }

  if(ZNumVec){

    ZNumVec->clear();
    delete ZNumVec;
  }

  if(ZNumVecFluor){

    ZNumVecFluor->clear();
    delete ZNumVecFluor;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4LowEnergyIonisation::SetCutForLowEnSecPhotons(G4double cut){

  CutForLowEnergySecondaryPhotons = cut;
}

void G4LowEnergyIonisation::SetCutForLowEnSecElectrons(G4double cut){

  CutForLowEnergySecondaryElectrons = cut;
  //  LowestKineticEnergy = 2*cut;
}

void G4LowEnergyIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
//  just call BuildLossTable+BuildLambdaTable
{
 
    BuildLossTable(aParticleType) ;

  if(&aParticleType==G4Electron::Electron())
  {
    RecorderOfElectronProcess[CounterOfElectronProcess] = (*this).theLossTable ;
    CounterOfElectronProcess++;
  }
  else
  {
    RecorderOfPositronProcess[CounterOfPositronProcess] = (*this).theLossTable ;
    CounterOfPositronProcess++;
  }
 
  BuildDEDXTable(aParticleType);
  
  BuildZVec();

  BuildShellCrossSectionTable();
  
  BuildFluorTransitionTable();
  
  BuildBindingEnergyTable();
  
  BuildSamplingCoeffTable();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
// Build tables for the ionization energy loss
//  the tables are built for *MATERIALS*

    const G4double twoln10 = 2.*log(10.);
    const G4double Factor = twopi_mc2_rcl2;

    G4double LowEdgeEnergy, ionloss;
    
    // material properties
    G4double ElectronDensity,Eexc,Eexcm2,Cden,Mden,Aden,X0den,X1den ;
    // some local variables
    G4double tau,Tmax,gamma,gamma2,bg2,beta2,d,d2,d3,d4,delta,x,y ;

    ParticleMass = aParticleType.GetPDGMass();
    G4double* ParticleCutInKineticEnergy = aParticleType.GetEnergyCuts() ;

    //  create table
    
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    G4int numOfMaterials = theMaterialTable->length();

     if (theLossTable) { theLossTable->clearAndDestroy();
                         delete theLossTable;
                       }
                       
    theLossTable = new G4PhysicsTable(numOfMaterials);

//  loop for materials

    for (G4int J=0; J<numOfMaterials; J++)
    {
      // create physics vector and fill it

      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);
  
      // get material parameters needed for the energy loss calculation   
      const G4Material* material= (*theMaterialTable)[J];

      ElectronDensity = material->GetElectronDensity();
      Eexc   = material->GetIonisation()->GetMeanExcitationEnergy();
      Eexc  /= ParticleMass; Eexcm2 = Eexc*Eexc;
      Cden   = material->GetIonisation()->GetCdensity();
      Mden   = material->GetIonisation()->GetMdensity();
      Aden   = material->GetIonisation()->GetAdensity();
      X0den  = material->GetIonisation()->GetX0density();
      X1den  = material->GetIonisation()->GetX1density();

      // now comes the loop for the kinetic energy values

      for (G4int i = 0 ; i < TotBin ; i++)
         {
          LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;      
          tau = LowEdgeEnergy/ParticleMass ;

          // Seltzer-Berger formula 
          gamma = tau + 1.; gamma2 = gamma*gamma; 
          bg2 = tau*(tau+2.);
          beta2 = bg2/gamma2;

          // electron
          if (&aParticleType==G4Electron::Electron())
            {
              Tmax = LowEdgeEnergy/2.;  
              d = min(ParticleCutInKineticEnergy[J], Tmax)/ParticleMass;
              ionloss = log(2.*(tau+2.)/Eexcm2)-1.-beta2
                       + log((tau-d)*d)+tau/(tau-d)
                       + (0.5*d*d+(2.*tau+1.)*log(1.-d/tau))/gamma2;
            }
          else        //positron
            {
              Tmax = LowEdgeEnergy ;  
              d = min(ParticleCutInKineticEnergy[J], Tmax)/ParticleMass;
              d2=d*d/2.; d3=d*d*d/3.; d4=d*d*d*d/4.;
              y=1./(1.+gamma);
              ionloss = log(2.*(tau+2.)/Eexcm2)+log(tau*d)
                       - beta2*(tau+2.*d-y*(3.*d2+y*(d-d3+y*(d2-tau*d3+d4))))/tau;
            } 

          //density correction
          x = log(bg2)/twoln10;
          if (x < X0den) delta = 0.;
          else { delta = twoln10*x - Cden;
                 if (x < X1den) delta += Aden*pow((X1den-x),Mden);
               } 

          //now you can compute the total ionization loss
          ionloss -= delta ;
          ionloss *= Factor*ElectronDensity/beta2 ;
          if (ionloss <= 0.) ionloss = 0.;
   
          aVector->PutValue(i,ionloss) ;
         }          

      theLossTable->insert(aVector);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::BuildShellCrossSectionTable(){

   if (allAtomShellCrossSec) {

    delete allAtomShellCrossSec;
   }

   allAtomShellCrossSec = new allAtomTable();
   G4int dataNum = 2;
 
   for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){

     G4int AtomInd = (G4int) (*ZNumVec)[TableInd];

     oneAtomTable* oneAtomShellCS = util.BuildSecondLevelTables(AtomInd, dataNum, "ioni/ion-ss-cs-");
     
     allAtomShellCrossSec->insert(oneAtomShellCS);
   
   }//end for on atoms
}
void G4LowEnergyIonisation::BuildBindingEnergyTable(){

  if (theBindingEnergyTable) {

    delete theBindingEnergyTable;
  }

  G4int dataNum = 2;
  theBindingEnergyTable = util.BuildSecondLevelTables(0,dataNum,"fluor/binding");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::BuildFluorTransitionTable(){

  if (theFluorTransitionTable) {
    
    delete theFluorTransitionTable;
  }
  
  theFluorTransitionTable = new allAtomTable();
  ZNumVecFluor = new G4Data(*ZNumVec);
  G4int dataNum = 3;
  
  for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    if(AtomInd > 5){
      
      oneAtomTable* oneAtomShellFL = util.BuildSecondLevelTables(AtomInd, dataNum, "fluor/fl-tr-pr-");
      theFluorTransitionTable->insert(oneAtomShellFL);
    }
    else{
      ZNumVecFluor->remove(AtomInd);
    }
  }//end for on atoms
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::BuildSamplingCoeffTable(){

  if (theSamplingCoeffTable) {
    
    delete theSamplingCoeffTable;
  }
  
  theSamplingCoeffTable = new allAtomTable();
  
  G4int dataNum = 12;

  for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    oneAtomTable* oneAtomShellSc = util.BuildSecondLevelTables(AtomInd, dataNum, "ioni/ion-co-");
    
    theSamplingCoeffTable->insert(oneAtomShellSc);

  }//end for on atoms
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::BuildZVec(){

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

G4double G4LowEnergyIonisation::ComputeCrossSection(const G4double AtomIndex,
						    const G4double IncEnergy){
  // calculates the microscopic cross section from subshell cross sections
  //(it is called for elements , AtomicNumber = Z )
 
  G4double TotalCrossSection(0.);

  const oneAtomTable* oneAtomCS
    = (*allAtomShellCrossSec)[ZNumVec->index(AtomIndex)];

  for(G4int ind = 0; ind < oneAtomCS->entries(); ind++){

    G4double crossSec = 0;
    G4Data* EnergyVector = (*(*oneAtomCS)[ind])[0];
    G4Data* CrossSecVector = (*(*oneAtomCS)[ind])[1];

    if(IncEnergy < (*EnergyVector)[1]){ // First element is the shell number

      crossSec = 0;
    }

    else{

      crossSec = util.DataLogInterpolation(IncEnergy, (*EnergyVector), (*CrossSecVector));

    }

    TotalCrossSection += crossSec;
  }

  return TotalCrossSection ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
 
G4VParticleChange* G4LowEnergyIonisation::AlongStepDoIt( const G4Track& trackData,
                                                 const G4Step&  stepData)
{

  static const G4double AbsCut = LowestKineticEnergy;

  // get particle energy
  G4double E = trackData.GetDynamicParticle()->GetKineticEnergy();

  aParticleChange.Initialize(trackData);

  if(E <= AbsCut)
  {
    aParticleChange.SetStatusChange(fStopAndKill);
    aParticleChange.SetEnergyChange(0.);
    aParticleChange.SetLocalEnergyDeposit(E);
  }
  return &aParticleChange;
}


G4VParticleChange* G4LowEnergyIonisation::PostStepDoIt( const G4Track& trackData,   
							const G4Step&  stepData){

  aParticleChange.Initialize(trackData);
  
  G4Material*               aMaterial = trackData.GetMaterial() ;
  const G4DynamicParticle*  aParticle = trackData.GetDynamicParticle() ;

  // select randomly one element constituing the material.
  G4Element* anElement = SelectRandomAtom(aParticle, aMaterial);
  G4int AtomIndex = (G4int) anElement->GetZ();
  G4double KineticEnergy = aParticle->GetKineticEnergy();

  // Select the subshell WARNING!!!!
  G4int subShellIndex = SelectRandomShell(AtomIndex, KineticEnergy);

  G4FirstLevel* theBindEnVec = (*theBindingEnergyTable)[AtomIndex-1];
  G4int thePrimaryShell = (G4int) (*(*theBindEnVec)[0])[subShellIndex];
  G4double BindingEn = (*(*theBindEnVec)[1])[subShellIndex];
  G4double theEnergyDeposit = BindingEn; 

  G4double charge = aParticle->GetDefinition()->GetPDGCharge();
  ParticleMass = aParticle->GetDefinition()->GetPDGMass();

  G4double Psquare = KineticEnergy*(KineticEnergy+2*ParticleMass);
  G4double TotalMomentum = sqrt(Psquare);
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection();

  //  get kinetic energy cut for the electron
  G4double* DeltaCutInKineticEnergy = G4Electron::Electron()->GetCutsInEnergy() ;
  G4double  DeltaThreshold = DeltaCutInKineticEnergy[aMaterial->GetIndex()];
 
  // some kinematics
  G4double MaxKineticEnergyTransfer;

  if (charge < 0.) MaxKineticEnergyTransfer = 0.5*(KineticEnergy);
  else             MaxKineticEnergyTransfer =     KineticEnergy;

  // sampling kinetic energy of the delta ray 
  if (MaxKineticEnergyTransfer <= 0 || MaxKineticEnergyTransfer <= LowestKineticEnergy/2){

    // pathological case (should not happen, there is no change at all)
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }

  // **** normal case ****

  //Energy Sampling
  G4double DeltaKineticEnergy = EnergySampling(AtomIndex, subShellIndex, KineticEnergy);

  // protection :do not produce a secondary with 0. kinetic energy !
  if (DeltaKineticEnergy <= 0.){

    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  } 

  if(DeltaKineticEnergy <= DeltaThreshold){

    aParticleChange.SetNumberOfSecondaries(0);
    aParticleChange.SetEnergyChange(KineticEnergy - DeltaKineticEnergy - BindingEn);
    aParticleChange.SetLocalEnergyDeposit(DeltaKineticEnergy+BindingEn);
    
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  } 

  //  G4double finalKineticEnergy = KineticEnergy - DeltaKineticEnergy;
  G4double finalKineticEnergy = KineticEnergy - DeltaKineticEnergy - BindingEn;
    
  if(thePrimShVec.length() != 0){
    
    thePrimShVec.clear();
  }
  
  thePrimShVec.insert(thePrimaryShell);
  
  // delta ray kinematics
  G4double DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
							   2. * electron_mass_c2 ));
  
  if(finalKineticEnergy > 0.){
    
    G4double finalMomentum=sqrt(finalKineticEnergy*
				(finalKineticEnergy+2.*ParticleMass));
    
    G4double costheta = (Psquare-(finalMomentum*finalMomentum)+
			 (DeltaTotalMomentum*DeltaTotalMomentum))/(2*DeltaTotalMomentum*TotalMomentum); 
 
    G4double costhetasc = (Psquare+(finalMomentum*finalMomentum)-
			   (DeltaTotalMomentum*DeltaTotalMomentum))/(2*finalMomentum*TotalMomentum); 

    if (costheta < -1.) costheta = -1.;
    if (costheta > +1.) costheta = +1.;
    
    //  direction of the delta electron
    G4double phi = twopi * G4UniformRand(); 
    G4double sintheta = sqrt((1.+costheta)*(1.-costheta));
    
    G4double dirx = sintheta * cos(phi), diry = sintheta * sin(phi), dirz = costheta;
    
    G4ThreeVector DeltaDirection(dirx,diry,dirz);
    DeltaDirection.rotateUz(ParticleDirection);
    
    // finalKineticEnergy and finalMomentum defined above 
    // because needed for costheta computation  
    
    G4double finalPx = (TotalMomentum*ParticleDirection.x()
			- DeltaTotalMomentum*DeltaDirection.x())/finalMomentum; 

    G4double finalPy = (TotalMomentum*ParticleDirection.y()
			- DeltaTotalMomentum*DeltaDirection.y())/finalMomentum; 

    G4double finalPz = (TotalMomentum*ParticleDirection.z()
			- DeltaTotalMomentum*DeltaDirection.z())/finalMomentum; 
    
    G4double momtot = sqrt(finalPx*finalPx + finalPy*finalPy + finalPz*finalPz);
    


    if(momtot-1. > 1e-6){
      
      finalPx /= momtot; finalPy /= momtot; finalPz /= momtot;

    }
    
    // Create lists of pointers to DynamicParticles (photons and electrons) 
    G4ParticleVector photvec;
    G4int photInd = 0; 
    G4ParticleVector elecvec;
    G4int elecInd = 0; 
    
    // create G4DynamicParticle object for delta ray
    G4DynamicParticle* theDeltaRay = new G4DynamicParticle;
    theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
    theDeltaRay->SetMomentumDirection(DeltaDirection.x(),
				      DeltaDirection.y(),
				      DeltaDirection.z()); 
    
    theDeltaRay->SetDefinition(G4Electron::Electron());
    elecvec.insert(theDeltaRay);
    
    // FLUORESCENCE
    // load the transition probability table for the element
    // theTable[i][j][k] 
    // i = subshell, j = type of information (second shell, transition energy , 
    // transition probability), k = previous vectors.
    
    // Fluorescence data start from element 6
    
    if(AtomIndex > 5){
      
      G4bool ThereAreShells = TRUE;
      G4int AtomInd = ZNumVecFluor->index(AtomIndex);
      oneAtomTable* oneAtomFluorTrans = (*theFluorTransitionTable)[AtomInd];
      
      while(ThereAreShells == TRUE){
	
	// Select the second transition from another subshell
	// fluorPar[0] = SubShell 
	// fluorPar[1] = Sec SubShell (if there is), 
	// fluorPar[2] = Transition Probability
	// fluorPar[3] = Transition Energy
	// the same for augerPar
	
	G4double fluorPar[3] = {0};
	
	// SelectRandomTransition argument is oneAtomTable loop on shells is inside
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
	newPartDirection.rotateUz(ParticleDirection);
	
	if(ThereAreShells != FALSE){
	  
	  thePrimaryShell = (G4int) fluorPar[0];
	  theEnergyDeposit -= fluorPar[2]*MeV;
	  
	  if(fluorPar[2] >= CutForLowEnergySecondaryPhotons){
	    
	    newPart = new G4DynamicParticle (G4Gamma::Gamma(), 
					     newPartDirection, 
					     fluorPar[2]);
	    
	    photvec.insert(newPart);
	  }
	}
	else{
	  
	  // last shell transition from continuum
	  G4int k = 0;
	  while(thePrimaryShell != (*(*theBindEnVec)[0])[k]){
	    k++;
	  }
	  
	  G4double lastTransEnergy = (*(*theBindEnVec)[1])[k];
	  thePrimaryShell = (G4int) fluorPar[0];
	  
	  if(fluorPar[2] >= CutForLowEnergySecondaryPhotons){
	    
	    theEnergyDeposit -= lastTransEnergy*MeV;
	    
	    newPart = new G4DynamicParticle(G4Gamma::Gamma(), 
					    newPartDirection, 
					    lastTransEnergy);
	    
	    photvec.insert(newPart);
	  }
	    
	  thePrimShVec.insert(thePrimaryShell);
	}
      }
    } //END OF THE CHECK ON ATOMIC NUMBER

    G4int numOfElec = elecvec.entries(), numOfPhot = photvec.entries();
    G4int numOfDau = numOfElec + numOfPhot;
    aParticleChange.SetNumberOfSecondaries(numOfDau);
    G4int l = 0;
    for(l = 0; l<numOfElec; l++ ){

      aParticleChange.AddSecondary(elecvec[l]);
   }
    
    for(l = 0; l < numOfPhot; l++) {
      
      aParticleChange.AddSecondary(photvec[l]); 
    }

    photvec.clear();
    elecvec.clear();
    
    // fill aParticleChange 
    // changed energy and momentum of the actual particle
    if(theEnergyDeposit < 0){
      
      theEnergyDeposit = 0;
    }
    
    aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
    aParticleChange.SetEnergyChange( finalKineticEnergy );
    aParticleChange.SetLocalEnergyDeposit (theEnergyDeposit);
  }
    
  else{
    
    finalKineticEnergy = 0.;
    
    if (charge < 0.) aParticleChange.SetStatusChange(fStopAndKill);
    else             aParticleChange.SetStatusChange(fStopButAlive);
  }
    
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::Print()
{
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4LowEnergyIonisation::SelectRandomShell(const G4int AtomIndex, const G4double IncEnergy){
  
  G4double partialSum = 0;
  G4double totalSum = ComputeCrossSection(AtomIndex,IncEnergy);

  G4double rval = totalSum*G4UniformRand();
  const oneAtomTable* oneAtomCS 
    = (*allAtomShellCrossSec)[ZNumVec->index(AtomIndex)];

  for(G4int ind = 0; ind < oneAtomCS->entries(); ind++){

    G4double crossSec;
    G4Data* EnergyVector = (*(*oneAtomCS)[ind])[0];
    G4Data* CrossSecVector = (*(*oneAtomCS)[ind])[1];
    if(IncEnergy < (*EnergyVector)[0]){ //First element is the shell number

      crossSec = 0;
    }

    else{

      crossSec = util.DataLogInterpolation(IncEnergy, (*EnergyVector), (*CrossSecVector));

    }
    
    partialSum += crossSec;

    if(rval <= partialSum) return ind;
  }
  
  G4Exception("LEIonisation: Cannot select a shell");
  return 0;
}


G4Element*
G4LowEnergyIonisation::SelectRandomAtom(const G4DynamicParticle* aDynamicParticle, 
					G4Material* aMaterial){

  // select randomly 1 element within the material
  G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  if (NumberOfElements == 1) return (*theElementVector)(0);

  const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();

  G4double PartialSumSigma = 0. ;

  // G4int materialIndex = aMaterial->GetIndex();

  G4double rval = G4UniformRand()/MeanFreePath;
 
  for ( G4int i=0 ; i < NumberOfElements ; i++ ){ 

    G4double crossSection;
    if (KineticEnergy <  LowestKineticEnergy)

      crossSection = 0. ;

    else {

      if (KineticEnergy > HighestKineticEnergy) KineticEnergy = 0.99*HighestKineticEnergy;

      G4int AtomIndex = (G4int) (*theElementVector)(i)->GetZ();
      crossSection = ComputeCrossSection(AtomIndex, KineticEnergy);
    }

    PartialSumSigma += theAtomNumDensityVector[i] * crossSection;

    if (rval <= PartialSumSigma) return ((*theElementVector)(i));

  }
  return (*theElementVector)(0);
}

G4bool G4LowEnergyIonisation::SelectRandomTransition(G4int thePrimShell, 
						     G4double* TransParam,
						     const oneAtomTable* TransitionTable){
  
  G4int SubShellCol = 0, ProbCol = 1, EnergyCol = 2;

  //transitionTable means  for one atom not for one shell

  // too check when the subshell are finished
  G4bool ColIsFull = TRUE;
  G4int ShellNum = 0;
  G4double TotalSum = 0; 
  G4int maxNumOfShells = TransitionTable->entries()-1;
  if(thePrimShell <= (*(*(*TransitionTable)[maxNumOfShells])[0])[0]){

    while(thePrimShell != (*(*(*TransitionTable)[ShellNum])[0])[0]){
  
      if(ShellNum == maxNumOfShells){
	break;
      }
      
      ShellNum++;
    }
    
    //TransProb start from 1 because the first element of the list is the primary shall id number
    G4int TransProb = 1;
    for(TransProb = 1; TransProb < (*(*TransitionTable)[ShellNum])[ProbCol]->length(); TransProb++){ 
      
      TotalSum += (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
    }
    
    G4double PartialProb = G4UniformRand()*TotalSum;
    G4double PartSum = 0;
    
    TransProb = 1; 
    while(TransProb < (*(*TransitionTable)[ShellNum])[ProbCol]->length()){
      
      PartSum += (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
      
      if(PartialProb <= PartSum){
	
	TransParam[0] = (*(*(*TransitionTable)[ShellNum])[SubShellCol])[TransProb];
	TransParam[1] = (*(*(*TransitionTable)[ShellNum])[ProbCol])[TransProb];
	TransParam[2] = (*(*(*TransitionTable)[ShellNum])[EnergyCol])[TransProb];
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

G4double G4LowEnergyIonisation::EnergySampling(const G4int AtomicNumber, 
					       const G4int ShellIndex, 
					       const G4double KinEn){

  //  1) Load Coefficients (I need Z number and the index of the shell)
  oneAtomTable* oneAtomCoeffTable = (*theSamplingCoeffTable)[ZNumVec->index(AtomicNumber)];

  oneShellTable* oneShellCoeffTable = (*oneAtomCoeffTable)[ShellIndex];
  G4double BindingEn = (*(*(*theBindingEnergyTable)[AtomicNumber-1])[1])[ShellIndex];

  // 2) Interpolate coefficients (I need the incoming electron kinetic energy)
  const G4int CoeffNumber = oneShellCoeffTable->entries();

  const G4Data* energyVec = (*oneShellCoeffTable)[0];
  const G4int LastPar = energyVec->length()-1;
  G4Data Parms;

  for(G4int ind = 1; ind < CoeffNumber-1; ind++){

    const G4Data* oneCoeffVec = (*oneShellCoeffTable)[ind]; 

    if(KinEn < (*energyVec)[0]){
      Parms.insert((*oneCoeffVec)[0]);
    }

    else if(KinEn > (*energyVec)[LastPar]){

      Parms.insert((*oneCoeffVec)[LastPar]);
    }

    else{

      G4double par = util.DataSemiLogInterpolation(KinEn,(*energyVec),(*oneCoeffVec));
      Parms.insert(par);
    }
  }

  // cut in energy is always the same
  Parms.insert((*(*oneShellCoeffTable)[CoeffNumber-1])[0]);

  // 2') order of parameters:
  //     * Parms[0] = par1 LET
  //     * Parms[1] = par2 LET
  //     * Parms[2] = par3 LET
  //     * Parms[3] = par4 LET
  //     * Parms[4] = par5 LET
  //     * Parms[5] = par6 LET
  //     * Parms[6] = par1 HET
  //     * Parms[7] = max rejection function: g(x)
  //     * Parms[8] = area1
  //     * Parms[9] = area2
  //     * Parms[10] = cut in energy

  // 3) Compute partial areas (with functions here the cut is used)

  // minimum energy that can take an ejected electron
  const G4double minEn = 0.1*eV;

  const G4double argmax = 1/(BindingEn+Parms[10]);
  const G4double argmin = 1/(minEn+BindingEn);
  const G4double area1 = Parms[8];
  //Parms[0]*log(argmin/argmax)+Parms[1]*(argmin-argmax)+
  //2*Parms[2]*(pow(argmin,2)-pow(argmax,2))+3*Parms[3]*(pow(argmin,3)-pow(argmax,3))+
  //4*Parms[4]*(pow(argmin,4)-pow(argmax,4))+5*Parms[5]*(pow(argmin,5)-pow(argmax,5));

  const G4double maxEn = (KinEn-BindingEn)/2;

  G4double area2;
  if(maxEn >= Parms[10]){

    area2 = Parms[9];
  }
  else{

    area2 = 0;
  }

  G4double areaTot = area1+area2;
  G4int which;
  // 4) Generate a random number .to select the region of work
  G4double rand1 = areaTot*G4UniformRand();

  // 5) Sampling
  G4double sample = 0;

  if(rand1 < area1){
    // Low energy transfer 
    G4double rejection = 0;
    which =1;
    do{
      
      G4double rand2 = G4UniformRand();
      G4double Ka = 0;

      if(Parms[10] < maxEn){

	Ka = (BindingEn + Parms[10])/(minEn+BindingEn);
      }

      else{

	Ka = (BindingEn + maxEn)/(minEn+BindingEn);
      }

      sample = (minEn + BindingEn)*pow(Ka,rand2)-BindingEn;

      G4double arg = sample + BindingEn;

      rejection = Parms[0]/arg+Parms[1]/pow(arg,2)+Parms[2]/pow(arg,3)+
	Parms[3]/pow(arg,4)+Parms[4]/pow(arg,5)+Parms[5]/pow(arg,6);

      rejection /= Parms[7];
	
    }while(rejection < G4UniformRand());
  }

  else if(area1 < rand1 && rand1 < areaTot){
    which = 2;
     // High energy transfer
      G4double Norm = (1/Parms[10])-(1/maxEn);
      G4double rand2 = Norm*G4UniformRand();
      sample = 1/((1/Parms[10])-rand2);
     
  }
  //cout<<"ShellIndex: "<<ShellIndex<<" ShellId: "<<(*(*(*theBindingEnergyTable)[AtomicNumber-1])[0])[ShellIndex]<<endl;

  return sample;
}







































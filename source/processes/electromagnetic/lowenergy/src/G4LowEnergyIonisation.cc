// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyIonisation.cc,v 1.25 1999-11-11 15:36:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4LowEnergyIonisation physics process --------
//                   by Michel Maire, April 1996
//      ---------- G4LowEnergyIonisation low energy modifications -----------
//                by Alessandra Forti May 1999  
// **************************************************************
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Added EnergySampling method A. Forti
// Modified PostStepDoIt to insert sampling with EEDL data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
//                                                                 
// --------------------------------------------------------------
 
// This Class Header
#include "G4LowEnergyIonisation.hh"

// Collaborating Class Headers
#include "G4EnergyLossTables.hh"
#include "G4Gamma.hh"
#include "G4UnitsTable.hh"
#include <fstream.h>

typedef G4RWTPtrOrderedVector<G4DynamicParticle> G4ParticleVector;

// constructor and destructor
G4LowEnergyIonisation::G4LowEnergyIonisation(const G4String& processName)
   : G4eEnergyLoss(processName),
     allAtomShellCrossSec(0),
     theBindingEnergyTable(0),
     theFluorTransitionTable(0),
     theSamplingCoeffTable(0),
     LowestKineticEnergy(10.*eV),
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

  PrintInfoDefinition();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// METHOD BELOW  FROM STANDARD E_M PROCESSES CODELEFT BUT AT THE MOMENT NOT USED
void G4LowEnergyIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType) {

    
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
  
  for (G4int J=0; J<numOfMaterials; J++){
    
      // create physics vector and fill it
      
      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy,
							   HighestKineticEnergy,
							   TotBin);
      
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

// CONSTRUCT THE CS TABLE FOR THE ELEMENTS MAPPED IN ZNUMVEC USING EEDL DATA
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
// CONSTRUCT THE CROSS SECTION TABLE FOR ALL THE ELEMENTS
void G4LowEnergyIonisation::BuildBindingEnergyTable(){

  if (theBindingEnergyTable) {

    delete theBindingEnergyTable;
  }

  G4int dataNum = 2;
  theBindingEnergyTable = util.BuildSecondLevelTables(0,dataNum,"fluor/binding");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// CONSTRUCT THE TABLE FOR THE ELEMENTS MAPPED IN ZNUMVECFLUOR DIFFERENT 
// FROM ZNUMVEC BECAUSE THERE IS NO FLUORESCENCE FOR THE FIRST 5 ELEMENTS. 
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
// BUILD THE TABLE OF COEFFICIENT FOR THE ELEMENTS MAPPED IN ZNUMVEC
void G4LowEnergyIonisation::BuildSamplingCoeffTable(){

  if (theSamplingCoeffTable) {
    
    delete theSamplingCoeffTable;
  }
  
  theSamplingCoeffTable = new allAtomTable();
  
  G4int dataNum = 17;

  for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    oneAtomTable* oneAtomShellSc = util.BuildSecondLevelTables(AtomInd, dataNum, "ioni/io-co-");
    
    theSamplingCoeffTable->insert(oneAtomShellSc);

  }//end for on atoms
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//
// vector mapping the elements of the material table 
// 
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
  //(it is called for elements , AtomicNumber = Z ) all the subshell
  // cross section are computed at the energy IncEnergy with a log-log 
  // interpolation as reported by stepanek and then summed to get the 
  // total cross section 
 
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
 
G4VParticleChange* G4LowEnergyIonisation::PostStepDoIt( const G4Track& trackData,   
							const G4Step&  stepData){

  // Fluorescence (as reported by stepanek):
  // J. Stepanek " A program to determine the radiation spectra due to a single atomic 
  // subshell ionisation by a particle or due to deexcitation or decay of radionuclides", 
  // Comp. Phys. Comm. 1206 pp 1-19 (1997)
  // 

  aParticleChange.Initialize(trackData);
  
  G4Material*               aMaterial = trackData.GetMaterial() ;
  const G4DynamicParticle*  aParticle = trackData.GetDynamicParticle() ;

  //  get kinetic energy cut for the electron
  G4double* DeltaCutInKineticEnergy = G4Electron::Electron()->GetCutsInEnergy() ;
  G4double  DeltaThreshold = DeltaCutInKineticEnergy[aMaterial->GetIndex()];
 
  G4double KineticEnergy = aParticle->GetKineticEnergy();

  // kill e- if Tkin is too small ......
  if(KineticEnergy <= 2.*DeltaThreshold)
  { 
    aParticleChange.SetStatusChange(fStopAndKill);
    aParticleChange.SetEnergyChange(0.);
    
 if(KineticEnergy < 0.)
   G4cout << " 1. negative deposit:" << KineticEnergy/eV << endl; 
    aParticleChange.SetLocalEnergyDeposit(KineticEnergy);

    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  }

  // select randomly one element constituing the material.
  G4Element* anElement = SelectRandomAtom(aParticle, aMaterial);
  G4int AtomIndex = (G4int) anElement->GetZ();

  // First Ionised subshell chosen basing on subshell integrated cross section EPDL97
  // partial sum method
  // Select the subshell WARNING!!!!: it returns the subshell index in the table.
  G4int subShellIndex = SelectRandomShell(AtomIndex, KineticEnergy);

  //mv ionsubShellIndex = 20; // WARNING: REMEMBER TO REMOVE THIS CONSTRAIN
  G4FirstLevel* theBindEnVec = (*theBindingEnergyTable)[AtomIndex-1];
  G4int thePrimaryShell = (G4int) (*(*theBindEnVec)[0])[subShellIndex];
  G4double BindingEn = (*(*theBindEnVec)[1])[subShellIndex];
  G4double theEnergyDeposit = BindingEn; 

  if(KineticEnergy <= BindingEn)
  {
    G4cout << " Tkin=" << KineticEnergy/eV << "  Ebind=" << BindingEn/eV
           << "   selection of subshell ???????" << endl;
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }

  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection();

  // some kinematics

  G4double MaxKineticEnergyTransfer = 0.5*KineticEnergy;

  // sampling kinetic energy of the delta ray 
  if (MaxKineticEnergyTransfer <= LowestKineticEnergy/2){

    // pathological case (should not happen, there is no change at all)
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }

  // **** normal case ****

  // Energy Sampling
  G4double DeltaKineticEnergy = EnergySampling(AtomIndex, subShellIndex, KineticEnergy);
  G4double finalKineticEnergy = KineticEnergy - DeltaKineticEnergy - BindingEn;

  if(finalKineticEnergy < 0.){

    G4cout << "Tkin=" << KineticEnergy/eV << "  Tdel=" << DeltaKineticEnergy/eV
	   << "  BindingEn=" << BindingEn/eV << "  ***********" << endl;
  }

  // deposit energy if delta energy is below cut
  if(DeltaKineticEnergy <= DeltaThreshold){

    if(finalKineticEnergy > 2.*DeltaThreshold){

      aParticleChange.SetEnergyChange(KineticEnergy - DeltaKineticEnergy - BindingEn);

      if((DeltaKineticEnergy+BindingEn) < 0.){

	G4cout << " 2. negative deposit:" << (DeltaKineticEnergy+BindingEn)/eV << endl; 
      }

      aParticleChange.SetLocalEnergyDeposit(DeltaKineticEnergy+BindingEn);
    }
    else{
      
      aParticleChange.SetStatusChange(fStopAndKill);
      aParticleChange.SetEnergyChange(0.);

      if(KineticEnergy < 0.){

	G4cout << " 3. negative deposit:" << KineticEnergy/eV << endl; 
      }

      aParticleChange.SetLocalEnergyDeposit(KineticEnergy);
  
    }

    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  } 
    
  if(thePrimShVec.length() != 0){
    
    thePrimShVec.clear();
  }
  
  thePrimShVec.insert(thePrimaryShell);
  
  // delta ray kinematics
  // step 1. costheta from assumption : atomic e- free and at rest 

    G4double costheta = sqrt(DeltaKineticEnergy*(KineticEnergy+2.*electron_mass_c2)/
                             KineticEnergy*(DeltaKineticEnergy+2.*electron_mass_c2));

    if (costheta < -1.) costheta = -1.;
    if (costheta > +1.) costheta = +1.;
    
    G4double phi = twopi * G4UniformRand(); 
    G4double sintheta = sqrt((1.+costheta)*(1.-costheta));
    
   // step 2. correction of the delta momentum direction 
    G4double dirx,diry,dirz ;
    static G4double cmom=2. ;

    G4double peinit=sqrt(2.*cmom*electron_mass_c2*BindingEn) ;
    G4double pdelta=sqrt(DeltaKineticEnergy*(DeltaKineticEnergy+2.*electron_mass_c2)) ;
    G4double ct=-1.*2.*G4UniformRand() ;
    if(ct > 1.) ct=1. ;
    if(ct <-1.) ct=-1. ;
    G4double st=sqrt(1.-ct*ct) ;
    G4double ph = twopi * G4UniformRand();
    dirx=pdelta*sintheta*cos(phi)+peinit*st*cos(ph) ;
    diry=pdelta*sintheta*sin(phi)+peinit*st*sin(ph) ;
    dirz=pdelta*costheta+peinit*ct ;
    G4double mom=sqrt(dirx*dirx+diry*diry+dirz*dirz) ;
    if(mom > 0.)
    {
     dirx=dirx/mom ;
     diry=diry/mom ;
     dirz=dirz/mom ;
    }
    else
    {
     dirx = sintheta * cos(phi), diry = sintheta * sin(phi), dirz = costheta;
    }

    G4ThreeVector DeltaDirection(dirx,diry,dirz);
    DeltaDirection.rotateUz(ParticleDirection);
    
    G4double finalPx=ParticleDirection.x() ;
    G4double finalPy=ParticleDirection.y() ;
    G4double finalPz=ParticleDirection.z() ;
    
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
    
    if(finalKineticEnergy > 2.*DeltaThreshold)
    {
      aParticleChange.SetMomentumChange(finalPx,finalPy,finalPz);
      aParticleChange.SetEnergyChange(finalKineticEnergy);
 if(theEnergyDeposit < 0.)
   G4cout << " 4. negative deposit:" << theEnergyDeposit/eV << endl; 
      aParticleChange.SetLocalEnergyDeposit (theEnergyDeposit);
    }
    else
    {
      theEnergyDeposit += finalKineticEnergy ;
 if(theEnergyDeposit < 0.)
 {
   G4cout << " 5. negative deposit:" << theEnergyDeposit/eV << endl; 
   G4cout << " finalKineticEnergy=" << finalKineticEnergy/eV << endl;
 }

      aParticleChange.SetLocalEnergyDeposit (theEnergyDeposit);
      aParticleChange.SetEnergyChange(0.);
      aParticleChange.SetStatusChange(fStopAndKill);
    }
    
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//
// selection of the first ionized subshell: partial sum method
// 
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

//
// Element selction in the material already in the standard processes
// 
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

//
// Select a random transition with the transition probabilities and the partial sum 
// method using EADL data (A. Forti)
//

G4bool G4LowEnergyIonisation::SelectRandomTransition(G4int thePrimShell, 
						     G4double* TransParam,
						     const oneAtomTable* TransitionTable){
  
  G4int SubShellCol = 0, ProbCol = 1, EnergyCol = 2;
  // transitionTable contains all the transition probabilities of one atom:
  // loop on subshell is inside the method.

  // when the last subshell is reached CollIsFull becomes FALSE.
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

    // TransProb is the index of the loop and of the table of transition. it starts from 1 
    // because the first element of the data table is the primary shell id number and not a 
    // transition probability: it must not be added to TotalSum. 
  
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

  // Sampling Method: acceptation - rejection method from
  // EGS4 W.R. Nelson et al. The EGS4 Code System. SLAC-Report-265 , December 1985 
  // Fit functions:
  //
  // first function: sum_i A_i/(en+bind)^i  i=2,7 used in all the cases/subshells (see below)
  //
  // second function: B_1*(en)**B_2 used in 9 cases/subshells (no references)
  //
  // third function: sum_i C_i/(en+bind)^i i=2,5 used in 462 cases/subshells (see below)
  //
  // Y.-K.Kim M.E.Rudd http://physics.nist.gov/PhysRefData/Ionization/intro.html
  // First and Third functions are a modified version: degree of monomial in the 
  // truncated sum is from 2 to 7 instead that from 2 to 6. 
  // method A. Forti.

  //  1) Load Coefficients (I need Z number and the index of the shell)
  oneAtomTable* oneAtomCoeffTable = (*theSamplingCoeffTable)[ZNumVec->index(AtomicNumber)];

  oneShellTable* oneShellCoeffTable = (*oneAtomCoeffTable)[ShellIndex];
  G4double BindingEn = (*(*(*theBindingEnergyTable)[AtomicNumber-1])[1])[ShellIndex];

  // 2) Interpolate coefficients (I need the incoming electron kinetic energy)
  const G4int CoeffNumber = oneShellCoeffTable->entries();

  const G4Data* energyVec = (*oneShellCoeffTable)[0];
  const G4int LastPar = energyVec->length()-1;
  G4Data Parms;

  for(G4int ind = 1; ind < CoeffNumber-2; ind++){

    const G4Data* oneCoeffVec = (*oneShellCoeffTable)[ind]; 

    if(KinEn < (*energyVec)[0]){
      Parms.insert((*oneCoeffVec)[0]);
    }

    else if(KinEn > (*energyVec)[LastPar]){

      Parms.insert((*oneCoeffVec)[LastPar]);
    }

    else{

      // NO REFERENCES for SEMI-LOG interpolation it is used because coefficients
      // can be negative. 

      G4double par = util.DataSemiLogInterpolation(KinEn,(*energyVec),(*oneCoeffVec));
      Parms.insert(par);
    }
  }

  // cuts in energy are always the same

  // First cut
  G4double fstCut = (*(*oneShellCoeffTable)[CoeffNumber-2])[0] == 0. ? DBL_MAX : (*(*oneShellCoeffTable)[CoeffNumber-2])[0]*MeV;

  // Second cut
  G4double sndCut = (*(*oneShellCoeffTable)[CoeffNumber-1])[0] == 0. ? DBL_MAX : (*(*oneShellCoeffTable)[CoeffNumber-1])[0]*MeV;

  // 2') order of parameters:
  //     * Parms[0] = par1 LET
  //     * Parms[1] = par2 LET
  //     * Parms[2] = par3 LET
  //     * Parms[3] = par4 LET
  //     * Parms[4] = par5 LET
  //     * Parms[5] = par6 LET
  //     * Parms[6] = max rej 1
  //     * Parms[7] = b1 HET      
  //     * Parms[8] = b2 HET
  //     * Parms[9] = b3 HET
  //     * Parms[10] = b4 HET
  //     * Parms[11] = max rej 2
  //     * Parms[12] = c1 MET
  //     * Parms[13] = c2 MET

  // 3) Compute partial areas (with functions here the cut is used)
  // **************************************************************
  // Sampling energy limits and areas variable declaration
  static const G4double minEn = 0.1*eV;
  G4double maxEn = (KinEn-BindingEn)/2;
  G4double area1 = 0.,  area2 = 0., area3 = 0.;

  // **************************************************************
  // AREAS
  // **************************************************************
  // compute here the weights of the two regions (areas)
  // temporary solution, but needed for consistent sampling
  
  // NOTATION the number order of the functions follows the energy increasing order
  // even if what is called second function should be the third following the number 
  // of cases in which it appears (only 9)

  // Area of the first function ALWAYS computed
  G4double high = maxEn+BindingEn;
  if(fstCut < maxEn){
    
    high = fstCut+BindingEn;
  }
  
  G4double low = BindingEn+minEn;

  //
  // First Function: sum_2(A_i/(energy+Binding)^i)  i = 2,3,4,5,6,7
  //
  G4double aa=1./low ;
  G4double bb=1./high ;
  G4double saa = aa, sbb = bb;
  G4double llow = low*low;
  G4double s1 = 0. ;
  
  for (G4int ii = 1; ii < 7; ii++){

    //
    // area1: integral of the normalized first function
    //
    area1 += Parms[ii-1]*(saa-sbb)/float(ii) ;
    saa *= aa ;
    sbb *= bb ;
    
    //
    //function itself at the minimum value (0.1*eV)
    //
    s1 += Parms[ii-1]/llow;
    llow *= low ;
  }
  
  // Area of the second function of the third function
  //
  // ATTENTION!!!! this flux works because in the three function case maxEn 
  // is always > than the second cut 
  G4double whatCut = 0.;
  if (maxEn >= sndCut){
    
    //
    // three functions used. 
    //

    // Second Function: B1*energy**B2
    //
    G4double expon = Parms[13]-1;
    
    // area2: integral of the normalized second function
    area2 = (Parms[12]/expon)*(pow(sndCut,expon)-pow(fstCut,expon));
    
    // Third Function: sum_2(Ci/(energy+Binding)^i)  i = 2,3,4,5 
    // it is the same as the first function but with two degree less
    //
    aa = 1./(sndCut+BindingEn); // equivalent to 1./low in the first function 
    bb = 1./(maxEn+BindingEn);  // equivalent to 1./high in the first function
    
    G4double saa = aa, sbb = bb;
    
    for(G4int kk = 1; kk < 5; kk++){ 
      
      // area3: integral of the normalized third function
      area3 += Parms[kk+6]*(saa-sbb)/float(kk);
      saa *= aa;
      sbb *= bb;     
    }
    
    whatCut = sndCut;
  }
  
  else if(maxEn >= fstCut){ 
    
    // two functions used.
    //
    // In this case the area2 keeps her initial value area2 = 0.
    //
    // Third Function: sum_2(Ci/(energy+Binding)^i)  i = 2,3,4,5 
    // it is the same as the first function but with two degree less
    //
    aa = 1./(fstCut+BindingEn); // equivalent to 1./low in the first function 
    bb = 1./(maxEn+BindingEn);  // equivalent to 1./high in the first function
    G4double saa = aa, sbb = bb;
    
    for(G4int kk = 1; kk < 5; kk++){

      // area3: integral of the normalized third function
      area3 += Parms[kk+6]*(saa-sbb)/float(kk);
      saa *= aa;
      sbb *= bb;     
    }
     
    whatCut = fstCut;
  }
   
  // **********************************************************************************
  // SAMPLING 
  // **********************************************************************************
  
  G4double areaTot = area1+area2+area3;
  G4double areaDue = area1+area2;
  G4int which;
  // 4) Generate a random number .to select the region of work
  G4double rand1 = areaTot*G4UniformRand();
  
  // 5) Sampling
  G4double sample = 0;
  
  if(rand1 <= area1){
    // Sampling from the first function
    
    // Low energy transfer 
    G4double rejection = 0;
    which = 1;
    do{
      
      G4double rand2 = G4UniformRand();
      G4double Ka = 0;
      
      if(fstCut < maxEn){
	
	Ka = (BindingEn + fstCut)/(minEn+BindingEn);
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
  
  else if(rand1 > area1 && rand1 <= areaDue){
    
    //Sampling from the second function only 9 subshells
    G4double expon = Parms[13]-1;
    G4double norm = (pow(sndCut,expon)-pow(fstCut,expon));
    G4double sum = norm*G4UniformRand()+pow(fstCut,expon);
    G4double exponInv = 1/expon;
    sample = pow(sum,exponInv);
    which = 2;
  }
  
  else if (rand1 > areaDue && rand1 <= areaTot){ 
    
    // Sampling from the third function 462 cases of which 9 with the second
    G4double rejection = 0;
    which = 3;
    do{
      
      G4double rand2 = G4UniformRand();
      G4double Ka = 0;
      
      if(whatCut < maxEn){
	
	Ka = (BindingEn + maxEn)/(whatCut+BindingEn);
      }
      
      sample = (whatCut + BindingEn)*pow(Ka,rand2)-BindingEn;
      
      G4double arg = sample + BindingEn;
      
      rejection = Parms[7]/arg+Parms[8]/pow(arg,2)+Parms[9]/pow(arg,3)+Parms[10]/pow(arg,4);
      
      rejection /= Parms[11];
      
    }while(rejection < G4UniformRand()); 
  }
  
  return sample;
}

void G4LowEnergyIonisation::PrintInfoDefinition()
{
  G4String comments = "First version of low energy ionisation code,";
           comments += "\n At present it can be used for electrons only ";
           comments += " in the energy range [250 eV,100 GeV]";
  G4cout << endl << GetProcessName() << ":  " << comments << endl;

}





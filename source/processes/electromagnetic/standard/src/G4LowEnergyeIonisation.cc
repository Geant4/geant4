// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyeIonisation.cc,v 1.1 1999-01-08 14:16:38 gunter Exp $
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
//      ---------- G4LowEnergyeIonisation physics process -----------
//                by Laszlo Urban, 20 March 1997 
// **************************************************************
// It is the first implementation of the NEW IONISATION PROCESS.
// It calculates the ionisation of e+/e-.
// **************************************************************
//
// 07-04-98: remove 'tracking cut' of the ionizing particle, MMa 
// --------------------------------------------------------------
 

#include "G4LowEnergyeIonisation.hh"
#include "G4EpdlData.hh"
#include "G4PhysicsTable.hh"
#include "CLHEP/String/Strings.h"
// constructor and destructor
 
G4LowEnergyeIonisation::G4LowEnergyeIonisation(const G4String& processName)
   : G4eIonisation("eIonisation"),
     LowestKineticEnergy(1*keV),
     HighestKineticEnergy(10.*keV),
     TotBin(100)
{
  theCrossSectionTable = 0;
  theMeanFreePathTable  = 0;
  theBindingEnergyTable = 0;
  theRecoilElectronSpectrumTable = 0;
  theRecoilElectronEnergiesTable = 0;

}

G4LowEnergyeIonisation::~G4LowEnergyeIonisation() 
{
  if(theCrossSectionTable){
  
    theCrossSectionTable->clearAndDestroy(); 
    delete theCrossSectionTable;
  }

  if(theBindingEnergyTable){

    theBindingEnergyTable->clearAndDestroy(); 
    delete theBindingEnergyTable;
  }

  if(theRecoilElectronSpectrumTable){

    theRecoilElectronSpectrumTable->clearAndDestroy(); 
    delete theRecoilElectronSpectrumTable;
  }

  if(theRecoilElectronEnergiesTable){

    theRecoilElectronEnergiesTable->clearAndDestroy(); 
    delete theRecoilElectronEnergiesTable;
  }

  if(theMeanFreePathTable){
    
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }

}
 
 
// methods.............................................

void G4LowEnergyeIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
//  build Cross Section Table, Recoil Electron Tables, call BuildLossTable+BuildLambdaTable
{
 
  // Build Cross Sections Table
  BuildCrossSectionTable(aParticleType);

  BuildLossTable(aParticleType) ;

  if(&aParticleType == theElectron)
  {
    RecorderOfElectronProcess[CounterOfElectronProcess] = (*this).theLossTable ;
    CounterOfElectronProcess++;
  }
  else
  {
    RecorderOfPositronProcess[CounterOfPositronProcess] = (*this).theLossTable ;
    CounterOfPositronProcess++;
  }
 
    BuildLambdaTable(aParticleType) ;
 
    BuildDEDXTable(aParticleType) ;
  
}

void G4LowEnergyeIonisation::BuildCrossSectionTable(const G4ParticleDefinition& aParticleType){

  // Subshells Cross Section Tables
  HepString file;
  G4double param[4] = {81, 21, 91, 0}; // check 21 number!!!!!!!!!!
  if(&aParticleType == theElectron) {
    file = "eedl.asc";
  }

  else {
    file = "epodl.asc";
  }

  //At the moment only for the K subshell
  for(G4int subShInd = 1; subShInd < 62; subShInd++){

    param[3] = subShInd;
    G4EpdlData c(file, param);
    
    c.FillDataTable(1., 1.);
    if(subShInd == 1) G4PhysicsTable* theCrossSectionTable = new G4PhysicsTable((*c.GetFstDataTable()));
    if(subShInd >1){ 

      // WARNING! no tables for other subshells

    }
  }
}

void G4LowEnergyeIonisation::BuildBindingEnergiesTable(const G4ParticleDefinition& aParticleType){

  // Build the binding energy table
  if(theBindingEnergyTable){
    theBindingEnergyTable->clearAndDestroy(); delete theBindingEnergyTable;
  }
  
  G4double param1[4] = {91, 913, 0, 0};
  HepString file1("eadl.asc");
  G4EpdlData c1(file1, param1);
  c1.FillDataTable( 1., 1.);
  theBindingEnergyTable = new G4PhysicsTable(*(c1.GetFstDataTable())) ;
}

void G4LowEnergyeIonisation::BuildRecoilElectronTables(const G4ParticleDefinition& aParticleType){

  HepString file;
  G4double param[4] = {81, 22, 91, 0};

  if(&aParticleType == theElectron) {
    file = "eedl.asc";
  }

  else {
    file = "epodl.asc";
  }

  G4PhysicsTable* tmpReST;
  G4PhysicsTable* tmpReET;
  
  // at the moment only for K subshell
  for(G4int subShInd = 1; subShInd < 62; subShInd++){
    
    param[7] = subShInd;
    G4EpdlData c(file, param);
    c.FillDataTable(1., 1.); // 
    
    if(subShInd == 1) {
      
      tmpReST = c.GetFstDataTable();
      tmpReET = c.GetSndDataTable();
    }
    else if(subShInd >1){ 
      // WARNING! no tables for other subshells
    }
  }

   theRecoilElectronSpectrumTable = tmpReST;
  theRecoilElectronEnergiesTable = tmpReET;
}

void G4LowEnergyeIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  // Build tables for the ionization energy loss
  //  the tables are built for MATERIALS
  //                           *********
  
  G4double LowEdgeEnergy , ionloss ;
  G4bool isOutRange ;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4double twoln10 = 2.*log(10.) ;
  const G4double Factor = twopi_mc2_rcl2 ;
  
  ParticleMass = aParticleType.GetPDGMass();
  
  ParticleCutInKineticEnergy = aParticleType.GetEnergyCuts() ;
  
  //  create table
  
  G4int numOfMaterials = theMaterialTable->length();
  
  if (theLossTable) {
    theLossTable->clearAndDestroy();
    delete theLossTable;
  }
  theLossTable = new G4PhysicsTable(numOfMaterials);
  
  //  loop for materials
  
  for (G4int J=0; J<numOfMaterials; J++){

    // create physics vector and fill it
    
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy, HighestKineticEnergy, TotBin);
  
    // get material parameters needed for the energy loss calculation
    
    G4double ElectronDensity,Eexc,Eexcm2,Cden,Mden,Aden,X0den,X1den ;
    
    const G4Material* material= (*theMaterialTable)[J];
    
    ElectronDensity = material->GetElectronDensity();
    Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
    Eexcm2 = Eexc/ParticleMass ;
    Eexcm2 *= Eexcm2 ;
    Cden = material->GetIonisation()->GetCdensity();
    Mden = material->GetIonisation()->GetMdensity();
    Aden = material->GetIonisation()->GetAdensity();
    X0den = material->GetIonisation()->GetX0density();
    X1den = material->GetIonisation()->GetX1density();
    
    // get particle cut in kin. energy for the material
    
    ParticleCutInKineticEnergyNow = ParticleCutInKineticEnergy[J] ;
    
    // some local variables -------------------
    G4double tau,Tmax,gamma,gamma2,bg2,beta2,d,d2,d3,d4,delta,x,y ;
    
    // now comes the loop for the kinetic energy values*****************
    
    for (G4int i = 0 ; i < TotBin ; i++) {
      
      LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
      
      tau = LowEdgeEnergy/ParticleMass ;
      
      // Seltzer-Berger formula 
      
      gamma = tau +1. ;
      bg2 = tau*(tau+2.) ;
      gamma2 = gamma*gamma ;
      beta2 = bg2/gamma2 ;
      
      // electron .................................
      if(&aParticleType==theElectron){
	
	Tmax = LowEdgeEnergy/2. ;
	
	d = min(ParticleCutInKineticEnergyNow, Tmax)/ParticleMass;
	
	ionloss = log(2.*(tau+2.)/Eexcm2)-1.-beta2 ;
	ionloss += log((tau-d)*d)+tau/(tau-d) ;
	ionloss += (0.5*d*d+(2.*tau+1.)*log(1.-d/tau))/gamma2 ;
      }
      //positron ...............................
      else{
	
	Tmax = LowEdgeEnergy ;
	
	d = min(ParticleCutInKineticEnergyNow, Tmax)/ParticleMass;
	    
	d2=d*d/2. ;
	d3=d*d*d/3. ;
	d4=d*d*d*d/4. ;
	
	y=1./(1.+gamma) ;
	
	ionloss = log(2.*(tau+2.)/Eexcm2)+log(tau*d) ;
	ionloss-= beta2*(tau+2.*d-y*(3.*d2+y*(d-d3+y*(d2-tau*d3+d4))))/tau;
	
      } 

      // density correction   ................................
      
      x = log(bg2)/twoln10 ;
      
      if ( x < X0den )
	delta = 0. ;
      else{
	delta = twoln10*x - Cden ;
	
	if ( x < X1den )
	  delta += Aden*pow((X1den-x),Mden) ;
      } 
      
      // now you can compute the total ionization loss
      
      ionloss -= delta ;
      
      ionloss *= Factor*ElectronDensity/beta2 ;
      
      if ( ionloss <= 0.)
	ionloss = 0. ;
      
      aVector->PutValue(i,ionloss) ;
	  
    }
          
    theLossTable->insert(aVector);
  }
}

void G4LowEnergyeIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{
  // Build mean free path tables for the delta ray production process
  //     tables are built for MATERIALS 

      G4double Value, sigma ;
      G4bool isOutRange ;

      //create table
      const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
      G4int numOfMaterials = theMaterialTable->length();

      if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
      }

      theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);

      // get electron and particle cuts in kinetic energy
      // DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
      // ParticleCutInKineticEnergy = aParticleType.GetEnergyCuts() ;

      // loop for materials 
      for (G4int J=0 ; J < numOfMaterials; J++){ // for each material

        //create physics vector then fill it ....
        G4PhysicsFreeVector* aVector = new G4PhysicsFreeVector(TotBin);
	
       // compute the (macroscopic) cross section first
         const G4Material* material= (*theMaterialTable)[J];
        
        const G4ElementVector* theElementVector = material->GetElementVector() ;
        const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
        const G4int NumberOfElements = material->GetNumberOfElements() ;
 
       // get the electron kinetic energy cut for the actual material,
       //  it will be used in ComputeMicroscopicCrossSection
       // ( it is the SAME for ALL the ELEMENTS in THIS MATERIAL )
       //   ------------------------------------------------------

       //        DeltaKineticEnergyCutNow = DeltaCutInKineticEnergy[J] ;

        for ( G4int i = 0 ; i < TotBin ; i++ ){ // for each energy
	  
	  sigma = 0. ;
	  
	  for (G4int iel=0; iel<NumberOfElements; iel++ ){ // for each element
	    
	    sigma +=  theAtomicNumDensityVector[iel]*
	      (*(*theCrossSectionTable)((*theElementVector)(iel)->GetZ()))(i); 
	  }
	  
	  // mean free path = 1./macroscopic cross section
	  
	  Value = sigma<=0 ? BIGSTEP : 1./sigma ;     
        // MUST BE ABSOLUTELY CORRECT BUT IT DEPENDS ON HOW I  CONSTRUCT THE TABLES
       // THE ENERGIES OF EACH ELEMENT PHYSICS VECTOR ARE NOT THE SAME SO AT THE MOMENT
       // I TAKE THE FIRST ELEMENT ONES the energies must be uniform I don't know how
       // NumbBinTable must be the same for all the elements
	  G4double theEnergy = (*theCrossSectionTable)(0)->GetLowEdgeEnergy(i);
	  aVector->PutValue(i, theEnergy,Value) ;
        }

        theMeanFreePathTable->insert(aVector);
      }
}

G4VParticleChange* G4LowEnergyeIonisation::PostStepDoIt(
                                              const G4Track& trackData,   
                                              const G4Step& stepData)         
{
  // Units are expressed in GEANT4 internal units.

  aParticleChange.Initialize(trackData) ;

  const G4DynamicParticle* aParticle ;
  G4Material* aMaterial;  
  aParticle = trackData.GetDynamicParticle() ;
  aMaterial = trackData.GetMaterial() ;

  // Incoming electron variables
  G4double Charge, KineticEnergy, TotalEnergy, Psquare, Esquare; // incoming electron
  Charge=aParticle->GetDefinition()->GetPDGCharge();
  KineticEnergy=aParticle->GetKineticEnergy();
  TotalEnergy=KineticEnergy + ParticleMass ;
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
  Esquare=TotalEnergy*TotalEnergy ;
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection() ;

  //  get kinetic energy cut for the electron....
  DeltaKineticEnergyCutNow = DeltaCutInKineticEnergy[aMaterial->GetIndex()];

  // some kinematics......................
  G4double MaxKineticEnergyTransfer;
  if(Charge<0.)
   {
     MaxKineticEnergyTransfer = 0.5*KineticEnergy ;
     ParticleCutInKineticEnergyNow = ((*theElectron).GetCutsInEnergy())[aMaterial->GetIndex()];
   }
  else
   {
     MaxKineticEnergyTransfer =     KineticEnergy ;
     ParticleCutInKineticEnergyNow = ((*thePositron).GetCutsInEnergy())[aMaterial->GetIndex()];
   }

  // Delta ray (recoil electron)  and scattered electron energies
  G4double DeltaKineticEnergy,  ScatteredKineticEnergy;

  // sampling kinetic energy of the delta ray 

  if( MaxKineticEnergyTransfer <= DeltaKineticEnergyCutNow )
  {
    // pathological case (it should not happen ,
    // there is no change at all).....

     //return &aParticleChange;
     return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  else{
    
    // normal case ......................................
    
    //
    // FROM HERE I HAVE TO PUT STEPANEK ALGORITHM
    //

    SpectrumFit();
    //
    // Spectrum fit: ??????????????????????????????????????????? needed to 
    // choose which distribution has to be used to sample the recoil electron energy ATTENZIONE !!!!!!! siamo in un materiale
    // E' UN CASINO PERCHE' PER OGNI ELEMENTO DI UN MATERIALE BISOGNA CALCOLARE LA DISTIRBUZIONE DI 
    // PROBABILITA' DELLO SPETTRO PER OGNI SHELL, PER OGNI ENERGIA INCIDENTE, PER OGNI ENERGIE DEL RAGGIO DELTA.............
    // E NON SONO SICURA DI AVERE CAPITO COME SI FA.
    //


    G4double p1 = 0, p2 = 0, p3 = 0, BindingEnergy;
    
    // After that
    G4double rand1 = G4UniformRand();
    G4double rand2 = G4UniformRand();

    G4double kinMinusBin = KineticEnergy -BindingEnergy;
    G4double kinPlusBin = KineticEnergy -BindingEnergy;
    
    if( rand1>= 0 && rand1<= p1)
      DeltaKineticEnergy =  (kinMinusBin*BindingEnergy*(1-rand2))/(rand2*kinMinusBin +2*BindingEnergy);

    else if( rand1 > p1 && rand1 <= (p1 +p2))
      DeltaKineticEnergy = (KineticEnergy*kinMinusBin*(1 - rand2))/(2*KineticEnergy - rand2*kinMinusBin);
    
    else if( rand1 > ( p1+p2 ) && rand1 <= 1){

      G4double var = exp( rand2*log((4*BindingEnergy*KineticEnergy)/(kinPlusBin*kinPlusBin)) +log((kinPlusBin*kinPlusBin)/4));
      DeltaKineticEnergy = 0.5*( kinMinusBin - sqrt(kinMinusBin*kinMinusBin + 4*(BindingEnergy*KineticEnergy - var)));
    }

    // Non c'e' nessuna funzione di rejezione ?????????????????????

    ScatteredKineticEnergy = kinMinusBin - DeltaKineticEnergy;
  } // end if(E > ...) else

    // protection :DO NOT PRODUCE a secondary with 0. kinetic energy !!!!!!
    // ********************************************************************
    if( DeltaKineticEnergy <= 0.)
      return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

    G4double ScatteredTotalMomentum, DeltaTotalMomentum, TotalMomentum;
    ScatteredTotalMomentum = sqrt(ScatteredKineticEnergy * (ScatteredKineticEnergy + 2. * electron_mass_c2 )) ; // momentum of the scattered electron
    TotalMomentum = sqrt(Psquare) ; // momentum of the primary electron
    DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy + 2. * electron_mass_c2 )) ; // momentum of the scattered electron

   //   costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2) 
   //       /(DeltaTotalMomentum * TotalMomentum) ;

    G4double cosThetaSca, cosThetaDel;
   cosThetaSca =  (TotalMomentum*TotalMomentum + ScatteredTotalMomentum*ScatteredTotalMomentum - DeltaTotalMomentum*DeltaTotalMomentum)/(2*TotalMomentum*ScatteredTotalMomentum);
 
   cosThetaDel =  (TotalMomentum*TotalMomentum - ScatteredTotalMomentum*ScatteredTotalMomentum + DeltaTotalMomentum*DeltaTotalMomentum)/(2*TotalMomentum*DeltaTotalMomentum);
   //  protection against costheta > 1 or < -1   --------------

   if ( cosThetaDel < -1. ) 
          cosThetaDel = -1. ;
   if ( cosThetaDel > +1. ) 
          cosThetaDel = +1. ;

   //  direction of the scattered electron  ........

   G4double phiDel, sinThetaDel, dirx, diry, dirz;
   phiDel = twopi * G4UniformRand() ; 
   sinThetaDel = sqrt(1.-cosThetaDel*cosThetaDel);
   dirx = sinThetaDel*cos(phiDel) ;
   diry = sinThetaDel*sin(phiDel) ;
   dirz = cosThetaDel ;

   G4ThreeVector DeltaDirection(dirx,diry,dirz) ;
   DeltaDirection.rotateUz(ParticleDirection) ;

   // create G4DynamicParticle object for delta ray

   G4DynamicParticle* theDeltaRay = new G4DynamicParticle;
   theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
   theDeltaRay->SetMomentumDirection(DeltaDirection.x(), DeltaDirection.y(), DeltaDirection.z()); 
   theDeltaRay->SetDefinition(G4Electron::Electron());
   
   // fill aParticleChange 
   // changed energy and momentum of the actual particle
   G4double ScatteredPx, ScatteredPy, ScatteredPz;
   if (ScatteredKineticEnergy > 0.){
     
     ScatteredPx = (TotalMomentum*ParticleDirection.x() - DeltaTotalMomentum*DeltaDirection.x())/ScatteredTotalMomentum ; 
     ScatteredPy = (TotalMomentum*ParticleDirection.y() - DeltaTotalMomentum*DeltaDirection.y())/ScatteredTotalMomentum ; 
     ScatteredPz = (TotalMomentum*ParticleDirection.z() - DeltaTotalMomentum*DeltaDirection.z())/ScatteredTotalMomentum ; 

     aParticleChange.SetMomentumChange( ScatteredPx, ScatteredPy, ScatteredPz );
   }

   else{
        
     ScatteredKineticEnergy = 0.;
     if (Charge < 0.) aParticleChange.SetStatusChange(fStopAndKill);
     else  aParticleChange.SetStatusChange(fStopButAlive);
   }

   aParticleChange.SetEnergyChange( ScatteredKineticEnergy );
   aParticleChange.SetNumberOfSecondaries(1);  
   aParticleChange.AddSecondary( theDeltaRay );
   aParticleChange.SetLocalEnergyDeposit (0.);
      
   return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

void G4LowEnergyeIonisation::SpectrumFit(){

  cout<<"SpectrumFit not available yet"<<endl;
}
  













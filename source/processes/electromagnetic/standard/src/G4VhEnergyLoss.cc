// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VhEnergyLoss.cc,v 1.12 2000-08-18 17:37:30 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// -----------------------------------------------------------
// 7/10/98: bug fixes + some cleanup , L.Urban 
// 22/10/98 : cleanup , L.Urban
// 07/12/98 : works for ions as well+ bug corrected, L.Urban
// 02/02/99 : several bugs fixed, L.Urban
// 01/03/99 : creation of sub-cutoff delta rays, L.Urban
// 28/04/99 : bug fixed in DoIt , L.Urban
// 10/02/00  modifications , new e.m. structure, L.Urban
// 18/07/00 : bug fix in AlongStepDoIt V.Ivanchenko
// 10/08/00 : V.Ivanchenko change AlongStepDoIt and
//            add EnergyLossFluctuation in order to simulate
//            energy losses of ions
// 17/08/00 : V.Ivanchenko change EnergyLossFluctuation 
// 18/08/00 : V.Ivanchenko bug fixed in GetConstrained 
// --------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VhEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4Poisson.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int            G4VhEnergyLoss::NbOfProcesses    = 1 ;

G4int            G4VhEnergyLoss::CounterOfProcess = 0 ;
G4PhysicsTable** G4VhEnergyLoss::RecorderOfProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4VhEnergyLoss::CounterOfpProcess = 0 ;
G4PhysicsTable** G4VhEnergyLoss::RecorderOfpProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4VhEnergyLoss::CounterOfpbarProcess = 0 ;
G4PhysicsTable** G4VhEnergyLoss::RecorderOfpbarProcess =
                                           new G4PhysicsTable*[10] ;

G4PhysicsTable* G4VhEnergyLoss::theDEDXpTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theDEDXpbarTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theRangepTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theRangepbarTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theInverseRangepTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theInverseRangepbarTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theLabTimepTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theLabTimepbarTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theProperTimepTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::theProperTimepbarTable = NULL ;

G4PhysicsTable* G4VhEnergyLoss::thepRangeCoeffATable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::thepRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::thepRangeCoeffCTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::thepbarRangeCoeffATable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::thepbarRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VhEnergyLoss::thepbarRangeCoeffCTable = NULL ;

G4PhysicsTable* G4VhEnergyLoss::theDEDXTable = NULL ;

const G4Proton* G4VhEnergyLoss::theProton=G4Proton::Proton() ;
const G4AntiProton* G4VhEnergyLoss::theAntiProton=G4AntiProton::AntiProton() ;

G4double G4VhEnergyLoss::ptableElectronCutInRange = 0.0*mm ;
G4double G4VhEnergyLoss::pbartableElectronCutInRange = 0.0*mm ;

G4double         G4VhEnergyLoss::MinDeltaCutInRange = 0.1*mm ;
G4double*        G4VhEnergyLoss::MinDeltaEnergy     = NULL   ;
G4bool		 G4VhEnergyLoss::setMinDeltaCutInRange = false ;

G4double         G4VhEnergyLoss::Charge ;   

G4double G4VhEnergyLoss::LowerBoundEloss = 1.*keV ;
G4double G4VhEnergyLoss::UpperBoundEloss = 100.*TeV;
G4int G4VhEnergyLoss::NbinEloss = 100 ;
G4double G4VhEnergyLoss::RTable,G4VhEnergyLoss::LOGRTable;

G4double G4VhEnergyLoss::c0N       = 9.0e-21*MeV*MeV*mm*mm ;
G4double G4VhEnergyLoss::c1N       = 25.0e-21*keV*mm*mm    ;
G4double G4VhEnergyLoss::c2N       = 13.25e-21*keV*mm*mm   ;
G4double G4VhEnergyLoss::c3N       = 0.500e-21*mm*mm       ;
G4int    G4VhEnergyLoss::Ndeltamax = 100                   ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4VhEnergyLoss::G4VhEnergyLoss(const G4String& processName)
   : G4VEnergyLoss (processName),
     theLossTable (NULL),
     MinKineticEnergy(1.*eV), 
     linLossLimit(0.05)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VhEnergyLoss::~G4VhEnergyLoss() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable; theLossTable = NULL;
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4VhEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{
  //  calculate data members LOGRTable,RTable first
  G4double lrate = log(UpperBoundEloss/LowerBoundEloss);
  LOGRTable=lrate/NbinEloss;
  RTable   =exp(LOGRTable);

  // create table if there is no table or there is a new cut value
  G4bool MakeTable = false ;
     
  G4double ElectronCutInRange = G4Electron::Electron()->GetCuts();

  // create/fill proton or antiproton tables depending on the charge 
  Charge = aParticleType.GetPDGCharge()/eplus;
  ParticleMass = aParticleType.GetPDGMass() ;

  if (Charge>0.) {theDEDXTable= theDEDXpTable;}
  else           {theDEDXTable= theDEDXpbarTable;}

  if(
     ((Charge>0.) && ((theDEDXTable==NULL) || 
     (ElectronCutInRange != ptableElectronCutInRange)))
     ||  
     ((Charge<0.) && ((theDEDXTable==NULL) || 
     (ElectronCutInRange != pbartableElectronCutInRange)))
    )
      MakeTable = true ;
  
  const G4MaterialTable* theMaterialTable=
                                   G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

  if( MakeTable )
  {

  // Build energy loss table as a sum of the energy loss due to the
  //              different processes.                                           
    if( Charge >0.)    
    {
      RecorderOfProcess=RecorderOfpProcess;
      CounterOfProcess=CounterOfpProcess;

      if(CounterOfProcess == NbOfProcesses)
      {
        if(theDEDXpTable)
        { theDEDXpTable->clearAndDestroy();
          delete theDEDXpTable; }
        theDEDXpTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXpTable;
        ptableElectronCutInRange = ElectronCutInRange ;
      }
    }
    else
    {
      RecorderOfProcess=RecorderOfpbarProcess;
      CounterOfProcess=CounterOfpbarProcess;

      if(CounterOfProcess == NbOfProcesses)
      {
        if(theDEDXpbarTable)
        { theDEDXpbarTable->clearAndDestroy();
          delete theDEDXpbarTable; }
        theDEDXpbarTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXpbarTable;
        pbartableElectronCutInRange = ElectronCutInRange ;
      }
    }

    if(CounterOfProcess == NbOfProcesses)
    {
      //  loop for materials
      G4double LowEdgeEnergy , Value ;
      G4bool isOutRange ;
      G4PhysicsTable* pointer ;

      for (G4int J=0; J<numOfMaterials; J++)
      { 
        // create physics vector and fill it
        G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowerBoundEloss, UpperBoundEloss, NbinEloss);   

        // loop for the kinetic energy
        for (G4int i=0; i<NbinEloss; i++)
        {
          LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;      
          Value = 0. ;
    
          // loop for the contributing processes
          for (G4int process=0; process < NbOfProcesses; process++)
          {
            pointer= RecorderOfProcess[process];
            Value += (*pointer)[J]->
                               GetValue(LowEdgeEnergy,isOutRange) ;
          }

          aVector->PutValue(i,Value) ; 
        }

        theDEDXTable->insert(aVector) ;
      }
      
      //  reset counter to zero ..................
      if( Charge >0.)    
        CounterOfpProcess=0 ;
      else
        CounterOfpbarProcess=0 ;
      ParticleMass = aParticleType.GetPDGMass() ;

      if(Charge > 0.)
      {
       // Build range table
       theRangepTable = BuildRangeTable(theDEDXpTable,
                        theRangepTable,
                        LowerBoundEloss,UpperBoundEloss,NbinEloss);
       // Build lab/proper time tables
       theLabTimepTable = BuildLabTimeTable(theDEDXpTable,
                         theLabTimepTable,
                         LowerBoundEloss,UpperBoundEloss,NbinEloss);
       theProperTimepTable = BuildProperTimeTable(theDEDXpTable,
                            theProperTimepTable,
                            LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build coeff tables for the energy loss calculation
       thepRangeCoeffATable = BuildRangeCoeffATable(theRangepTable,
                             thepRangeCoeffATable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepRangeCoeffBTable = BuildRangeCoeffBTable(theRangepTable,
                             thepRangeCoeffBTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepRangeCoeffCTable = BuildRangeCoeffCTable(theRangepTable,
                             thepRangeCoeffCTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // invert the range table
       theInverseRangepTable = BuildInverseRangeTable(theRangepTable,
                              thepRangeCoeffATable,
                              thepRangeCoeffBTable,
                              thepRangeCoeffCTable,
                              theInverseRangepTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);
 
      }
      else
      {
        // Build range table
       theRangepbarTable = BuildRangeTable(theDEDXpbarTable,
                        theRangepbarTable,
                        LowerBoundEloss,UpperBoundEloss,NbinEloss);
       // Build lab/proper time tables
       theLabTimepbarTable = BuildLabTimeTable(theDEDXpbarTable,
                         theLabTimepbarTable,
                         LowerBoundEloss,UpperBoundEloss,NbinEloss);
       theProperTimepbarTable = BuildProperTimeTable(theDEDXpbarTable,
                            theProperTimepbarTable,
                            LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build coeff tables for the energy loss calculation
       thepbarRangeCoeffATable = BuildRangeCoeffATable(theRangepbarTable,
                             thepbarRangeCoeffATable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepbarRangeCoeffBTable = BuildRangeCoeffBTable(theRangepbarTable,
                             thepbarRangeCoeffBTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepbarRangeCoeffCTable = BuildRangeCoeffCTable(theRangepbarTable,
                             thepbarRangeCoeffCTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // invert the range table
       theInverseRangepbarTable = BuildInverseRangeTable(theRangepbarTable,
                              thepbarRangeCoeffATable,
                              thepbarRangeCoeffBTable,
                              thepbarRangeCoeffCTable,
                              theInverseRangepbarTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);
 
      }

    }

  }
  
  // make the energy loss and the range table available

  G4EnergyLossTables::Register(&aParticleType,  
    (Charge>0)?
      theDEDXpTable: theDEDXpbarTable,
    (Charge>0)?
      theRangepTable: theRangepbarTable,
    (Charge>0)?
      theInverseRangepTable: theInverseRangepbarTable,
    (Charge>0)?
      theLabTimepTable: theLabTimepbarTable,
    (Charge>0)?
      theProperTimepTable: theProperTimepbarTable,
    LowerBoundEloss, UpperBoundEloss,
    proton_mass_c2/aParticleType.GetPDGMass(),NbinEloss);

   // create array for the min. delta cuts in kinetic energy
   if(!setMinDeltaCutInRange)
     MinDeltaCutInRange = G4Electron::Electron()->GetCuts()/10.;

  if(aParticleType.GetParticleName()=="proton")
  {
    G4cout << G4endl;
    G4cout.precision(5) ;
    G4cout << " hIoni    Minimum Delta cut in range=" << MinDeltaCutInRange/mm
           << "  mm." << G4endl;
    G4cout << G4endl;
    G4cout << "           material        min.delta energy(keV) " << G4endl;
    G4cout << G4endl;
  }

    if(MinDeltaEnergy) delete MinDeltaEnergy ;
    MinDeltaEnergy = new G4double [numOfMaterials] ;
    G4double Tlowerlimit = 1.*keV ;
    for(G4int mat=0; mat<numOfMaterials; mat++)
    {
      MinDeltaEnergy[mat] = G4EnergyLossTables::GetPreciseEnergyFromRange(
                            G4Electron::Electron(),MinDeltaCutInRange,
                                       (*theMaterialTable)(mat)) ;
      if(MinDeltaEnergy[mat]<Tlowerlimit) MinDeltaEnergy[mat]=Tlowerlimit ;

      if(MinDeltaEnergy[mat]>G4Electron::Electron()->GetCutsInEnergy()[mat])
        MinDeltaEnergy[mat]=G4Electron::Electron()->GetCutsInEnergy()[mat] ;

     if(aParticleType.GetParticleName()=="proton")
       G4cout << G4std::setw(20) << (*theMaterialTable)(mat)->GetName()
              << G4std::setw(15) << MinDeltaEnergy[mat]/keV << G4endl;
    }
}
      
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VhEnergyLoss::GetConstraints(const G4DynamicParticle *aParticle,
                                              G4Material *aMaterial)
{
  // returns the Step limit
  // dRoverRange is the max. allowed relative range loss in one step
  // it calculates dEdx and the range as well....

  G4double KineticEnergy,StepLimit;
  G4bool isOut ;

  Charge = aParticle->GetDefinition()->GetPDGCharge()/eplus ;

  KineticEnergy = aParticle->GetKineticEnergy();

  G4double massratio=proton_mass_c2/
           aParticle->GetDefinition()->GetPDGMass() ;

  G4double Tscaled= KineticEnergy*massratio ; 
  G4double ChargeSquare = Charge*Charge ;

     if(Charge>0.)
     {
       fRangeNow = G4EnergyLossTables::GetRange( theProton,
                                            Tscaled,aMaterial) ;
        fdEdx     = G4EnergyLossTables::GetDEDX( theProton,
                                            Tscaled,aMaterial) ;
     }
     else
     {
       fRangeNow = G4EnergyLossTables::GetRange( theAntiProton,
                                             Tscaled,aMaterial) ;
       fdEdx     = G4EnergyLossTables::GetDEDX( theAntiProton,
                                             Tscaled,aMaterial) ;
     }
     fdEdx     *= ChargeSquare ;
     fRangeNow /= (ChargeSquare*massratio) ;

  // compute the (random) Step limit ..............
  if(fRangeNow > finalRange)
  {
    StepLimit = (c1lim*fRangeNow+c2lim+c3lim/fRangeNow) ;

    //  randomise this value
    if(rndmStepFlag) StepLimit = 
                finalRange+(StepLimit-finalRange)*G4UniformRand() ;
    if(StepLimit > fRangeNow) StepLimit = fRangeNow ;
  }
  else StepLimit = fRangeNow ;

  return StepLimit ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VhEnergyLoss::AlongStepDoIt( 
                              const G4Track& trackData,const G4Step& stepData) 
 // compute the energy loss after a step 
{
  const G4DynamicParticle* aParticle;
  G4Material* aMaterial;
  G4double E,finalT,Step,ChargeSquare,MeanLoss ;

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;
  
  // get the actual (true) Step length from stepData 
  Step = stepData.GetStepLength() ;

  aParticle = trackData.GetDynamicParticle() ;
  ChargeSquare = Charge*Charge ;

  G4int index = aMaterial->GetIndex() ;
  E = aParticle->GetKineticEnergy() ;

  if(E < MinKineticEnergy) MeanLoss = E ;
  else
  {
    if(Step >= fRangeNow ) MeanLoss = E ;

    else if(( E > UpperBoundEloss)||( E <= LowerBoundEloss))
              MeanLoss = Step*fdEdx ;
     
    else
    {
      if(Step>linLossLimit*fRangeNow)
      {
        G4double massratio=proton_mass_c2/
                 aParticle->GetDefinition()->GetPDGMass() ;

        G4double rscaled= fRangeNow*massratio*ChargeSquare ;
        G4double sscaled=   Step   *massratio*ChargeSquare ;

        if(Charge>0.)
        {
          MeanLoss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theProton,
                                         rscaled        ,aMaterial) -
                     G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theProton,
                                         rscaled-sscaled,aMaterial) ;
        }
        else
        {
          MeanLoss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theAntiProton,
                                         rscaled        ,aMaterial) -
                     G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theAntiProton,
                                         rscaled-sscaled,aMaterial) ;
        }
        MeanLoss /= massratio ;
      }
      else MeanLoss = Step*fdEdx ;
    }
  } 
  finalT = E - MeanLoss ;

  //   subcutoff delta ray production start                          
 if(subSecFlag)
 {
  G4double MinDeltaEnergyNow,Tc,TmintoProduceDelta,w,ww ;
  G4double rcut,T0,presafety,postsafety,safety,
           delta,fragment,Tmax,mass ;
  G4double frperstep,x1,y1,z1,dx,dy,dz,dTime,time0,DeltaTime;
  G4double epsil = MinKineticEnergy/2. ;

  MinDeltaEnergyNow = MinDeltaEnergy[index] ;
  Tc=G4Electron::Electron()->GetCutsInEnergy()[index];
  const G4ParticleDefinition* aParticleType=aParticle->GetDefinition() ;
  mass=aParticleType->GetPDGMass() ;
  w=mass+electron_mass_c2 ;
  ww=2.*mass-MinDeltaEnergyNow ;
  TmintoProduceDelta=0.5*(sqrt(ww*ww+2.*w*w*MinDeltaEnergyNow/
                       electron_mass_c2)-ww) ;

  if((E > TmintoProduceDelta) && (MeanLoss > MinDeltaEnergyNow)
                                   && (finalT > MinKineticEnergy))
  {
    // max. possible delta energy 
    Tmax = 2.*electron_mass_c2*E*(E+2.*mass)/
           (mass*mass+2.*electron_mass_c2*(E+mass)+
            electron_mass_c2*electron_mass_c2) ;

    rcut=G4Electron::Electron()->GetCuts();

    if(Tc > Tmax) Tc=Tmax ;
    // generate subcutoff delta rays only if Tc>MinDeltaEnergyNow!
    if((Tc > MinDeltaEnergyNow) && (Tmax > MinDeltaEnergyNow))
    {
      presafety  = stepData.GetPreStepPoint()->GetSafety() ;
     // postsafety = stepData.GetPostStepPoint()->GetSafety() ;

      G4Navigator *navigator=
         G4TransportationManager::GetTransportationManager()
                                   ->GetNavigatorForTracking();
      postsafety =
          navigator->ComputeSafety(stepData.GetPostStepPoint()->GetPosition());

      safety = G4std::min(presafety,postsafety) ;

      if(safety < rcut)
     {

        x1=stepData.GetPreStepPoint()->GetPosition().x();
        y1=stepData.GetPreStepPoint()->GetPosition().y();
        z1=stepData.GetPreStepPoint()->GetPosition().z();
        dx=stepData.GetPostStepPoint()->GetPosition().x()-x1 ;
        dy=stepData.GetPostStepPoint()->GetPosition().y()-y1 ;
        dz=stepData.GetPostStepPoint()->GetPosition().z()-z1 ;
        time0=stepData.GetPreStepPoint()->GetGlobalTime();
        dTime=stepData.GetPostStepPoint()->GetGlobalTime()-time0;

        if((presafety<rcut)&&(postsafety<rcut))
        {
          fragment = Step ;
          frperstep=1. ;
        }
        else if(presafety<rcut)
        {
          delta=presafety*Step/(postsafety-presafety) ;
          fragment=rcut*(Step+delta)/postsafety-delta ;
          frperstep=fragment/Step;
        }
        else if(postsafety<rcut)
        {
          delta=postsafety*Step/(presafety-postsafety) ;
          fragment=rcut*(Step+delta)/presafety-delta ;
          x1 += dx;
          y1 += dy;
          z1 += dz;  
          time0 += dTime ;

          frperstep=-fragment/Step;
        }

      if(fragment>0.)
      {
        T0=G4EnergyLossTables::GetPreciseEnergyFromRange(
                                             G4Electron::Electron(),
                                             G4std::min(presafety,postsafety),
                                             aMaterial) ;

        // absolute lower limit for T0
        if(T0<MinDeltaEnergyNow) T0=MinDeltaEnergyNow ;

        // compute nb of delta rays to be generated
        G4int N=int(fragment*(c0N/(E*T0)+c1N/T0-(c2N+c3N*T0)/Tc)* 
                (aMaterial->GetTotNbOfElectPerVolume())+0.5) ;

        G4double Px,Py,Pz ;
        G4ThreeVector ParticleDirection ;
        ParticleDirection=stepData.GetPostStepPoint()->
                                   GetMomentumDirection() ;
        Px =ParticleDirection.x() ;
        Py =ParticleDirection.y() ;
        Pz =ParticleDirection.z() ;
     
        G4int subdelta = 0;

        if(N > 0)
        {
          G4double Tkin,Etot,P,T,p,costheta,sintheta,phi,dirx,diry,dirz,
                   Pnew,delToverTc,
                   delTkin,delLoss,rate,
                   urandom ;
          G4StepPoint *point ;
  
          Tkin = E ;
          Etot = Tkin+mass ;
          P    = sqrt(Tkin*(Etot+mass)) ;

          aParticleChange.SetNumberOfSecondaries(N);
          do {
               subdelta += 1 ;

               Tmax = 2.*electron_mass_c2*Tkin*(Tkin+2.*mass)/
                      (mass*mass+2.*electron_mass_c2*(Tkin+mass)+
                        electron_mass_c2*electron_mass_c2) ;

               if(Tc>Tmax) Tc = Tmax ;

               //check if there is enough energy ....
               if((Tkin>TmintoProduceDelta)&&(Tc > T0)&&(MeanLoss>0.))
               {
                 delToverTc=1.-T0/Tc ;
                 T=T0/(1.-delToverTc*G4UniformRand()) ;
                 if(T > MeanLoss) T=MeanLoss ;
                 MeanLoss -= T ;
                 p=sqrt(T*(T+2.*electron_mass_c2)) ;

                 costheta = T*(Etot+electron_mass_c2)/(P*p) ;
                 if(costheta<-1.) costheta=-1.;
                 if(costheta> 1.) costheta= 1.;

                 phi=twopi*G4UniformRand() ;
                 sintheta=sqrt(1.-costheta*costheta);
                 dirx=sintheta*cos(phi);
                 diry=sintheta*sin(phi);
                 dirz=costheta;
                
               urandom = G4UniformRand() ;
               // distribute x,y,z along Pre-Post !
               G4double xd,yd,zd ;
               xd=x1+frperstep*dx*urandom ;
               yd=y1+frperstep*dy*urandom ;
               zd=z1+frperstep*dz*urandom ;
               G4ThreeVector DeltaPosition(xd,yd,zd) ;
               DeltaTime=time0+frperstep*dTime*urandom ;

               G4ThreeVector DeltaDirection(dirx,diry,dirz) ;
               DeltaDirection.rotateUz(ParticleDirection);

               G4DynamicParticle* theDelta = new G4DynamicParticle ;
               theDelta->SetDefinition(G4Electron::Electron());
               theDelta->SetKineticEnergy(T);

               theDelta->SetMomentumDirection(DeltaDirection.x(),
                              DeltaDirection.y(),DeltaDirection.z());

               // update initial particle,fill ParticleChange
               Tkin -= T ;
               Px =(P*ParticleDirection.x()-p*DeltaDirection.x()) ;
               Py =(P*ParticleDirection.y()-p*DeltaDirection.y()) ;
               Pz =(P*ParticleDirection.z()-p*DeltaDirection.z()) ;
               Pnew = sqrt(Px*Px+Py*Py+Pz*Pz) ;
               Px /= Pnew ;
               Py /= Pnew ;
               Pz /= Pnew ;
               P  = Pnew ;
               G4ThreeVector ParticleDirectionnew(Px,Py,Pz) ;
               ParticleDirection = ParticleDirectionnew;

               G4Track* deltaTrack =
                        new G4Track(theDelta,DeltaTime,DeltaPosition);
               deltaTrack->
                SetTouchable(stepData.GetPostStepPoint()->GetTouchable()) ;

               deltaTrack->SetParentID(trackData.GetTrackID()) ;

               aParticleChange.AddSecondary(deltaTrack) ;

               }

             } while (subdelta<N) ;

             // update the particle direction and kinetic energy
             if(subdelta > 0)
               aParticleChange.SetMomentumChange(Px,Py,Pz) ;
             E = Tkin ;
           }
          }
         }
       }
     }
  }
  //   end of subcutoff business 

  finalT = E - MeanLoss ;
  if(finalT < MinKineticEnergy) finalT = 0. ;

  //  now the loss with fluctuation
  if((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowerBoundEloss))
  {
    finalT = E -
        EnergyLossFluctuation(aParticle,aMaterial,ChargeSquare,MeanLoss,Step) ;
    if (finalT < 0.) finalT = 0.  ;
  }

  //  kill the particle if the kinetic energy <= 0  
  if (finalT <= 0. )
  {
    finalT = 0.;
    if(aParticle->GetDefinition()->GetParticleName() == "proton")
      aParticleChange.SetStatusChange(fStopAndKill);
    else  
      aParticleChange.SetStatusChange(fStopButAlive); 
  } 

  // aParticleChange.SetNumberOfSecondaries(0);
  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VhEnergyLoss::EnergyLossFluctuation(
                                   const G4DynamicParticle *aParticle,
                                         G4Material *aMaterial,
                                         G4double ChargeSquare,
                                         G4double MeanLoss,
                                         G4double Step)
{
  G4double loss = GetLossWithFluct(aParticle,aMaterial,
                                   ChargeSquare,MeanLoss,Step) ;
  return loss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

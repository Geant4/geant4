// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VeEnergyLossPlus.cc,v 1.1 2000-04-25 14:33:10 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//  
// -----------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4VeEnergyLossPlus physics process -----------
//                by Laszlo Urban, 20 March 1997 
// **************************************************************
// It is the first implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of e+/e-.
// --------------------------------------------------------------
// 18/11/98  , L. Urban
//  It is a modified version of G4VeEnergyLoss:
//  continuous energy loss with generation of subcutoff delta rays
// 02/02/99  important correction in AlongStepDoIt , L.Urban
// 28/04/99  bug fixed (unit independece now),L.Urban
// 10/02/00  modifications , new e.m. structure, L.Urban
// --------------------------------------------------------------
 
#include "G4VeEnergyLossPlus.hh"
#include "G4EnergyLossMessenger.hh"
#include "G4Poisson.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Initialisation of static data members
// -------------------------------------

G4int            G4VeEnergyLossPlus::NbOfProcesses  = 2;
G4int            G4VeEnergyLossPlus::CounterOfElectronProcess = 0;
G4int            G4VeEnergyLossPlus::CounterOfPositronProcess = 0;
G4PhysicsTable** G4VeEnergyLossPlus::RecorderOfElectronProcess =
                                           new G4PhysicsTable*[10];
G4PhysicsTable** G4VeEnergyLossPlus::RecorderOfPositronProcess =
                                           new G4PhysicsTable*[10];
                                           
G4double         G4VeEnergyLossPlus::MinDeltaCutInRange = 0.010*mm ;
G4double*        G4VeEnergyLossPlus::MinDeltaEnergy     = NULL   ;

G4PhysicsTable*  G4VeEnergyLossPlus::theDEDXElectronTable         = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theDEDXPositronTable         = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theRangeElectronTable        = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theRangePositronTable        = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theInverseRangeElectronTable = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theInverseRangePositronTable = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theLabTimeElectronTable      = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theLabTimePositronTable      = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theProperTimeElectronTable   = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theProperTimePositronTable   = NULL;

G4PhysicsTable*  G4VeEnergyLossPlus::theeRangeCoeffATable         = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theeRangeCoeffBTable         = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::theeRangeCoeffCTable         = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::thepRangeCoeffATable         = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::thepRangeCoeffBTable         = NULL;
G4PhysicsTable*  G4VeEnergyLossPlus::thepRangeCoeffCTable         = NULL;

G4double G4VeEnergyLossPlus::LowerBoundEloss = 1.*keV ;
G4double G4VeEnergyLossPlus::UpperBoundEloss = 100.*TeV ;
G4int    G4VeEnergyLossPlus::NbinEloss = 100 ;
G4double G4VeEnergyLossPlus::RTable,G4VeEnergyLossPlus::LOGRTable;

G4EnergyLossMessenger* G4VeEnergyLossPlus::eLossMessenger         = NULL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// constructor and destructor
 
G4VeEnergyLossPlus::G4VeEnergyLossPlus(const G4String& processName)
   : G4VEnergyLoss (processName),
     theLossTable(NULL),
     theDEDXTable(NULL),
     Charge(-1.),lastCharge(0.),
     MinKineticEnergy(1.*eV),
     linLossLimit(0.05),
     c1N(2.86e-23*MeV*mm*mm),
     c2N(c1N*MeV/10.),
     Ndeltamax(100)
{
 //create (only once) EnergyLoss messenger 
 if(!eLossMessenger) eLossMessenger = new G4EnergyLossMessenger();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VeEnergyLossPlus::~G4VeEnergyLossPlus() 
{
     if (theLossTable) 
       {
         theLossTable->clearAndDestroy();
         delete theLossTable;
         if(MinDeltaEnergy) delete MinDeltaEnergy ;
       }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

void G4VeEnergyLossPlus::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{
  ParticleMass = aParticleType.GetPDGMass(); 

  //  calculate data members LOGRTable,RTable first
  G4double lrate = log(UpperBoundEloss/LowerBoundEloss);
  LOGRTable=lrate/NbinEloss;
  RTable   =exp(LOGRTable);

  // Build energy loss table as a sum of the energy loss due to the
  // different processes.                                           
  //

  const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
  
  // create table for the total energy loss

  if (&aParticleType==G4Electron::Electron())
    {
      RecorderOfProcess=RecorderOfElectronProcess;
      CounterOfProcess=CounterOfElectronProcess;
      if (CounterOfProcess == NbOfProcesses)
        {
         if (theDEDXElectronTable)
           { 
             theDEDXElectronTable->clearAndDestroy();
             delete theDEDXElectronTable; 
           }
         theDEDXElectronTable = new G4PhysicsTable(numOfMaterials);
         theDEDXTable = theDEDXElectronTable;
        }
    }
  if (&aParticleType==G4Positron::Positron())
    {
     RecorderOfProcess=RecorderOfPositronProcess;
     CounterOfProcess=CounterOfPositronProcess;
     if (CounterOfProcess == NbOfProcesses)
       {
        if (theDEDXPositronTable)
          { 
            theDEDXPositronTable->clearAndDestroy();
            delete theDEDXPositronTable; 
          }
        theDEDXPositronTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXPositronTable;
       }
    }

  if (CounterOfProcess == NbOfProcesses)
    {
     // fill the tables
     // loop for materials
     G4double LowEdgeEnergy , Value;
     G4bool isOutRange;
     G4PhysicsTable* pointer;

     for (G4int J=0; J<numOfMaterials; J++)
        {
         // create physics vector and fill it
 
         G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowerBoundEloss, UpperBoundEloss, NbinEloss);   

         // loop for the kinetic energy
   
         for (G4int i=0; i<NbinEloss; i++) 
            {
              LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;      
              //here comes the sum of the different tables created by the  
              //processes (ionisation,bremsstrahlung,etc...)              
              Value = 0.;    
              for (G4int process=0; process < NbOfProcesses; process++)
                 {
                   pointer= RecorderOfProcess[process];
                   Value += (*pointer)[J]->GetValue(LowEdgeEnergy,isOutRange);
                 }

              aVector->PutValue(i,Value) ; 
            }

         theDEDXTable->insert(aVector) ;

        }

 
     //reset counter to zero
     if (&aParticleType==G4Electron::Electron()) CounterOfElectronProcess=0;
     if (&aParticleType==G4Positron::Positron()) CounterOfPositronProcess=0;

     ParticleMass = aParticleType.GetPDGMass();

     if (&aParticleType==G4Electron::Electron())
     {
       // Build range table
       theRangeElectronTable = BuildRangeTable(theDEDXElectronTable,
                                               theRangeElectronTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build lab/proper time tables
       theLabTimeElectronTable = BuildLabTimeTable(theDEDXElectronTable,
                         theLabTimeElectronTable,
                         LowerBoundEloss,UpperBoundEloss,NbinEloss);
       theProperTimeElectronTable = BuildProperTimeTable(theDEDXElectronTable,
                            theProperTimeElectronTable,
                            LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build coeff tables for the energy loss calculation
       theeRangeCoeffATable = BuildRangeCoeffATable(theRangeElectronTable,
                             theeRangeCoeffATable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       theeRangeCoeffBTable = BuildRangeCoeffBTable(theRangeElectronTable,
                             theeRangeCoeffBTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       theeRangeCoeffCTable = BuildRangeCoeffCTable(theRangeElectronTable,
                             theeRangeCoeffCTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // invert the range table
       theInverseRangeElectronTable = BuildInverseRangeTable(theRangeElectronTable,
                              theeRangeCoeffATable,
                              theeRangeCoeffBTable,
                              theeRangeCoeffCTable,
                              theInverseRangeElectronTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);
     }
     if (&aParticleType==G4Positron::Positron())
     {
        // Build range table
       theRangePositronTable = BuildRangeTable(theDEDXPositronTable,
                                               theRangePositronTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);


       // Build lab/proper time tables
       theLabTimePositronTable = BuildLabTimeTable(theDEDXPositronTable,
                         theLabTimePositronTable,
                         LowerBoundEloss,UpperBoundEloss,NbinEloss);
       theProperTimePositronTable = BuildProperTimeTable(theDEDXPositronTable,
                            theProperTimePositronTable,
                            LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build coeff tables for the energy loss calculation
       thepRangeCoeffATable = BuildRangeCoeffATable(theRangePositronTable,
                             thepRangeCoeffATable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepRangeCoeffBTable = BuildRangeCoeffBTable(theRangePositronTable,
                             thepRangeCoeffBTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepRangeCoeffCTable = BuildRangeCoeffCTable(theRangePositronTable,
                             thepRangeCoeffCTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // invert the range table
       theInverseRangePositronTable = BuildInverseRangeTable(theRangePositronTable,
                              thepRangeCoeffATable,
                              thepRangeCoeffBTable,
                              thepRangeCoeffCTable,
                              theInverseRangePositronTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

     }

     // make the energy loss and the range table available
     G4EnergyLossTables::Register(&aParticleType,  
       (&aParticleType==G4Electron::Electron())?
       theDEDXElectronTable: theDEDXPositronTable,
       (&aParticleType==G4Electron::Electron())?
       theRangeElectronTable: theRangePositronTable,
       (&aParticleType==G4Electron::Electron())?
       theInverseRangeElectronTable: theInverseRangePositronTable,
       (&aParticleType==G4Electron::Electron())?
       theLabTimeElectronTable: theLabTimePositronTable,
       (&aParticleType==G4Electron::Electron())?
       theProperTimeElectronTable: theProperTimePositronTable,
       LowerBoundEloss, UpperBoundEloss, 1.,NbinEloss);

     if(&aParticleType==G4Electron::Electron())
     {
       // create array for the min. delta cuts in kinetic energy 
       G4double absLowerLimit = 1.*keV ;
       G4cout << G4endl;
       G4cout.precision(5) ;
       G4cout << " eIoni+ Minimum Delta cut in range=" << MinDeltaCutInRange/mm
              << "  mm." << G4endl;
       G4cout << " min. delta energies (keV) " << G4endl;
       G4cout << "   material         min.delta energy " << G4endl;
       G4cout << G4endl;

       if(MinDeltaEnergy) delete MinDeltaEnergy ;
       MinDeltaEnergy = new G4double [numOfMaterials] ; 
       G4double Tlowerlimit = 1.*keV ;
       for(G4int mat=0; mat<numOfMaterials; mat++)
       {
         MinDeltaEnergy[mat] = G4EnergyLossTables::GetPreciseEnergyFromRange(
                               G4Electron::Electron(),MinDeltaCutInRange,
                                          (*theMaterialTable)(mat)) ;

         if(MinDeltaEnergy[mat]<absLowerLimit)
           MinDeltaEnergy[mat] = absLowerLimit ; 

         if(MinDeltaEnergy[mat]<Tlowerlimit) MinDeltaEnergy[mat]=Tlowerlimit ;
         G4cout << G4std::setw(20) << (*theMaterialTable)(mat)->GetName() 
                << G4std::setw(15) << MinDeltaEnergy[mat]/keV << G4endl;
       }
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
      
G4VParticleChange* G4VeEnergyLossPlus::AlongStepDoIt( const G4Track& trackData,
                                                 const G4Step&  stepData)
{                              
 // compute the energy loss after a Step

  static const G4double faclow = 1.5 ;

  // get particle and material pointers from trackData
  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle();
  G4double E      = aParticle->GetKineticEnergy() ;
 
  G4Material* aMaterial = trackData.GetMaterial();
  G4int index = aMaterial->GetIndex();
  G4double Step = stepData.GetStepLength();
 
  aParticleChange.Initialize(trackData);

  G4double MeanLoss, finalT;
 
  if (E < MinKineticEnergy)   finalT = 0.;
 
  else if (E<faclow*LowerBoundEloss)
  {
    if (Step >= fRangeNow)  finalT = 0.;
   // else finalT = E*(1.-Step/fRangeNow) ;
    else finalT = E*(1.-sqrt(Step/fRangeNow)) ;
  }
   
  else if (E>=UpperBoundEloss) finalT = E - Step*fdEdx;

  else if (Step >= fRangeNow)  finalT = 0.;
 
  else
  {
    if(Step/fRangeNow < linLossLimit) finalT = E-Step*fdEdx ;
    else
    {
      if (Charge<0.) finalT = G4EnergyLossTables::GetPreciseEnergyFromRange
                             (G4Electron::Electron(),fRangeNow-Step,aMaterial);
      else           finalT = G4EnergyLossTables::GetPreciseEnergyFromRange
                             (G4Positron::Positron(),fRangeNow-Step,aMaterial);
     }
  }

  if(finalT < MinKineticEnergy) finalT = 0. ;

  MeanLoss = E - finalT ;
 
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  G4double MinDeltaEnergyNow = MinDeltaEnergy[index] ;
  G4double TmintoProduceDelta=0.5*(3.-Charge)*MinDeltaEnergyNow ;
  if((E > TmintoProduceDelta) && (MeanLoss > MinDeltaEnergyNow)
                                   && (finalT > MinKineticEnergy)) 
  {
    G4double rcut,Tc,T0,presafety,postsafety,safety,
             delta,fragment ;
    G4double frperstep,x1,y1,z1,dx,dy,dz,dTime,time0,DeltaTime;
    G4double epsil= MinKineticEnergy/2. ;

    if(Charge < 0.)
    {
      rcut=G4Electron::Electron()->GetCuts();
      Tc=G4Electron::Electron()->GetCutsInEnergy()[index];
      // threshold !
      if(Tc > 0.5*E) Tc=0.5*E ;
    }
    else
    {
      rcut=G4Positron::Positron()->GetCuts();
      Tc=G4Positron::Positron()->GetCutsInEnergy()[index];
      // threshold !
      if(Tc > E) Tc=E ;
    }
    // generate subcutoff delta rays only if Tc>MinDeltaEnergy!
    if(Tc > MinDeltaEnergyNow)
    {
      presafety  = stepData.GetPreStepPoint()->GetSafety() ;
       
      // postsafety = stepData.GetPostStepPoint()->GetSafety() ;
     
      G4Navigator *navigator=
         G4TransportationManager::GetTransportationManager()
                                   ->GetNavigatorForTracking();
      postsafety =
          navigator->ComputeSafety(stepData.GetPostStepPoint()->GetPosition());
        
      safety=G4std::min(presafety,postsafety);

      if(safety<rcut)
      {
        T0=G4EnergyLossTables::GetPreciseEnergyFromRange(
                     G4Electron::Electron(),safety,aMaterial) ;   

        // absolute lower limit for T0 
        if(T0<MinDeltaEnergyNow) T0=MinDeltaEnergyNow ;
 // ..................................................................

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
        // compute nb of delta rays to be generated
        G4int N=int(fragment*(c1N*(1.-T0/Tc)+c2N/E)*
                (aMaterial->GetTotNbOfElectPerVolume())/T0+0.5) ;
        if(N > Ndeltamax)
           N = Ndeltamax ;
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
                   sumT,delTkin,delLoss,rate,
                   urandom ;
          G4StepPoint *point ;
   
          sumT=0.;

          Tkin = E ;
          Etot = Tkin+electron_mass_c2 ;
          P    = sqrt(Tkin*(Etot+electron_mass_c2)) ;

          aParticleChange.SetNumberOfSecondaries(N);
          do {
               subdelta += 1 ;

               if((Charge<0.)&&(Tc>0.5*Tkin)) Tc=0.5*Tkin ;
               if((Charge>0.)&&(Tc>    Tkin)) Tc=    Tkin ;
  
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

               sumT += T ;

               urandom = G4UniformRand() ;
               // distribute x,y,z along Pre-Post !
               G4double xd,yd,zd ;
               xd=x1+frperstep*dx*urandom ;
               yd=y1+frperstep*dy*urandom ;
               zd=z1+frperstep*dz*urandom ;
               G4ThreeVector DeltaPosition(xd,yd,zd) ;
               DeltaTime=time0+frperstep*dTime*urandom ;
               ParticleDirection=stepData.GetPostStepPoint()->
                                     GetMomentumDirection() ;    
                    
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
 // ................................................................
     }
    }
  }

  finalT = E - MeanLoss ;
  if(finalT < MinKineticEnergy) finalT = 0. ;

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
  //now the loss with fluctuation
  if ((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowerBoundEloss))
  {
    finalT = E-GetLossWithFluct(aParticle,aMaterial,MeanLoss);
    if (finalT < 0.) finalT = 0. ;
  }

  // kill the particle if the kinetic energy <= 0  
   if (finalT <= 0. )
   {
     finalT = 0.;
     if (Charge < 0.) aParticleChange.SetStatusChange(fStopAndKill);
     else             aParticleChange.SetStatusChange(fStopButAlive); 
   } 

  aParticleChange.SetEnergyChange(finalT);
  aParticleChange.SetLocalEnergyDeposit(E-finalT);

  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


   

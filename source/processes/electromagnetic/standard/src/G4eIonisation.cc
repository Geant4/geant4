// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4eIonisation.cc,v 1.4 1999-12-15 14:51:52 gunter Exp $
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
//      ---------- G4eIonisation physics process -----------
//                by Laszlo Urban, 20 March 1997 
// **************************************************************
// It is the first implementation of the NEW IONISATION PROCESS.
// It calculates the ionisation of e+/e-.
// **************************************************************
//
// 07-04-98: remove 'tracking cut' of the ionizing particle, MMa 
// 04-09-98: new methods SetBining() PrintInfo()
// 07-09-98: Cleanup
// 02/02/99: correction inDoIt , L.Urban
// --------------------------------------------------------------
 

#include "G4eIonisation.hh"
#include "G4EnergyLossTables.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

// constructor and destructor
 
G4eIonisation::G4eIonisation(const G4String& processName)
   : G4eEnergyLoss(processName),
     theMeanFreePathTable(NULL),
     LowestKineticEnergy(1.*keV),
     HighestKineticEnergy(100.*TeV),
     TotBin(100)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisation::~G4eIonisation() 
{
     if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
     }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisation::SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins)
{
  LowestKineticEnergy = lowE;  HighestKineticEnergy = highE; TotBin = nBins;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
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
 
    BuildLambdaTable(aParticleType) ;
 
    BuildDEDXTable(aParticleType) ;
                                             
  if(&aParticleType==G4Electron::Electron())
    PrintInfoDefinition();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
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
              d = G4std::min(ParticleCutInKineticEnergy[J], Tmax)/ParticleMass;
              ionloss = log(2.*(tau+2.)/Eexcm2)-1.-beta2
                       + log((tau-d)*d)+tau/(tau-d)
                       + (0.5*d*d+(2.*tau+1.)*log(1.-d/tau))/gamma2;
            }
          else        //positron
            {
              Tmax = LowEdgeEnergy ;  
              d = G4std::min(ParticleCutInKineticEnergy[J], Tmax)/ParticleMass;
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

void G4eIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{
  // Build mean free path tables for the delta ray production process
  //     tables are built for MATERIALS 

  G4double LowEdgeEnergy, Value, SIGMA;

  //create table
  const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

  if (theMeanFreePathTable) { theMeanFreePathTable->clearAndDestroy();
                              delete theMeanFreePathTable;
                            }

  theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);

  // get electron  cuts in kinetic energy
  G4double* DeltaCutInKineticEnergy = G4Electron::Electron()->GetCutsInEnergy() ;
 

  // loop for materials 

 for (G4int J=0 ; J < numOfMaterials; J++)
    { 
     //create physics vector then fill it ....

     G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
               LowestKineticEnergy, HighestKineticEnergy, TotBin);

     // compute the (macroscopic) cross section first
 
     const G4Material* material= (*theMaterialTable)[J];        
     const 
     G4ElementVector* theElementVector = material->GetElementVector();
     const
     G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
     const
     G4int NumberOfElements = material->GetNumberOfElements() ;
 
     // get the electron kinetic energy cut for the actual material,
     // it will be used in ComputeMicroscopicCrossSection
     // (--> it will be the same for all the elements in this material )
     G4double DeltaThreshold = DeltaCutInKineticEnergy[J] ;

     for (G4int i = 0 ; i < TotBin ; i++)
        {
          LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
          SIGMA = 0.;           
          for (G4int iel=0; iel<NumberOfElements; iel++ )
             {
                SIGMA += theAtomicNumDensityVector[iel]*
                         ComputeMicroscopicCrossSection( aParticleType,
                                                         LowEdgeEnergy,
                                      (*theElementVector)(iel)->GetZ(),
                                                       DeltaThreshold);
             }

          // mean free path = 1./macroscopic cross section
          Value = SIGMA > DBL_MIN ? 1./SIGMA : DBL_MAX;     
          aVector->PutValue(i, Value) ;
        }
     theMeanFreePathTable->insert(aVector);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eIonisation::ComputeMicroscopicCrossSection(
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber ,
                                 G4double DeltaThreshold)
{
  // calculates the microscopic cross section
  //(it is called for elements , AtomicNumber = Z )
 
  G4double MaxKineticEnergyTransfer, TotalCrossSection(0.);
  
  ParticleMass = aParticleType.GetPDGMass();
  G4double TotalEnergy = KineticEnergy + ParticleMass;

  G4double betasquare = KineticEnergy*(TotalEnergy+ParticleMass)
                       /(TotalEnergy*TotalEnergy);
  G4double gamma = TotalEnergy/ParticleMass, gamma2 = gamma*gamma;
  G4double x=DeltaThreshold/KineticEnergy, x2 = x*x;

  if (&aParticleType==G4Electron::Electron())
                            MaxKineticEnergyTransfer = 0.5*KineticEnergy;
  else                      MaxKineticEnergyTransfer =     KineticEnergy;

  // now you can calculate the total cross section

 if (MaxKineticEnergyTransfer > DeltaThreshold)
   {
    if (&aParticleType==G4Electron::Electron())   //Moller (e-e-) scattering
      {
        TotalCrossSection  = (gamma-1.)*(gamma-1.)*(0.5-x)/gamma2 + 1./x
                            - 1./(1.-x)-(2.*gamma-1.)*log((1.-x)/x)/gamma2; 
        TotalCrossSection /= betasquare; 
      }
    else                                         //Bhabha (e+e-) scattering
      {
       G4double y=1./(1.+gamma), y2 =y*y, y12=1.-2.*y;
       G4double b1=2.-y2, b2=y12*(3.+y2), b4=y12*y12*y12, b3=b4+y12*y12;
       TotalCrossSection = (1./x-1.)/betasquare+b1*log(x)+b2*(1.-x)
                          - b3*(1.-x2)/2.+b4*(1.-x2*x)/3.;
      }   
    TotalCrossSection *= (twopi_mc2_rcl2*AtomicNumber/KineticEnergy);
   }
 return TotalCrossSection ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
 
G4VParticleChange* G4eIonisation::PostStepDoIt( const G4Track& trackData,   
                                                const G4Step&  stepData)         
{
  aParticleChange.Initialize(trackData) ;
  
  G4Material*               aMaterial = trackData.GetMaterial() ;
  const G4DynamicParticle*  aParticle = trackData.GetDynamicParticle() ;

  ParticleMass = aParticle->GetDefinition()->GetPDGMass();
  G4double KineticEnergy = aParticle->GetKineticEnergy();
  G4double TotalEnergy = KineticEnergy + ParticleMass;
  G4double Psquare = KineticEnergy*(TotalEnergy+ParticleMass);
  G4double TotalMomentum = sqrt(Psquare);
  G4double Esquare=TotalEnergy*TotalEnergy;
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection();

  //  get kinetic energy cut for the electron
  G4double* DeltaCutInKineticEnergy = G4Electron::Electron()->GetCutsInEnergy() ;
  G4double DeltaThreshold = DeltaCutInKineticEnergy[aMaterial->GetIndex()];

  // some kinematics
  G4double MaxKineticEnergyTransfer;
  if (Charge < 0.) MaxKineticEnergyTransfer = 0.5*KineticEnergy;
  else             MaxKineticEnergyTransfer =     KineticEnergy;

  // sampling kinetic energy of the delta ray 

  if (MaxKineticEnergyTransfer <= DeltaThreshold) // pathological case (should not happen,
                                                  // there is no change at all)
     return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);


  // normal case
  G4double cc,y,y2,c2,b0,b1,b2,b3,b4,x,x1,grej,grejc;
 
  G4double tau = KineticEnergy/ParticleMass;
  G4double gamma = tau+1., gamma2=gamma*gamma;
  G4double xc = DeltaThreshold/KineticEnergy, xc1=1.-xc;

  if (Charge < 0.)  // Moller (e-e-) scattering
    { 
      b1=4./(9.*gamma2-10.*gamma+5.);
      b2=tau*tau*b1; b3=(2.*gamma2+2.*gamma-1.)*b1;
      cc=1.-2.*xc;       
      do { 
           x    = xc/(1.-cc*G4UniformRand()); x1 = 1.-x;
           grej = b2*x*x-b3*x/x1+b1*gamma2/(x1*x1);
         } while (G4UniformRand()>grej) ;
    }
  else             // Bhabha (e+e-) scattering
    {
      y=1./(gamma+1.); y2=y*y; cc=1.-2.*y;
      b1=2.-y2; b2=cc*(3.+y2);
      c2=cc*cc; b4=c2*cc; b3=c2+b4;
      b0=gamma2/(gamma2-1.);
      grejc=(((b4*xc-b3)*xc+b2)*xc-b1)*xc+b0;
      do {
           x    = xc/(1.-xc1*G4UniformRand());
           grej = ((((b4*x-b3)*x+b2)*x-b1)*x+b0)/grejc;
         } while (G4UniformRand()>grej);
    }
 
    G4double DeltaKineticEnergy = x * KineticEnergy;

  // protection :do not produce a secondary with 0. kinetic energy !
  if (DeltaKineticEnergy <= 0.)
      return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  G4double DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
                                                       2. * electron_mass_c2 ));
   
  G4double costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
                      /(DeltaTotalMomentum * TotalMomentum);

  if (costheta < -1.) costheta = -1.;
  if (costheta > +1.) costheta = +1.;

  //  direction of the delta electron

  G4double phi = twopi * G4UniformRand(); 
  G4double sintheta = sqrt((1.+costheta)*(1.-costheta));
  G4double dirx = sintheta * cos(phi), diry = sintheta * sin(phi), dirz = costheta;

  G4ThreeVector DeltaDirection(dirx,diry,dirz);
  DeltaDirection.rotateUz(ParticleDirection);

  // create G4DynamicParticle object for delta ray

  G4DynamicParticle* theDeltaRay = new G4DynamicParticle;
  theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
  theDeltaRay->SetMomentumDirection(
                   DeltaDirection.x(),DeltaDirection.y(),DeltaDirection.z()); 
  theDeltaRay->SetDefinition(G4Electron::Electron());
   
  // fill aParticleChange 
  // changed energy and momentum of the actual particle
  G4double finalKineticEnergy = KineticEnergy - DeltaKineticEnergy;
  
  G4double Edep = 0. ;

  if (finalKineticEnergy > MinKineticEnergy)
    {
      G4double finalPx = TotalMomentum*ParticleDirection.x()
                        - DeltaTotalMomentum*DeltaDirection.x(); 
      G4double finalPy = TotalMomentum*ParticleDirection.y()
                        - DeltaTotalMomentum*DeltaDirection.y(); 
      G4double finalPz = TotalMomentum*ParticleDirection.z()
                        - DeltaTotalMomentum*DeltaDirection.z(); 
      G4double finalMomentum =
                sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz) ;
      finalPx /= finalMomentum ;
      finalPy /= finalMomentum ;
      finalPz /= finalMomentum ;

      aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
    }
  else
    {
      finalKineticEnergy = 0.;
      Edep = finalKineticEnergy ;
      if (Charge < 0.) aParticleChange.SetStatusChange(fStopAndKill);
      else             aParticleChange.SetStatusChange(fStopButAlive);
    }
      
  aParticleChange.SetEnergyChange( finalKineticEnergy );
  aParticleChange.SetNumberOfSecondaries(1);  
  aParticleChange.AddSecondary( theDeltaRay );
  aParticleChange.SetLocalEnergyDeposit (Edep);
      
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisation::PrintInfoDefinition()
{
  G4String comments = "delta cross sections from Moller+Bhabha. ";
           comments += "Good description from 1 KeV to 100 GeV.\n";
           comments += "        delta ray energy sampled from  differential Xsection.";
                     
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,"Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy") 
         << " in " << TotBin << " bins. \n";
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
// $Id: G4hIonisation52.cc,v 1.4 2004/12/01 19:37:16 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//---------------- G4hIonisation52 physics process -------------------------------
//                 by Laszlo Urban, 30 May 1997
//------------------------------------------------------------------------------
//
// corrected by L.Urban on 24/09/97
// several bugs corrected by L.Urban on 13/01/98
// 07-04-98 remove 'tracking cut' of the ionizing particle, mma
// 22-10-98 cleanup L.Urban
// 02-02-99 bugs fixed , L.Urban
// 29-07-99 correction in BuildLossTable for low energy, L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 10-08-00 V.Ivanchenko change BuildLambdaTable, in order to
//          simulate energy losses of ions; correction to
//          cross section for particles with spin 1 is inserted as well
// 28-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 14-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 29-08-01 PostStepDoIt: correction for spin 1/2 (instead of 1) (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 25-09-01 completion of RetrievePhysicsTable() (mma)
// 29-10-01 all static functions no more inlined
// 08-11-01 Charge renamed zparticle; added to the dedx
// 27-03-02 Bug fix in scaling of lambda table (V.Ivanchenko)
// 09-04-02 Update calculation of tables for GenericIons (V.Ivanchenko)
// 10-06-02 bug fixed for stopping hadrons (V.Ivanchenko)
// 15-01-03 Migrade to cut per region (V.Ivanchenko)
// 10-03-03 Use SubType for GenericIons (V.Ivanchenko)
// 07-04-03 Fix problem of several runs (V.Ivanchenko)
// 08-04-03 finalRange is region aware (V.Ivanchenko)
// 17-04-03 fix problem of hadron tests (V.Ivanchenko)
// 26-04-03 fix problems of retrieve tables (V.Ivanchenko)
// 08-08-03 This class is frozen at the release 5.2 (V.Ivanchenko)
// 08-11-04 Remove of Store/Retrieve tables (V.Ivantchenko)
//
//------------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4hIonisation52.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4ProductionCutsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hIonisation52::LowerBoundLambda = 1.*keV;
G4double G4hIonisation52::UpperBoundLambda = 100.*TeV;
G4int	 G4hIonisation52::NbinLambda = 100;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4hIonisation52::G4hIonisation52(const G4String& processName)
   : G4VhEnergyLoss(processName),
     theMeanFreePathTable(0),
     Tmincut(1*keV)
{
  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4hIonisation52::~G4hIonisation52()
{
  if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hIonisation52::SetLowerBoundLambda(G4double val)
     {LowerBoundLambda = val;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hIonisation52::SetUpperBoundLambda(G4double val)
     {UpperBoundLambda = val;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hIonisation52::SetNbinLambda(G4int n)
     {NbinLambda = n;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hIonisation52::GetLowerBoundLambda()
         {return LowerBoundLambda;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hIonisation52::GetUpperBoundLambda()
         {return UpperBoundLambda;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4hIonisation52::GetNbinLambda()
      {return NbinLambda;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hIonisation52::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
// just call BuildLossTable+BuildLambdaTable
{

  if(verboseLevel > 0) {
    G4cout << "G4hIonisation52::BuildPhysicsTable for "
           << aParticleType.GetParticleName()
           << " mass(MeV)= " << aParticleType.GetPDGMass()/MeV
           << " charge= " << aParticleType.GetPDGCharge()/eplus
           << " type= " << aParticleType.GetParticleType()
           << G4endl;

    if(verboseLevel > 1) {
      G4cout << " MFPtable= " << theMeanFreePathTable
             << " DEDXtable= " << theDEDXpTable
             << " iniMass= " << initialMass
             << G4endl;
    }
  }

  if(aParticleType.GetParticleType() == "nucleus" &&
     aParticleType.GetParticleName() != "GenericIon" &&
     aParticleType.GetParticleSubType() == "generic")
  {

     G4EnergyLossTables::Register(&aParticleType,
              theDEDXpTable,
              theRangepTable,
              theInverseRangepTable,
              theLabTimepTable,
              theProperTimepTable,
              LowestKineticEnergy, HighestKineticEnergy,
              proton_mass_c2/aParticleType.GetPDGMass(),
              TotBin);

     return;
  }

  // get bining from EnergyLoss
  LowestKineticEnergy  = GetLowerBoundEloss();
  HighestKineticEnergy = GetUpperBoundEloss();
  TotBin               = GetNbinEloss();
  const G4ParticleDefinition* theProton = G4Proton::Proton();
  G4bool makeTables = false;

  if (aParticleType.GetPDGCharge() > 0.)
   {
     if( CutsWhereModified() || !theDEDXpTable )
    {
      BuildLossTable(*theProton);
      RecorderOfpProcess[0] = (*this).theLossTable;
//      CounterOfpProcess++;
      makeTables = true;
    }
   }
  else
   {
     if( CutsWhereModified() || !theDEDXpbarTable )
    {
      BuildLossTable(*(G4AntiProton::AntiProton())) ;
      RecorderOfpProcess[0] = (*this).theLossTable;
//      CounterOfpbarProcess++;
      makeTables = true;
    }
   }

  BuildLambdaTable(aParticleType);

  if( makeTables ) BuildDEDXTable(aParticleType);

  if(2 < verboseLevel) {
    G4cout << "MeanFreePathTable is built for "
           << aParticleType.GetParticleName() << G4endl;
    G4cout << (*theMeanFreePathTable) << G4endl;
  }


  if (&aParticleType == theProton)  PrintInfoDefinition();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hIonisation52::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  if(0 < verboseLevel) {
    G4cout << "G4hIonisation52::BuildLossTable() for process "
           << GetProcessName() << " and particle "
           << aParticleType.GetParticleName() << G4endl;
  }

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if (theLossTable) {theLossTable->clearAndDestroy(); delete theLossTable;}
  theLossTable = new G4PhysicsTable(numOfCouples);

  secondaryEnergyCuts = theCoupleTable->GetEnergyCutsVector(1);

  //  loop for materials
  //
  for (size_t J=0; J<numOfCouples; J++)
   {
    // create physics vector and fill it
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                        LowestKineticEnergy, HighestKineticEnergy, TotBin);

    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(J);
    const G4Material* material= couple->GetMaterial();

    // get  electron cut in kinetic energy for the material
    G4double DeltaThreshold = SecondaryEnergyThreshold(J);

    // now comes the loop for the kinetic energy values
    //
    for (G4int i = 0 ; i < TotBin ; i++)
     {
      G4double dEdx = ComputeRestrictedMeandEdx(aParticleType,
	                                        aVector->GetLowEdgeEnergy(i),
	                                        material,
	                                        DeltaThreshold);
      aVector->PutValue(i,dEdx);
     }
    theLossTable->insert(aVector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hIonisation52::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{

  if(0 < verboseLevel) {
    G4cout << "G4hIonisation52::BuildLambdaTable() for process "
           << GetProcessName() << " and particle "
           << aParticleType.GetParticleName() << G4endl;
  }

  //create table
  //
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if (theMeanFreePathTable)
    {theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable;}

  theMeanFreePathTable = new G4PhysicsTable(numOfCouples);

  // get electron cut in kinetic energy
  secondaryEnergyCuts = theCoupleTable->GetEnergyCutsVector(1);

  // loop for materials

  for (size_t J=0 ; J < numOfCouples; J++)
    {
     //create physics vector then fill it ....
     G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
               LowerBoundLambda,UpperBoundLambda,NbinLambda);

     // compute the (macroscopic) cross section first
     const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(J);
     const G4Material* material= couple->GetMaterial();

      // get  electron cut in kinetic energy for the material
     G4double DeltaThreshold = SecondaryEnergyThreshold(J);

     const G4ElementVector* theElementVector = material->GetElementVector();
     const G4double* NbOfAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
     G4int NumberOfElements = material->GetNumberOfElements();

     if(1 < verboseLevel) {
       G4cout << "### For material " << material->GetName()
              << " Tcut(MeV)= " << DeltaThreshold/MeV
              << " Tmin(MeV)= " << LowerBoundLambda/MeV
              << " Tmax(MeV)= " << UpperBoundLambda/MeV
              << " nbins= "  << NbinLambda
              << G4endl;
     }


     for ( G4int i = 0 ; i < NbinLambda ; i++ )
        {
          G4double LowEdgeEnergy = aVector->GetLowEdgeEnergy(i);
          G4double sigma = 0.;
          for (G4int iel=0; iel<NumberOfElements; iel++ )
            {
             sigma +=  NbOfAtomsPerVolume[iel]*
                       ComputeCrossSectionPerAtom(aParticleType,
                                                  LowEdgeEnergy,
		               (*theElementVector)[iel]->GetZ(),
		                               DeltaThreshold);
            }

          // mean free path = 1./macroscopic cross section
          G4double Value = sigma > DBL_MIN ? 1./sigma : DBL_MAX;
          aVector->PutValue(i, Value) ;
        }
     theMeanFreePathTable->insert(aVector);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hIonisation52::ComputeRestrictedMeandEdx (
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
				 const G4Material* material,
				 G4double DeltaThreshold)
{
  // calculate the dE/dx due to the ionization process (Geant4 internal units)
  // Bethe-Bloch formula
  //
  G4double particleMass   = proton_mass_c2;

  G4double ElectronDensity = material->GetElectronDensity();
  G4double Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double Eexc2 = Eexc*Eexc;

  G4double tau = KineticEnergy/particleMass;
  G4double gamma = tau + 1., bg2 = tau*(tau+2.), beta2 = bg2/(gamma*gamma);
  G4double RateMass = electron_mass_c2/particleMass;
  G4double Tmax=2.*electron_mass_c2*bg2/(1.+2.*gamma*RateMass+RateMass*RateMass);

  G4double taul = material->GetIonisation()->GetTaul();

  G4double dEdx = 0.;
  //
  // high energy part , Bethe-Bloch formula
  //
  if (tau > taul)
    {
     G4double rcut = min(DeltaThreshold/Tmax, 1.);
     dEdx = log(2.*electron_mass_c2*bg2*Tmax/Eexc2)
            +log(rcut)-(1.+rcut)*beta2;

     //density correction
     G4double Cden   = material->GetIonisation()->GetCdensity();
     G4double Mden   = material->GetIonisation()->GetMdensity();
     G4double Aden   = material->GetIonisation()->GetAdensity();
     G4double X0den  = material->GetIonisation()->GetX0density();
     G4double X1den  = material->GetIonisation()->GetX1density();

     const G4double twoln10 = 2.*log(10.);
     G4double  x = log(bg2)/twoln10;
     G4double delta;
     if (x < X0den) delta = 0.;
     else          {delta = twoln10*x - Cden;
                    if (x < X1den) delta += Aden*pow((X1den-x),Mden);
                   }

     // shell correction
     G4double* ShellCorrectionVector = material->GetIonisation()->
                                       GetShellCorrectionVector();
     const G4double bg2lim = 0.0169, taulim = 8.4146e-3;
     G4double sh = 0., xs = 1.;
     if (bg2 > bg2lim) for (G4int k=0; k<3; k++)
                          {xs *= bg2; sh += ShellCorrectionVector[k]/xs;}
     else { for (G4int k=0; k<3; k++)
                       {xs *= bg2lim; sh += ShellCorrectionVector[k]/xs;}
            sh *= log(tau/taul)/log(taulim/taul);
          }

     // now you can compute the total ionization loss
     dEdx -= (delta + sh);
     dEdx *= twopi_mc2_rcl2*ElectronDensity/beta2;
     if (dEdx < 0.) dEdx = 0.;
   }
  //
  //  low energy part , parametrized energy loss formulae
  //
  if (tau <= taul)
   {
     // get elements in the actual material,
     const G4ElementVector* theElementVector = material->GetElementVector();
     const G4double* NbOfAtomsPerVolume=material->GetVecNbOfAtomsPerVolume();
     G4int NumberOfElements = material->GetNumberOfElements();

     //  loop for the elements in the material
     dEdx = 0.;
     for (G4int iel=0; iel<NumberOfElements; iel++)
        {
          const G4Element* element = (*theElementVector)[iel];
          if (tau < element->GetIonisation()->GetTau0())
            dEdx += NbOfAtomsPerVolume[iel]
                       *(element->GetIonisation()->GetAlow()*sqrt(tau)
                       + element->GetIonisation()->GetBlow()*tau);
          else
            dEdx += NbOfAtomsPerVolume[iel]
                       * element->GetIonisation()->GetClow()/sqrt(tau);
        }
     G4double deltaloss = 0.;
     if (DeltaThreshold < Tmax)
       {
         deltaloss = log(Tmax/DeltaThreshold)-
                      beta2*(1.-DeltaThreshold/Tmax) ;
         if (aParticleType.GetPDGSpin() == 0.5)
            deltaloss += 0.25*(Tmax-DeltaThreshold)*(Tmax-DeltaThreshold)/
                 (KineticEnergy*KineticEnergy+proton_mass_c2*proton_mass_c2);
         deltaloss *= twopi_mc2_rcl2*ElectronDensity/beta2;
       }
     dEdx -= deltaloss;
     if (dEdx < 0.) dEdx = 0.;
   }
 return dEdx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hIonisation52::ComputeCrossSectionPerAtom(
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber,
				 G4double DeltaThreshold)
{
  // calculates the totalcross section per atom in GEANT4 internal units
  //    ( it is called for elements , AtomicNumber = Z )
  //
  // nb: cross section formula is OK for spin=0 and 1/2 only !

  initialMass = aParticleType.GetPDGMass();
  G4double particleMass = initialMass;

  G4double TotalEnergy = KineticEnergy + particleMass;

  G4double betasquare = KineticEnergy*(TotalEnergy+particleMass)
                      /(TotalEnergy*TotalEnergy);
  G4double tempvar = particleMass+electron_mass_c2;
  G4double MaxKineticEnergyTransfer = 2.*electron_mass_c2*KineticEnergy
                     *(TotalEnergy+particleMass)
                     /(tempvar*tempvar+2.*electron_mass_c2*KineticEnergy);

  G4double TotalCrossSection = 0.;
  if (MaxKineticEnergyTransfer > DeltaThreshold)
   {
     tempvar = DeltaThreshold/MaxKineticEnergyTransfer;
     TotalCrossSection = (1.-tempvar*(1.-betasquare*log(tempvar)))
                           /DeltaThreshold;

     G4double spin = aParticleType.GetPDGSpin();
     if (spin == 0.5)  TotalCrossSection +=  0.5
                       *(MaxKineticEnergyTransfer-DeltaThreshold)
                       /(TotalEnergy*TotalEnergy);
     if (spin == 1.)   TotalCrossSection +=
                       -log(tempvar)/(3.0*DeltaThreshold) +
	               (MaxKineticEnergyTransfer - DeltaThreshold) *
                       ((5.0+ 1.0/tempvar)*0.25 / (TotalEnergy*TotalEnergy) -
	               betasquare /
                       (MaxKineticEnergyTransfer * DeltaThreshold)) / 3.0;

     TotalCrossSection *= twopi_mc2_rcl2*AtomicNumber/betasquare;
   }
  return TotalCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4hIonisation52::PostStepDoIt(const G4Track& trackData,
                                               const G4Step&  stepData)
{
  aParticleChange.Initialize(trackData);

  const G4MaterialCutsCouple* couple = trackData.GetMaterialCutsCouple();
  const G4DynamicParticle*  aParticle = trackData.GetDynamicParticle();

  G4double particleMass = aParticle->GetMass();
  G4double KineticEnergy = aParticle->GetKineticEnergy();
  G4double TotalEnergy = KineticEnergy + particleMass;
  G4double Psquare = KineticEnergy*(TotalEnergy+particleMass);
  G4double Esquare = TotalEnergy*TotalEnergy;
  G4double betasquare=Psquare/Esquare;
  G4double summass = particleMass + electron_mass_c2;
  G4double MaxKineticEnergyTransfer = 2.*electron_mass_c2*Psquare
                      /(summass*summass+2.*electron_mass_c2*KineticEnergy);
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection();

  // get electron cut in kinetic energy

  G4double DeltaThreshold = SecondaryEnergyThreshold(couple->GetIndex());

  // sampling kinetic energy of the delta ray
  //
  if (MaxKineticEnergyTransfer <= DeltaThreshold)
    // pathological case (it should not happen, there is no change at all)
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  // normal case
  G4double xc = DeltaThreshold/MaxKineticEnergyTransfer;
  G4double rate = MaxKineticEnergyTransfer/TotalEnergy;
  G4double te2 = 0.;
  if (aParticle->GetDefinition()->GetPDGSpin() == 0.5) te2=0.5*rate*rate;

  // sampling follows ...
  G4double x,grej;
  G4double grejc=1.-betasquare*xc+te2*xc*xc;
  do { x=xc/(1.-(1.-xc)*G4UniformRand());
       grej=(1.-x*(betasquare-x*te2))/grejc;
     } while(G4UniformRand() > grej);

  G4double  DeltaKineticEnergy = x * MaxKineticEnergyTransfer;

  if (DeltaKineticEnergy <= 0.)
   return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  G4double DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
                                               2. * electron_mass_c2 ));
  G4double TotalMomentum = sqrt(Psquare);
  G4double costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
            /(DeltaTotalMomentum * TotalMomentum);

  if (costheta < -1.) costheta = -1.;
  if (costheta > +1.) costheta = +1.;

  //  direction of the delta electron
  //
  G4double phi = twopi*G4UniformRand();
  G4double sintheta = sqrt((1.+costheta)*(1.-costheta));
  G4double dirx = sintheta*cos(phi), diry = sintheta*sin(phi), dirz = costheta;

  G4ThreeVector DeltaDirection(dirx,diry,dirz);
  DeltaDirection.rotateUz(ParticleDirection);

  // create G4DynamicParticle object for delta ray
  //
  G4DynamicParticle *theDeltaRay = new G4DynamicParticle;
  theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
  theDeltaRay->SetMomentumDirection(
                   DeltaDirection.x(),DeltaDirection.y(),DeltaDirection.z());
  theDeltaRay->SetDefinition(G4Electron::Electron());

  // fill aParticleChange
  //
  G4double finalKineticEnergy = KineticEnergy - DeltaKineticEnergy;
  G4double Edep = 0;

  if (finalKineticEnergy > MinKineticEnergy)
   {
    G4double finalPx = TotalMomentum*ParticleDirection.x()
                      - DeltaTotalMomentum*DeltaDirection.x();
    G4double finalPy = TotalMomentum*ParticleDirection.y()
                      - DeltaTotalMomentum*DeltaDirection.y();
    G4double finalPz = TotalMomentum*ParticleDirection.z()
                      - DeltaTotalMomentum*DeltaDirection.z();
    G4double finalMomentum =
              sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
    finalPx /= finalMomentum;
    finalPy /= finalMomentum;
    finalPz /= finalMomentum;

    aParticleChange.ProposeMomentumDirection( finalPx,finalPy,finalPz );
   }
  else
   {
     Edep = finalKineticEnergy;
     finalKineticEnergy = 0.;
     if (!aParticle->GetDefinition()->GetProcessManager()->GetAtRestProcessVector()->size())
           aParticleChange.ProposeTrackStatus(fStopAndKill);
     else  aParticleChange.ProposeTrackStatus(fStopButAlive);
   }

  aParticleChange.ProposeEnergy( finalKineticEnergy );
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(theDeltaRay);
  aParticleChange.ProposeLocalEnergyDeposit (Edep);

  //ResetNumberOfInteractionLengthLeft();
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hIonisation52::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . "
            "\n         Good description above the mean excitation energy.\n"
            "         delta ray energy sampled from  differential Xsection.";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
                                                  "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. "
	 << "\n        Step function: finalRange(mm)= " << finalRange
	 << ",  dRoverRange= " << dRoverRange
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VEnergyLoss.cc,v 1.46 2006-06-29 19:55:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// --------------------------------------------------------------
//
//  bug fixed in fluct., L.Urban 01/02/01
//  bug fixed in fluct., L.Urban 26/05/00
//  bug fixed in fluct., L.Urban 22/11/00
//  bugfix in fluct.
//  (some variables are doubles instead of ints now),L.Urban 23/03/01
//  18/05/01 V.Ivanchenko Clean up againist Linux ANSI compilation
//  17-09-01 migration of Materials to pure STL (mma)
//  26-10-01 static inline functions moved from .hh file (mma)
//  08.11.01 some static methods,data members are not static L.Urban
//  11.02.02 subSecFlag = false --> No sucutoff generation (mma)
//  14.02.02 initial value of data member finalRange has been changed L.Urban
//  26.02.02 initial value of data member finalRange = 1 mm (mma)
//  21.07.02 V.Ivanchenko Fix at low energies - if tmax below ionisation
//           potential then only Gaussian fluctuations are sampled.
//  15.01.03 Migrade to cut per region (V.Ivanchenko)
//  05.02.03 Minor fix for several region case (V.Ivanchenko)
//  25.03.03 add finalRangeRequested (mma)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEnergyLoss.hh"
#include "G4EnergyLossMessenger.hh"
#include "G4ProductionCutsTable.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool       G4VEnergyLoss::rndmStepFlag   = false;
G4bool       G4VEnergyLoss::EnlossFlucFlag = true;

G4bool       G4VEnergyLoss::subSecFlag     = false;
G4double     G4VEnergyLoss::MinDeltaCutInRange = 0.100*mm;
G4double*    G4VEnergyLoss::MinDeltaEnergy	    = 0;
G4bool*      G4VEnergyLoss::LowerLimitForced        = 0;
G4bool	     G4VEnergyLoss::setMinDeltaCutInRange = false;

G4double     G4VEnergyLoss::dRoverRange    = 20*perCent;
G4double     G4VEnergyLoss::finalRange     = 1*mm;
G4double     G4VEnergyLoss::finalRangeRequested = -1*mm;
G4double     G4VEnergyLoss::c1lim = dRoverRange;
G4double     G4VEnergyLoss::c2lim = 2.*(1.-dRoverRange)*finalRange;
G4double     G4VEnergyLoss::c3lim = -(1.-dRoverRange)*finalRange*finalRange;

G4EnergyLossMessenger* G4VEnergyLoss::ELossMessenger  = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLoss::G4VEnergyLoss()
                   :G4VContinuousDiscreteProcess("No Name Loss Process"),
     lastMaterial(NULL),
     nmaxCont1(4),
     nmaxCont2(16)
{
  G4Exception("G4VEnergyLoss:: default constructor is called");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLoss::G4VEnergyLoss(const G4String& aName , G4ProcessType aType)
                  : G4VContinuousDiscreteProcess(aName, aType),
     lastMaterial(NULL),
     nmaxCont1(4),
     nmaxCont2(16)
{
  //create (only once) EnergyLoss messenger
//  if(!ELossMessenger) ELossMessenger = new G4EnergyLossMessenger();
  if(!ELossMessenger) ELossMessenger = G4LossTableManager::Instance()->GetMessenger();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLoss::~G4VEnergyLoss()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLoss::G4VEnergyLoss(G4VEnergyLoss& right)
                  : G4VContinuousDiscreteProcess(right),
     lastMaterial(NULL),
     nmaxCont1(4),
     nmaxCont2(16)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::SetRndmStep(G4bool value) {rndmStepFlag = value;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::SetEnlossFluc(G4bool value) {EnlossFlucFlag = value;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::SetSubSec(G4bool value) {subSecFlag = value;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::SetMinDeltaCutInRange(G4double value)
{MinDeltaCutInRange = value; setMinDeltaCutInRange = true;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::SetStepFunction(G4double c1, G4double c2)
{
 dRoverRange = c1; finalRangeRequested = c2;
 c1lim=dRoverRange;
 c2lim=2.*(1-dRoverRange)*finalRangeRequested;
 c3lim=-(1.-dRoverRange)*finalRangeRequested*finalRangeRequested;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLoss::BuildRangeTable(
        G4PhysicsTable* theDEDXTable,G4PhysicsTable* theRangeTable,
        G4double LowestKineticEnergy,G4double HighestKineticEnergy,G4int TotBin)
// Build range table from the energy loss table
{
   size_t numOfMaterials = theDEDXTable->length();

   if(theRangeTable)
   { theRangeTable->clearAndDestroy();
     delete theRangeTable; }
   theRangeTable = new G4PhysicsTable(numOfMaterials);

   // loop for materials

   for (size_t J=0;  J<numOfMaterials; J++)
   {
     G4PhysicsLogVector* aVector;
     aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                              HighestKineticEnergy,TotBin);

     BuildRangeVector(theDEDXTable,LowestKineticEnergy,HighestKineticEnergy,
                      TotBin,J,aVector);
     theRangeTable->insert(aVector);
   }
   return theRangeTable ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::BuildRangeVector(G4PhysicsTable* theDEDXTable,
                G4double,G4double HighestKineticEnergy,G4int TotBin,
                                 G4int materialIndex,G4PhysicsLogVector* rangeVector)
//  create range vector for a material
{
  G4int nbin=100,i;
  G4bool isOut;

  const G4double small = 1.e-6 ;
  const G4double masslimit = 0.52*MeV ;

  G4double tlim=2.*MeV,t1=0.1*MeV,t2=0.025*MeV ;
  G4double tlime=0.2*keV,factor=2.*electron_mass_c2 ;
  G4double loss1,loss2,ca,cb,cba ;
  G4double losslim,clim ;
  G4double taulim,rangelim,ltaulim,ltaumax,
           LowEdgeEnergy,tau,Value,tau1,sqtau1 ;
  G4double oldValue,tauold ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];

  // cure 'accidental' 0. dE/dx vales first .....
  G4double lossmin = +1.e10 ;
  for (G4int i1=0; i1<TotBin; i1++)
  {
    LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i1) ;
    Value = physicsVector->GetValue(LowEdgeEnergy,isOut);
    if((Value < lossmin)&&(Value > 0.))
      lossmin = Value ;
  }
  for (G4int i2=0; i2<TotBin; i2++)
  {
    LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i2) ;
    Value = physicsVector->GetValue(LowEdgeEnergy,isOut);
    if(Value < lossmin)
      physicsVector->PutValue(i2,small*lossmin) ;
  }
        
  // low energy part first...
 // heavy particle
 if(ParticleMass > masslimit)
 {
  loss1 = physicsVector->GetValue(t1,isOut);
  loss2 = physicsVector->GetValue(t2,isOut);
  tau1 = t1/ParticleMass ;
  sqtau1 = std::sqrt(tau1) ;
  ca = (4.*loss2-loss1)/sqtau1 ;
  cb = (2.*loss1-4.*loss2)/tau1 ;
  cba = cb/ca ;
  taulim = tlim/ParticleMass ;
  ltaulim = std::log(taulim) ;
  ltaumax = std::log(HighestKineticEnergy/ParticleMass) ;

  i=-1;
  oldValue = 0. ;
 
  do
  {
    i += 1 ;
    LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);
    tau = LowEdgeEnergy/ParticleMass;
    if ( tau <= tau1 )
    {
      Value = 2.*ParticleMass*std::log(1.+cba*std::sqrt(tau))/cb ;
    }
    else
    {
      Value = 2.*ParticleMass*std::log(1.+cba*sqtau1)/cb ;
      if(tau<=taulim)
      {
        taulow = tau1 ;
        tauhigh = tau ;
        Value += RangeIntLin(physicsVector,nbin);
      }
      else
      {

        taulow = tau1 ;
        tauhigh = taulim ;
        Value += RangeIntLin(physicsVector,nbin) ;
        ltaulow = ltaulim ;
        ltauhigh = std::log(tau) ;
        Value += RangeIntLog(physicsVector,nbin);
      }
    }

    rangeVector->PutValue(i,Value);
    oldValue = Value ;
    tauold = tau ;
  } while (tau<=taulim) ;
 }
 else
 // electron/positron 
 {
  losslim = physicsVector->GetValue(tlime,isOut) ;

  taulim  = tlime/electron_mass_c2;
  clim    = losslim;
  ltaulim = std::log(taulim);
  ltaumax = std::log(HighestKineticEnergy/electron_mass_c2);

  i=-1;
  oldValue = 0.;

  do
    {
     i += 1 ;
     LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);
     tau = LowEdgeEnergy/electron_mass_c2;
     if (tau <= taulim) Value = factor*std::sqrt(tau*taulim)/clim;
     else {
           rangelim = factor*taulim/losslim ;
           ltaulow  = std::log(taulim);
           ltauhigh = std::log(tau);
           Value    = rangelim+RangeIntLog(physicsVector,nbin);
          }
     rangeVector->PutValue(i,Value);
     oldValue = Value;
     tauold   = tau;

    } while (tau<=taulim);

 }

  i += 1 ;
  for (G4int j=i; j<TotBin; j++)
  {
    LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(j);
    tau = LowEdgeEnergy/ParticleMass;
    ltaulow = std::log(tauold);
    ltauhigh = std::log(tau);
    Value = oldValue+RangeIntLog(physicsVector,nbin);
    rangeVector->PutValue(j,Value);
    oldValue = Value ;

    tauold = tau ;
  }
}   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLoss::RangeIntLin(G4PhysicsVector* physicsVector,
                                    G4int nbin)
//  num. integration, linear binning
{
  G4double dtau,Value,taui,ti,lossi,ci;
  G4bool isOut;
  dtau = (tauhigh-taulow)/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
  {
    taui = taulow + dtau*i ;
    ti = ParticleMass*taui;
    lossi = physicsVector->GetValue(ti,isOut);
    if(i==0)
      ci=0.5;
    else
    {
      if(i<nbin)
        ci=1.;
      else
        ci=0.5;
    }
    Value += ci/lossi;
  }
  Value *= ParticleMass*dtau;
  return Value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLoss::RangeIntLog(G4PhysicsVector* physicsVector,
                                    G4int nbin)
//  num. integration, logarithmic binning
{
  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;
  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
  {
    ui = ltaulow+dltau*i;
    taui = std::exp(ui);
    ti = ParticleMass*taui;
    lossi = physicsVector->GetValue(ti,isOut);
    if(i==0)
      ci=0.5;
    else
    {
      if(i<nbin)
        ci=1.;
      else
        ci=0.5;
    }
    Value += ci*taui/lossi;
  }
  Value *= ParticleMass*dltau;
  return Value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLoss::BuildLabTimeTable(G4PhysicsTable* theDEDXTable,
                                     G4PhysicsTable* theLabTimeTable,
                                     G4double LowestKineticEnergy,
                                     G4double HighestKineticEnergy,G4int TotBin)
                            
{
  size_t numOfMaterials = theDEDXTable->length();

  if(theLabTimeTable)
  { theLabTimeTable->clearAndDestroy();
    delete theLabTimeTable; }
  theLabTimeTable = new G4PhysicsTable(numOfMaterials);


  for (size_t J=0;  J<numOfMaterials; J++)
  {
    G4PhysicsLogVector* aVector;

    aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                            HighestKineticEnergy,TotBin);

    BuildLabTimeVector(theDEDXTable,
              LowestKineticEnergy,HighestKineticEnergy,TotBin,J,aVector);
    theLabTimeTable->insert(aVector);


  }
  return theLabTimeTable ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLoss::BuildProperTimeTable(G4PhysicsTable* theDEDXTable,
                                     G4PhysicsTable* theProperTimeTable,
                                     G4double LowestKineticEnergy,
                                     G4double HighestKineticEnergy,G4int TotBin)

{
  size_t numOfMaterials = theDEDXTable->length();

  if(theProperTimeTable)
  { theProperTimeTable->clearAndDestroy();
    delete theProperTimeTable; }
  theProperTimeTable = new G4PhysicsTable(numOfMaterials);


  for (size_t J=0;  J<numOfMaterials; J++)
  {
    G4PhysicsLogVector* aVector;

    aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                            HighestKineticEnergy,TotBin);

    BuildProperTimeVector(theDEDXTable,
              LowestKineticEnergy,HighestKineticEnergy,TotBin,J,aVector);
    theProperTimeTable->insert(aVector);


  }
  return theProperTimeTable ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::BuildLabTimeVector(G4PhysicsTable* theDEDXTable,
                                    G4double,
                                    G4double HighestKineticEnergy,G4int TotBin,
                                    G4int materialIndex, G4PhysicsLogVector* timeVector)
//  create lab time vector for a material
{

  G4int nbin=100;
  G4bool isOut;
  G4double tlim=5.*keV,parlowen=0.4,ppar=0.5-parlowen ;
  G4double losslim,clim,taulim,timelim,ltaulim,ltaumax,
           LowEdgeEnergy,tau,Value ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];

  // low energy part first...
  losslim = physicsVector->GetValue(tlim,isOut);
  taulim=tlim/ParticleMass ;
  clim=std::sqrt(ParticleMass*tlim/2.)/(c_light*losslim*ppar) ;
  ltaulim = std::log(taulim);
  ltaumax = std::log(HighestKineticEnergy/ParticleMass) ;

  G4int i=-1;
  G4double oldValue = 0. ;
  G4double tauold ;
  do
  {
    i += 1 ;
    LowEdgeEnergy = timeVector->GetLowEdgeEnergy(i);
    tau = LowEdgeEnergy/ParticleMass ;
    if ( tau <= taulim )
    {
      Value = clim*std::exp(ppar*std::log(tau/taulim)) ;
    }
    else
    {
      timelim=clim ;
      ltaulow = std::log(taulim);
      ltauhigh = std::log(tau);
      Value = timelim+LabTimeIntLog(physicsVector,nbin);
    }
    timeVector->PutValue(i,Value);
    oldValue = Value ;
    tauold = tau ;
  } while (tau<=taulim) ;
  i += 1 ;
  for (G4int j=i; j<TotBin; j++)
  {
    LowEdgeEnergy = timeVector->GetLowEdgeEnergy(j);
    tau = LowEdgeEnergy/ParticleMass ;
    ltaulow = std::log(tauold);
    ltauhigh = std::log(tau);
    Value = oldValue+LabTimeIntLog(physicsVector,nbin);
    timeVector->PutValue(j,Value);
    oldValue = Value ;
    tauold = tau ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::BuildProperTimeVector(G4PhysicsTable* theDEDXTable,
                                    G4double,
                                    G4double HighestKineticEnergy,G4int TotBin,
                                    G4int materialIndex, G4PhysicsLogVector* timeVector)
//  create proper time vector for a material
{
  G4int nbin=100;
  G4bool isOut;
  G4double tlim=5.*keV,parlowen=0.4,ppar=0.5-parlowen ;
  G4double losslim,clim,taulim,timelim,ltaulim,ltaumax,
           LowEdgeEnergy,tau,Value ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];

  // low energy part first...
  losslim = physicsVector->GetValue(tlim,isOut);
  taulim=tlim/ParticleMass ;
  clim=std::sqrt(ParticleMass*tlim/2.)/(c_light*losslim*ppar) ;
  ltaulim = std::log(taulim);
  ltaumax = std::log(HighestKineticEnergy/ParticleMass) ;

  G4int i=-1;
  G4double oldValue = 0. ;
  G4double tauold ;
  do
  {
    i += 1 ;
    LowEdgeEnergy = timeVector->GetLowEdgeEnergy(i);
    tau = LowEdgeEnergy/ParticleMass ;
    if ( tau <= taulim )
    {
      Value = clim*std::exp(ppar*std::log(tau/taulim)) ;
    }
    else
    {
      timelim=clim ;
      ltaulow = std::log(taulim);
      ltauhigh = std::log(tau);
      Value = timelim+ProperTimeIntLog(physicsVector,nbin);
    }
    timeVector->PutValue(i,Value);
    oldValue = Value ;
    tauold = tau ;
  } while (tau<=taulim) ;
  i += 1 ;
  for (G4int j=i; j<TotBin; j++)
  {
    LowEdgeEnergy = timeVector->GetLowEdgeEnergy(j);
    tau = LowEdgeEnergy/ParticleMass ;
    ltaulow = std::log(tauold);
    ltauhigh = std::log(tau);
    Value = oldValue+ProperTimeIntLog(physicsVector,nbin);
    timeVector->PutValue(j,Value);
    oldValue = Value ;
    tauold = tau ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLoss::LabTimeIntLog(G4PhysicsVector* physicsVector,
                                    G4int nbin)
//  num. integration, logarithmic binning
{
  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;
  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
  {
    ui = ltaulow+dltau*i;
    taui = std::exp(ui);
    ti = ParticleMass*taui;
    lossi = physicsVector->GetValue(ti,isOut);
    if(i==0)
      ci=0.5;
    else
    {
      if(i<nbin)
        ci=1.;
      else
        ci=0.5;
    }
    Value += ci*taui*(ti+ParticleMass)/(std::sqrt(ti*(ti+2.*ParticleMass))*lossi);
  }
  Value *= ParticleMass*dltau/c_light;
  return Value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLoss::ProperTimeIntLog(G4PhysicsVector* physicsVector,
                                    G4int nbin)
//  num. integration, logarithmic binning
{
  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;
  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
  {
    ui = ltaulow+dltau*i;
    taui = std::exp(ui);
    ti = ParticleMass*taui;
    lossi = physicsVector->GetValue(ti,isOut);
    if(i==0)
      ci=0.5;
    else
    {
      if(i<nbin)
        ci=1.;
      else
        ci=0.5;
    }
    Value += ci*taui*ParticleMass/(std::sqrt(ti*(ti+2.*ParticleMass))*lossi);
  }
  Value *= ParticleMass*dltau/c_light;
  return Value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLoss::BuildInverseRangeTable(G4PhysicsTable* theRangeTable,
                                   G4PhysicsTable* theRangeCoeffATable,
                                   G4PhysicsTable* theRangeCoeffBTable,
                                   G4PhysicsTable* theRangeCoeffCTable,
                                   G4PhysicsTable* theInverseRangeTable,
                                   G4double LowestKineticEnergy,
                                   G4double HighestKineticEnergy,G4int TotBin)
// Build inverse table of the range table
{
  G4double SmallestRange,BiggestRange ;
  G4bool isOut ;
  size_t numOfMaterials = theRangeTable->length();

    if(theInverseRangeTable)
    { theInverseRangeTable->clearAndDestroy();
      delete theInverseRangeTable; }
    theInverseRangeTable = new G4PhysicsTable(numOfMaterials);

  // loop for materials
  for (size_t J=0;  J<numOfMaterials; J++)
  {
    SmallestRange = (*theRangeTable)(J)->
                       GetValue(LowestKineticEnergy,isOut) ;
    BiggestRange = (*theRangeTable)(J)->
                       GetValue(HighestKineticEnergy,isOut) ;
    G4PhysicsLogVector* aVector;
    aVector = new G4PhysicsLogVector(SmallestRange,
                            BiggestRange,TotBin);

    InvertRangeVector(theRangeTable,
                      theRangeCoeffATable,
                      theRangeCoeffBTable,
                      theRangeCoeffCTable,
         LowestKineticEnergy,HighestKineticEnergy,TotBin,J, aVector);

    theInverseRangeTable->insert(aVector);
  }
  return theInverseRangeTable ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLoss::InvertRangeVector(G4PhysicsTable* theRangeTable,
                              G4PhysicsTable* theRangeCoeffATable,
                              G4PhysicsTable* theRangeCoeffBTable,
                              G4PhysicsTable* theRangeCoeffCTable,
                              G4double LowestKineticEnergy,
                              G4double HighestKineticEnergy,G4int TotBin,
                              G4int  materialIndex,G4PhysicsLogVector* aVector)
//  invert range vector for a material
{
  G4double LowEdgeRange,A,B,C,discr,KineticEnergy ;
  G4double RTable = std::exp(std::log(HighestKineticEnergy/LowestKineticEnergy)/TotBin) ;
  G4double Tbin = LowestKineticEnergy/RTable ;
  G4double rangebin = 0.0 ;
  G4int binnumber = -1 ;
  G4bool isOut ;

  //loop for range values
  for( G4int i=0; i<TotBin; i++)
  {
    LowEdgeRange = aVector->GetLowEdgeEnergy(i) ;  //i.e. GetLowEdgeValue(i)
    if( rangebin < LowEdgeRange )
    {
      do
      {
        binnumber += 1 ;
        Tbin *= RTable ;
        rangebin = (*theRangeTable)(materialIndex)->GetValue(Tbin,isOut) ;
      }
      while ((rangebin < LowEdgeRange) && (binnumber < TotBin )) ;
    }

    if(binnumber == 0)
      KineticEnergy = LowestKineticEnergy ;
    else if(binnumber == TotBin-1)
      KineticEnergy = HighestKineticEnergy ;
    else
    {
      A = (*(*theRangeCoeffATable)(materialIndex))(binnumber-1) ;
      B = (*(*theRangeCoeffBTable)(materialIndex))(binnumber-1) ;
      C = (*(*theRangeCoeffCTable)(materialIndex))(binnumber-1) ;
      if(A==0.)
         KineticEnergy = (LowEdgeRange -C )/B ;
      else
      {
         discr = B*B - 4.*A*(C-LowEdgeRange);
         discr = discr>0. ? std::sqrt(discr) : 0.;
         KineticEnergy = 0.5*(discr-B)/A ;
      }
    }

    aVector->PutValue(i,KineticEnergy) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLoss::BuildRangeCoeffATable(G4PhysicsTable* theRangeTable,
                                      G4PhysicsTable* theRangeCoeffATable,
                                      G4double LowestKineticEnergy,
                             G4double HighestKineticEnergy,G4int TotBin)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "A"
{
  G4int numOfMaterials = theRangeTable->length();

  if(theRangeCoeffATable)
  { theRangeCoeffATable->clearAndDestroy();
    delete theRangeCoeffATable; }
  theRangeCoeffATable = new G4PhysicsTable(numOfMaterials);

  G4double RTable = std::exp(std::log(HighestKineticEnergy/LowestKineticEnergy)/TotBin) ;
  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  G4double w1 = RTable/w , w2 = -RTable*R1/w , w3 = R2/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4int binmax=TotBin ;
    G4PhysicsLinearVector* aVector =
                           new G4PhysicsLinearVector(0.,binmax, TotBin);
    Ti = LowestKineticEnergy ;
    G4PhysicsVector* rangeVector= (*theRangeTable)[J];

    for ( G4int i=0; i<TotBin; i++)
    {
      Ri = rangeVector->GetValue(Ti,isOut) ;
      if ( i==0 )
        Rim = 0. ;
      else
      {
        Tim = Ti/RTable ;
        Rim = rangeVector->GetValue(Tim,isOut);
      }
      if ( i==(TotBin-1))
        Rip = Ri ;
      else
      {
        Tip = Ti*RTable ;
        Rip = rangeVector->GetValue(Tip,isOut);
      }
      Value = (w1*Rip + w2*Ri + w3*Rim)/(Ti*Ti) ;

      aVector->PutValue(i,Value);
      Ti = RTable*Ti ;
    }
 
    theRangeCoeffATable->insert(aVector);
  }
  return theRangeCoeffATable ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4PhysicsTable* G4VEnergyLoss::BuildRangeCoeffBTable(G4PhysicsTable* theRangeTable,
                                      G4PhysicsTable* theRangeCoeffBTable,
                                      G4double LowestKineticEnergy,
                             G4double HighestKineticEnergy,G4int TotBin)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "B"
{
  G4int numOfMaterials = theRangeTable->length();

  if(theRangeCoeffBTable)
  { theRangeCoeffBTable->clearAndDestroy();
    delete theRangeCoeffBTable; }
  theRangeCoeffBTable = new G4PhysicsTable(numOfMaterials);

  G4double RTable = std::exp(std::log(HighestKineticEnergy/LowestKineticEnergy)/TotBin) ;
  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  G4double w1 = -R1/w , w2 = R1*(R2+1.)/w , w3 = -R2*R1/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4int binmax=TotBin ;
    G4PhysicsLinearVector* aVector =
                        new G4PhysicsLinearVector(0.,binmax, TotBin);
    Ti = LowestKineticEnergy ;
    G4PhysicsVector* rangeVector= (*theRangeTable)[J];
  
    for ( G4int i=0; i<TotBin; i++)
    {
      Ri = rangeVector->GetValue(Ti,isOut) ;
      if ( i==0 )
         Rim = 0. ;
      else
      {
        Tim = Ti/RTable ;
        Rim = rangeVector->GetValue(Tim,isOut);
      }
      if ( i==(TotBin-1))
        Rip = Ri ;
      else
      {
        Tip = Ti*RTable ;
        Rip = rangeVector->GetValue(Tip,isOut);
      }
      Value = (w1*Rip + w2*Ri + w3*Rim)/Ti;

      aVector->PutValue(i,Value);
      Ti = RTable*Ti ;
    }
    theRangeCoeffBTable->insert(aVector);
  }
  return theRangeCoeffBTable ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLoss::BuildRangeCoeffCTable(G4PhysicsTable* theRangeTable,
                                      G4PhysicsTable* theRangeCoeffCTable,
                                      G4double LowestKineticEnergy,
                             G4double HighestKineticEnergy,G4int TotBin)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "C"
{
  G4int numOfMaterials = theRangeTable->length();

  if(theRangeCoeffCTable)
  { theRangeCoeffCTable->clearAndDestroy();
    delete theRangeCoeffCTable; }
  theRangeCoeffCTable = new G4PhysicsTable(numOfMaterials);

  G4double RTable = std::exp(std::log(HighestKineticEnergy/LowestKineticEnergy)/TotBin) ;
  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  G4double w1 = 1./w , w2 = -RTable*R1/w , w3 = RTable*R2/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4int binmax=TotBin ;
    G4PhysicsLinearVector* aVector =
                      new G4PhysicsLinearVector(0.,binmax, TotBin);
    Ti = LowestKineticEnergy ;
    G4PhysicsVector* rangeVector= (*theRangeTable)[J];
  
    for ( G4int i=0; i<TotBin; i++)
    {
      Ri = rangeVector->GetValue(Ti,isOut) ;
      if ( i==0 )
        Rim = 0. ;
      else
      {
        Tim = Ti/RTable ;
        Rim = rangeVector->GetValue(Tim,isOut);
      }
      if ( i==(TotBin-1))
        Rip = Ri ;
      else
      {
        Tip = Ti*RTable ;
        Rip = rangeVector->GetValue(Tip,isOut);
      }
      Value = w1*Rip + w2*Ri + w3*Rim ;

      aVector->PutValue(i,Value);
      Ti = RTable*Ti ;
    }
    theRangeCoeffCTable->insert(aVector);
  }
  return theRangeCoeffCTable ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLoss::GetLossWithFluct(const G4DynamicParticle* aParticle,
                                         const G4MaterialCutsCouple* couple,
                                               G4double ChargeSquare,
                                               G4double    MeanLoss,
                                               G4double step )
//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is essentially the same as in Glandz in Geant3.
{
  const G4double minLoss = 1.*eV ;
  const G4double probLim = 0.01 ;
  const G4double sumaLim = -std::log(probLim) ;
  const G4double alim=10.;
  const G4double kappa = 10. ;
  const G4double factor = twopi_mc2_rcl2 ;
  const G4Material* aMaterial = couple->GetMaterial();

  // check if the material has changed ( cache mechanism)

  if (aMaterial != lastMaterial)
    {
      lastMaterial = aMaterial;
      imat         = couple->GetIndex();
      f1Fluct      = aMaterial->GetIonisation()->GetF1fluct();
      f2Fluct      = aMaterial->GetIonisation()->GetF2fluct();
      e1Fluct      = aMaterial->GetIonisation()->GetEnergy1fluct();
      e2Fluct      = aMaterial->GetIonisation()->GetEnergy2fluct();
      e1LogFluct   = aMaterial->GetIonisation()->GetLogEnergy1fluct();
      e2LogFluct   = aMaterial->GetIonisation()->GetLogEnergy2fluct();
      rateFluct    = aMaterial->GetIonisation()->GetRateionexcfluct();
      ipotFluct    = aMaterial->GetIonisation()->GetMeanExcitationEnergy();
      ipotLogFluct = aMaterial->GetIonisation()->GetLogMeanExcEnergy();
    }
  G4double threshold,w1,w2,C,
           beta2,suma,e0,loss,lossc ,w,electronDensity;
  G4double a1,a2,a3;
  G4double p1,p2,p3 ;
  G4int nb;
  G4double Corrfac, na,alfa,rfac,namean,sa,alfa1,ea,sea;
  G4double dp3;
  G4double siga ;

  // shortcut for very very small loss
  if(MeanLoss < minLoss) return MeanLoss ;

  // get particle data
  G4double Tkin   = aParticle->GetKineticEnergy();
  ParticleMass = aParticle->GetMass() ;

  threshold = (*((G4ProductionCutsTable::GetProductionCutsTable())
                ->GetEnergyCutsVector(1)))[imat];
  G4double rmass = electron_mass_c2/ParticleMass;
  G4double tau   = Tkin/ParticleMass, tau1 = tau+1., tau2 = tau*(tau+2.);
  G4double Tm    = 2.*electron_mass_c2*tau2/(1.+2.*tau1*rmass+rmass*rmass);

  if(Tm > threshold) Tm = threshold;

  beta2 = tau2/(tau1*tau1);

  // Gaussian fluctuation ?
  if(MeanLoss >= kappa*Tm || Tm < kappa*ipotFluct)
  {
    electronDensity = aMaterial->GetElectronDensity() ;
    siga = std::sqrt(Tm*(1.0-0.5*beta2)*step*
                factor*electronDensity*ChargeSquare/beta2) ;
    do {
      loss = G4RandGauss::shoot(MeanLoss,siga) ;
    } while (loss < 0. || loss > 2.*MeanLoss);
    return loss;
  }

  w1 = Tm/ipotFluct;
  w2 = std::log(2.*electron_mass_c2*tau2);

  C = MeanLoss*(1.-rateFluct)/(w2-ipotLogFluct-beta2);

  a1 = C*f1Fluct*(w2-e1LogFluct-beta2)/e1Fluct;
  a2 = C*f2Fluct*(w2-e2LogFluct-beta2)/e2Fluct;
  a3 = rateFluct*MeanLoss*(Tm-ipotFluct)/(ipotFluct*Tm*std::log(w1));
  if(a1 < 0.) a1 = 0.;
  if(a2 < 0.) a2 = 0.;
  if(a3 < 0.) a3 = 0.;

  suma = a1+a2+a3;

  loss = 0. ;

  if(suma < sumaLim)             // very small Step
    {
      e0 = aMaterial->GetIonisation()->GetEnergy0fluct();

      if(Tm == ipotFluct)
      {
        a3 = MeanLoss/e0;

        if(a3>alim)
        {
          siga=std::sqrt(a3) ;
          p3 = std::max(0.,G4RandGauss::shoot(a3,siga)+0.5);
        }
        else
          p3 = G4float(G4Poisson(a3));

        loss = p3*e0 ;

        if(p3 > 0)
          loss += (1.-2.*G4UniformRand())*e0 ;

      }
      else
      {
        Tm = Tm-ipotFluct+e0 ;
        a3 = MeanLoss*(Tm-e0)/(Tm*e0*std::log(Tm/e0));

        if(a3>alim)
        {
          siga=std::sqrt(a3) ;
          p3 = std::max(0.,G4RandGauss::shoot(a3,siga)+0.5);
        }
        else
          p3 = G4float(G4Poisson(a3));

        if(p3 > 0)
        {
          w = (Tm-e0)/Tm ;
          if(p3 > G4float(nmaxCont2))
          {
            dp3 = p3 ;
            Corrfac = dp3/G4float(nmaxCont2) ;
            p3 = G4float(nmaxCont2) ;
          }
          else
            Corrfac = 1. ;

          for(G4int i=0; i<G4int(p3); i++) loss += 1./(1.-w*G4UniformRand()) ;
          loss *= e0*Corrfac ;  
        }        
      }
    }
    
  else                              // not so small Step
    {
      // excitation type 1
      if(a1>alim)
      {
        siga=std::sqrt(a1) ;
        p1 = std::max(0.,G4RandGauss::shoot(a1,siga)+0.5);
      }
      else
       p1 = G4float(G4Poisson(a1));

      // excitation type 2
      if(a2>alim)
      {
        siga=std::sqrt(a2) ;
        p2 = std::max(0.,G4RandGauss::shoot(a2,siga)+0.5);
      }
      else
        p2 = G4float(G4Poisson(a2));

      loss = p1*e1Fluct+p2*e2Fluct;
      // smearing to avoid unphysical peaks
      if(p2 > 0)
        loss += (1.-2.*G4UniformRand())*e2Fluct;
      else if (loss>0.)
        loss += (1.-2.*G4UniformRand())*e1Fluct;   

      // ionisation .......................................
     if(a3 > 0.)
     {
      if(a3>alim)
      {
        siga=std::sqrt(a3) ;
        p3 = std::max(0.,G4RandGauss::shoot(a3,siga)+0.5);
      }
      else
        p3 = G4float(G4Poisson(a3));

      lossc = 0.;
      if(p3 > 0)
      {
        na = 0.; 
        alfa = 1.;
        if (p3 > G4float(nmaxCont2))
        {
          dp3        = p3;
          rfac       = dp3/(G4float(nmaxCont2)+dp3);
          namean     = p3*rfac;
          sa         = G4float(nmaxCont1)*rfac;
          na         = G4RandGauss::shoot(namean,sa);
          if (na > 0.)
          {
            alfa   = w1*(G4float(nmaxCont2)+p3)/(w1*G4float(nmaxCont2)+p3);
            alfa1  = alfa*std::log(alfa)/(alfa-1.);
            ea     = na*ipotFluct*alfa1;
            sea    = ipotFluct*std::sqrt(na*(alfa-alfa1*alfa1));
            lossc += G4RandGauss::shoot(ea,sea);
          }
        }

        nb = G4int(p3-na);
        if (nb > 0)
        {
          w2 = alfa*ipotFluct;
          w  = (Tm-w2)/Tm;      
          for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
        }
      }        
      loss += lossc;  
     }
    } 

  if( loss < 0.)
    loss = 0.;

  return loss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
   
G4bool G4VEnergyLoss::EqualCutVectors( G4double* vec1, G4double* vec2 )
{
  if ( (vec1==0 ) || (vec2==0) ) return false;
  
  G4bool flag = true;
   
  for (size_t j=0; flag && j<G4Material::GetNumberOfMaterials(); j++){
    flag = (vec1[j] == vec2[j]);
  }

  return flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double* G4VEnergyLoss::CopyCutVectors( G4double* dest, G4double* source )
{
  if ( dest != 0) delete [] dest;
  dest = new G4double [G4Material::GetNumberOfMaterials()];
  for (size_t j=0; j<G4Material::GetNumberOfMaterials(); j++){
    dest[j] = source[j];
  }
  return dest;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VEnergyLoss::CutsWhereModified()
{
  G4bool wasModified = false;
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for (size_t j=0; j<numOfCouples; j++){
    if (theCoupleTable->GetMaterialCutsCouple(j)->IsRecalcNeeded()) {
      wasModified = true;
      break;
    }
  }
  return wasModified;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

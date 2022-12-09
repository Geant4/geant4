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
// -------------------------------------------------------------------
//
// GEANT4 Class
// File name:     G4PAIPhotData.cc
//
// Author:        V.Grichine based on G4PAIModelData MT
//
// Creation date: 07.10.2013
//
// Modifications:
//

#include "G4PAIPhotData.hh"
#include "G4PAIPhotModel.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SandiaTable.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"

////////////////////////////////////////////////////////////////////////

using namespace std;

G4PAIPhotData::G4PAIPhotData(G4double tmin, G4double tmax, G4int ver)
{ 
  const G4int nPerDecade     = 10; 
  const G4double lowestTkin  = 50*keV;
  const G4double highestTkin = 10*TeV;

  // fPAIxSection.SetVerbose(ver);

  fLowestKineticEnergy  = std::max(tmin, lowestTkin);
  fHighestKineticEnergy = tmax;

  if(tmax < 10*fLowestKineticEnergy) 
  { 
    fHighestKineticEnergy = 10*fLowestKineticEnergy;
  } 
  else if(tmax > highestTkin) 
  {
    fHighestKineticEnergy = std::max(highestTkin, 10*fLowestKineticEnergy);
  }
  fTotBin = (G4int)(nPerDecade*
		    std::log10(fHighestKineticEnergy/fLowestKineticEnergy));

  fParticleEnergyVector = new G4PhysicsLogVector(fLowestKineticEnergy,
						 fHighestKineticEnergy,
						 fTotBin);
  if(0 < ver) {
    G4cout << "### G4PAIPhotData: Nbins= " << fTotBin
	   << " Tmin(MeV)= " << fLowestKineticEnergy/MeV
	   << " Tmax(GeV)= " << fHighestKineticEnergy/GeV 
	   << "  tmin(keV)= " << tmin/keV << G4endl;
  }
}

////////////////////////////////////////////////////////////////////////////

G4PAIPhotData::~G4PAIPhotData()
{
  //G4cout << "G4PAIPhotData::~G4PAIPhotData() " << this << G4endl;
  std::size_t n = fPAIxscBank.size();
  if(0 < n) 
  {
    for(std::size_t i=0; i<n; ++i) 
    {
      if(fPAIxscBank[i]) 
      {
	fPAIxscBank[i]->clearAndDestroy();
	delete fPAIxscBank[i];
	fPAIxscBank[i] = nullptr;
      }
      if(fPAIdEdxBank[i]) 
      {
	fPAIdEdxBank[i]->clearAndDestroy();
	delete fPAIdEdxBank[i];
	fPAIdEdxBank[i] = nullptr;
      }
      delete fdEdxTable[i];
      delete fdNdxCutTable[i];
      fdEdxTable[i] = nullptr;
      fdNdxCutTable[i] = nullptr;
    }
  }
  delete fParticleEnergyVector;
  fParticleEnergyVector = nullptr;
  //G4cout << "G4PAIPhotData::~G4PAIPhotData() done for " << this << G4endl;  
}

///////////////////////////////////////////////////////////////////////////////

void G4PAIPhotData::Initialise(const G4MaterialCutsCouple* couple,
                                G4double cut, G4PAIPhotModel* model)
{
  G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
  G4int jMatCC;

  for (jMatCC = 0; jMatCC < numOfCouples; ++jMatCC )
  {
    if( couple == theCoupleTable->GetMaterialCutsCouple(jMatCC) ) break;
  }
  if( jMatCC == numOfCouples && jMatCC > 0 ) --jMatCC;

  const vector<G4double>*  deltaCutInKineticEnergy = theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut);
  const vector<G4double>*  photonCutInKineticEnergy = theCoupleTable->GetEnergyCutsVector(idxG4GammaCut);
  G4double deltaCutInKineticEnergyNow  = (*deltaCutInKineticEnergy)[jMatCC];
  G4double photonCutInKineticEnergyNow = (*photonCutInKineticEnergy)[jMatCC];

  G4cout<<"G4PAIPhotData::Initialise: "<<"cut = "<<cut/keV<<" keV; cutEl = "
        <<deltaCutInKineticEnergyNow/keV<<" keV; cutPh = "
	<<photonCutInKineticEnergyNow/keV<<" keV"<<G4endl;

  // if( deltaCutInKineticEnergyNow != cut ) deltaCutInKineticEnergyNow = cut; // exception??

  auto dEdxCutVector =
    new G4PhysicsLogVector(fLowestKineticEnergy,
			   fHighestKineticEnergy,
			   fTotBin);

  auto dNdxCutVector =
    new G4PhysicsLogVector(fLowestKineticEnergy,
			   fHighestKineticEnergy,
			   fTotBin);
  auto dNdxCutPhotonVector =
    new G4PhysicsLogVector(fLowestKineticEnergy,
			   fHighestKineticEnergy,
			   fTotBin);
  auto dNdxCutPlasmonVector =
    new G4PhysicsLogVector(fLowestKineticEnergy,
			   fHighestKineticEnergy,
			   fTotBin);

  const G4Material* mat = couple->GetMaterial();     
  fSandia.Initialize(const_cast<G4Material*>(mat));

  auto PAItransferTable = new G4PhysicsTable(fTotBin+1);
  auto PAIphotonTable = new G4PhysicsTable(fTotBin+1);
  auto PAIplasmonTable = new G4PhysicsTable(fTotBin+1);

  auto PAIdEdxTable = new G4PhysicsTable(fTotBin+1);
  auto dEdxMeanVector =
    new G4PhysicsLogVector(fLowestKineticEnergy,
			   fHighestKineticEnergy,
			   fTotBin);

  // low energy Sandia interval
  G4double Tmin = fSandia.GetSandiaMatTablePAI(0,0); 

  // energy safety
  const G4double deltaLow = 100.*eV; 

  for (G4int i = 0; i <= fTotBin; ++i) 
  {
    G4double kinEnergy = fParticleEnergyVector->Energy(i);
    G4double Tmax = model->ComputeMaxEnergy(kinEnergy);
    G4double tau = kinEnergy/proton_mass_c2;
    G4double bg2 = tau*( tau + 2. );

    if ( Tmax < Tmin + deltaLow ) Tmax = Tmin + deltaLow; 

    fPAIxSection.Initialize( mat, Tmax, bg2, &fSandia);

    //G4cout << i << ". TransferMax(keV)= "<< Tmax/keV << "  cut(keV)= " 
    //	   << cut/keV << "  E(MeV)= " << kinEnergy/MeV << G4endl;

    G4int n = fPAIxSection.GetSplineSize();

    auto transferVector = new G4PhysicsFreeVector(n);
    auto photonVector   = new G4PhysicsFreeVector(n);
    auto plasmonVector  = new G4PhysicsFreeVector(n);

    auto dEdxVector     = new G4PhysicsFreeVector(n);

    for( G4int k = 0; k < n; k++ )
    {
      G4double t = fPAIxSection.GetSplineEnergy(k+1);

      transferVector->PutValue(k , t, 
                               t*fPAIxSection.GetIntegralPAIxSection(k+1));
      photonVector->PutValue(k , t, 
                               t*fPAIxSection.GetIntegralCerenkov(k+1));
      plasmonVector->PutValue(k , t, 
                               t*fPAIxSection.GetIntegralPlasmon(k+1));

      dEdxVector->PutValue(k, t, fPAIxSection.GetIntegralPAIdEdx(k+1));
    }
    // G4cout << *transferVector << G4endl;

    G4double ionloss = fPAIxSection.GetMeanEnergyLoss();//  total <dE/dx>

    if(ionloss < 0.0) ionloss = 0.0; 

    dEdxMeanVector->PutValue(i,ionloss);

    G4double dNdxCut = transferVector->Value(deltaCutInKineticEnergyNow)/deltaCutInKineticEnergyNow;
    G4double dNdxCutPhoton = photonVector->Value(photonCutInKineticEnergyNow)/photonCutInKineticEnergyNow;
    G4double dNdxCutPlasmon = plasmonVector->Value(deltaCutInKineticEnergyNow)/deltaCutInKineticEnergyNow;

    G4double dEdxCut = dEdxVector->Value(cut)/cut;
    //G4cout << "i= " << i << " x= " << dNdxCut << G4endl;

    if(dNdxCut < 0.0) { dNdxCut = 0.0; }
    if(dNdxCutPhoton < 0.0) { dNdxCutPhoton = 0.0; }
    if(dNdxCutPlasmon < 0.0) { dNdxCutPlasmon = 0.0; }

    dNdxCutVector->PutValue(i, dNdxCut);
    dNdxCutPhotonVector->PutValue(i, dNdxCutPhoton);
    dNdxCutPlasmonVector->PutValue(i, dNdxCutPlasmon);

    dEdxCutVector->PutValue(i, dEdxCut);

    PAItransferTable->insertAt(i,transferVector);
    PAIphotonTable->insertAt(i,photonVector);
    PAIplasmonTable->insertAt(i,plasmonVector);
    PAIdEdxTable->insertAt(i,dEdxVector);

  } // end of Tkin loop

  fPAIxscBank.push_back(PAItransferTable);
  fPAIphotonBank.push_back(PAIphotonTable);
  fPAIplasmonBank.push_back(PAIplasmonTable);

  fPAIdEdxBank.push_back(PAIdEdxTable);
  fdEdxTable.push_back(dEdxMeanVector);

  fdNdxCutTable.push_back(dNdxCutVector);
  fdNdxCutPhotonTable.push_back(dNdxCutPhotonVector);
  fdNdxCutPlasmonTable.push_back(dNdxCutPlasmonVector);

  fdEdxCutTable.push_back(dEdxCutVector);
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIPhotData::DEDXPerVolume(G4int coupleIndex, G4double scaledTkin,
			 G4double cut) const
{
  // VI: iPlace is the low edge index of the bin
  // iPlace is in interval from 0 to (N-1)
  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

  G4bool one = true;
  if(scaledTkin >= fParticleEnergyVector->Energy(nPlace)) { iPlace = nPlace; }
  else if(scaledTkin > fParticleEnergyVector->Energy(0)) { 
    one = false; 
  }

  // VI: apply interpolation of the vector
  G4double dEdx = fdEdxTable[coupleIndex]->Value(scaledTkin);
  G4double del  = (*(fPAIdEdxBank[coupleIndex]))(iPlace)->Value(cut);
  if(!one) {
    G4double del2 = (*(fPAIdEdxBank[coupleIndex]))(iPlace+1)->Value(cut);
    G4double E1 = fParticleEnergyVector->Energy(iPlace); 
    G4double E2 = fParticleEnergyVector->Energy(iPlace+1);
    G4double W  = 1.0/(E2 - E1);
    G4double W1 = (E2 - scaledTkin)*W;
    G4double W2 = (scaledTkin - E1)*W;
    del *= W1;
    del += W2*del2;
  }
  dEdx -= del;

  if( dEdx < 0.) { dEdx = 0.; }
  return dEdx;
}

/////////////////////////////////////////////////////////////////////////

G4double 
G4PAIPhotData::CrossSectionPerVolume(G4int coupleIndex, 
				      G4double scaledTkin,
				      G4double tcut, G4double tmax) const
{
  G4double cross, xscEl, xscEl2, xscPh, xscPh2;

  cross=tcut+tmax;

  // iPlace is in interval from 0 to (N-1)

  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

  G4bool one = true;

  if(      scaledTkin >= fParticleEnergyVector->Energy(nPlace))  iPlace = nPlace; 
  else if( scaledTkin > fParticleEnergyVector->Energy(0)      )   one   = false; 
  

  xscEl2 = (*fdNdxCutPlasmonTable[coupleIndex])(iPlace);
  xscPh2 = (*fdNdxCutPhotonTable[coupleIndex])(iPlace);

  xscPh = xscPh2;
  xscEl = xscEl2;

  cross  = xscPh + xscEl;
 
  if( !one ) 
  {
    xscEl2 = (*fdNdxCutPlasmonTable[coupleIndex])(iPlace+1);

    G4double E1 = fParticleEnergyVector->Energy(iPlace); 
    G4double E2 = fParticleEnergyVector->Energy(iPlace+1);

    G4double W  = 1.0/(E2 - E1);
    G4double W1 = (E2 - scaledTkin)*W;
    G4double W2 = (scaledTkin - E1)*W;

    xscEl *= W1;
    xscEl += W2*xscEl2;

    xscPh2 = (*fdNdxCutPhotonTable[coupleIndex])(iPlace+1);

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);

    W  = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;

    xscPh *= W1;
    xscPh += W2*xscPh2;

    cross = xscEl + xscPh;
  }
  if( cross < 0.0)  cross = 0.0; 

  return cross;
}

/////////////////////////////////////////////////////////////////////////

G4double 
G4PAIPhotData::GetPlasmonRatio(G4int coupleIndex, G4double scaledTkin) const
{
  G4double cross, xscEl, xscEl2, xscPh, xscPh2, plRatio;
  // iPlace is in interval from 0 to (N-1)

  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

  G4bool one = true;

  if(      scaledTkin >= fParticleEnergyVector->Energy(nPlace))  iPlace = nPlace; 
  else if( scaledTkin > fParticleEnergyVector->Energy(0)      )   one   = false; 
  

  xscEl2 = (*fdNdxCutPlasmonTable[coupleIndex])(iPlace);
  xscPh2 = (*fdNdxCutPhotonTable[coupleIndex])(iPlace);

  xscPh = xscPh2;
  xscEl = xscEl2;

  cross  = xscPh + xscEl;
 
  if( !one ) 
  {
    xscEl2 = (*fdNdxCutPlasmonTable[coupleIndex])(iPlace+1);

    G4double E1 = fParticleEnergyVector->Energy(iPlace); 
    G4double E2 = fParticleEnergyVector->Energy(iPlace+1);

    G4double W  = 1.0/(E2 - E1);
    G4double W1 = (E2 - scaledTkin)*W;
    G4double W2 = (scaledTkin - E1)*W;

    xscEl *= W1;
    xscEl += W2*xscEl2;

    xscPh2 = (*fdNdxCutPhotonTable[coupleIndex])(iPlace+1);

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);

    W  = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;

    xscPh *= W1;
    xscPh += W2*xscPh2;

    cross = xscEl + xscPh;
  }
  if( cross <= 0.0)  
  {
    plRatio = 2.0; 
  }
  else
  {
    plRatio = xscEl/cross;

    if( plRatio > 1. || plRatio < 0.) plRatio = 2.0;
  }
  return plRatio;
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIPhotData::SampleAlongStepTransfer(G4int coupleIndex, 
                                                 G4double kinEnergy,
						 G4double scaledTkin,
						 G4double stepFactor) const
{
  G4double loss = 0.0;
  G4double omega; 
  G4double position, E1, E2, W1, W2, W, dNdxCut1, dNdxCut2, meanNumber;

  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;
 
  G4bool one = true;

  if     (scaledTkin >= fParticleEnergyVector->Energy(nPlace)) iPlace = nPlace; 
  else if(scaledTkin > fParticleEnergyVector->Energy(0))          one = false; 

  G4PhysicsLogVector* vcut = fdNdxCutTable[coupleIndex];
  G4PhysicsVector*      v1 = (*(fPAIxscBank[coupleIndex]))(iPlace);
  G4PhysicsVector*      v2 = nullptr;

  dNdxCut1    = (*vcut)[iPlace];
  G4double e1 = v1->Energy(0);
  G4double e2 = e1;

  G4double meanN1 = ((*v1)[0]/e1 - dNdxCut1)*stepFactor;

  meanNumber = meanN1;

  // G4cout<<"iPlace = "<<iPlace<< " meanN1= " << meanN1 
  //	<<"    (*v1)[0]/e1 = "<<(*v1)[0]/e1<< " dNdxCut1= " << dNdxCut1 << G4endl;

  dNdxCut2 = dNdxCut1;
  W1 = 1.0;
  W2 = 0.0;
  if(!one) 
  {
    v2 = (*(fPAIxscBank[coupleIndex]))(iPlace+1);
    dNdxCut2 = (*vcut)[iPlace+1];
    e2 = v2->Energy(0);

    G4double meanN2 = ((*v2)[0]/e2 - dNdxCut2)*stepFactor;

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;
    meanNumber = W1*meanN1 + W2*meanN2;

    //G4cout<<"meanN= " <<  meanNumber << " meanN2= " << meanN2 
    //	  << " dNdxCut2= " << dNdxCut2 << G4endl;
  }
  if( meanNumber <= 0.0) return 0.0; 

  G4int numOfCollisions = (G4int)G4Poisson(meanNumber);

  //G4cout << "N= " << numOfCollisions << G4endl;

  if( 0 == numOfCollisions) return 0.0; 

  for(G4int i=0; i< numOfCollisions; ++i) 
  {
    G4double rand = G4UniformRand();
    position = dNdxCut1 + ((*v1)[0]/e1 - dNdxCut1)*rand;
    omega = GetEnergyTransfer(coupleIndex, iPlace, position);

    //G4cout << "omega(keV)= " << omega/keV << G4endl;

    if(!one) 
    {
      position = dNdxCut2 + ((*v2)[0]/e2 - dNdxCut2)*rand;
      G4double omega2 = GetEnergyTransfer(coupleIndex, iPlace+1, position);
      omega = omega*W1 + omega2*W2;
    }
    //G4cout << "omega(keV)= " << omega/keV << G4endl;

    loss += omega;
    if( loss > kinEnergy) { break; }
  }
  
  // G4cout<<"PAIPhotData AlongStepLoss = "<<loss/keV<<" keV, on step = "
  //<<step/mm<<" mm"<<G4endl; 

  if     ( loss > kinEnergy) loss = kinEnergy; 
  else if( loss < 0.)        loss = 0.;
 
  return loss;
}

////////////////////////////////////////////////////////////////////////

G4double G4PAIPhotData::SampleAlongStepPhotonTransfer(G4int coupleIndex, 
                                                 G4double kinEnergy,
						 G4double scaledTkin,
						 G4double stepFactor) const
{
  G4double loss = 0.0;
  G4double omega; 
  G4double position, E1, E2, W1, W2, W, dNdxCut1, dNdxCut2, meanNumber;

  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;
 
  G4bool one = true;

  if     (scaledTkin >= fParticleEnergyVector->Energy(nPlace)) iPlace = nPlace; 
  else if(scaledTkin > fParticleEnergyVector->Energy(0))          one = false; 

  G4PhysicsLogVector* vcut = fdNdxCutPhotonTable[coupleIndex];
  G4PhysicsVector*      v1 = (*(fPAIphotonBank[coupleIndex]))(iPlace);
  G4PhysicsVector*      v2 = nullptr;

  dNdxCut1    = (*vcut)[iPlace];
  G4double e1 = v1->Energy(0);
  G4double e2 = e1;

  G4double meanN1 = ((*v1)[0]/e1 - dNdxCut1)*stepFactor;

  meanNumber = meanN1;

  // G4cout<<"iPlace = "<<iPlace<< " meanN1= " << meanN1 
  //	<<"    (*v1)[0]/e1 = "<<(*v1)[0]/e1<< " dNdxCut1= " << dNdxCut1 << G4endl;

  dNdxCut2 = dNdxCut1;
  W1 = 1.0;
  W2 = 0.0;
  if(!one) 
  {
    v2 = (*(fPAIphotonBank[coupleIndex]))(iPlace+1);
    dNdxCut2 = (*vcut)[iPlace+1];
    e2 = v2->Energy(0);

    G4double meanN2 = ((*v2)[0]/e2 - dNdxCut2)*stepFactor;

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;
    meanNumber = W1*meanN1 + W2*meanN2;

    //G4cout<<"meanN= " <<  meanNumber << " meanN2= " << meanN2 
    //	  << " dNdxCut2= " << dNdxCut2 << G4endl;
  }
  if( meanNumber <= 0.0) return 0.0; 

  G4int numOfCollisions = (G4int)G4Poisson(meanNumber);

  //G4cout << "N= " << numOfCollisions << G4endl;

  if( 0 == numOfCollisions) return 0.0; 

  for(G4int i=0; i< numOfCollisions; ++i) 
  {
    G4double rand = G4UniformRand();
    position = dNdxCut1 + ((*v1)[0]/e1 - dNdxCut1)*rand;
    omega = GetEnergyPhotonTransfer(coupleIndex, iPlace, position);

    //G4cout << "omega(keV)= " << omega/keV << G4endl;

    if(!one) 
    {
      position = dNdxCut2 + ((*v2)[0]/e2 - dNdxCut2)*rand;
      G4double omega2 = GetEnergyPhotonTransfer(coupleIndex, iPlace+1, position);
      omega = omega*W1 + omega2*W2;
    }
    //G4cout << "omega(keV)= " << omega/keV << G4endl;

    loss += omega;
    if( loss > kinEnergy) { break; }
  }
  
  // G4cout<<"PAIPhotData AlongStepLoss = "<<loss/keV<<" keV, on step = "
  //<<step/mm<<" mm"<<G4endl; 

  if     ( loss > kinEnergy) loss = kinEnergy; 
  else if( loss < 0.)        loss = 0.;
 
  return loss;
}

//////////////////////////////////////////////////////////////////

G4double G4PAIPhotData::SampleAlongStepPlasmonTransfer(G4int coupleIndex, 
                                                 G4double kinEnergy,
						 G4double scaledTkin,
						 G4double stepFactor) const
{
  G4double loss = 0.0;
  G4double omega; 
  G4double position, E1, E2, W1, W2, W, dNdxCut1, dNdxCut2, meanNumber;

  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;
 
  G4bool one = true;

  if     (scaledTkin >= fParticleEnergyVector->Energy(nPlace)) iPlace = nPlace; 
  else if(scaledTkin > fParticleEnergyVector->Energy(0))          one = false; 

  G4PhysicsLogVector* vcut = fdNdxCutPlasmonTable[coupleIndex];
  G4PhysicsVector*      v1 = (*(fPAIplasmonBank[coupleIndex]))(iPlace);
  G4PhysicsVector*      v2 = nullptr;

  dNdxCut1    = (*vcut)[iPlace];
  G4double e1 = v1->Energy(0);
  G4double e2 = e1;

  G4double meanN1 = ((*v1)[0]/e1 - dNdxCut1)*stepFactor;

  meanNumber = meanN1;

  // G4cout<<"iPlace = "<<iPlace<< " meanN1= " << meanN1 
  //	<<"    (*v1)[0]/e1 = "<<(*v1)[0]/e1<< " dNdxCut1= " << dNdxCut1 << G4endl;

  dNdxCut2 = dNdxCut1;
  W1 = 1.0;
  W2 = 0.0;
  if(!one) 
  {
    v2 = (*(fPAIplasmonBank[coupleIndex]))(iPlace+1);
    dNdxCut2 = (*vcut)[iPlace+1];
    e2 = v2->Energy(0);

    G4double meanN2 = ((*v2)[0]/e2 - dNdxCut2)*stepFactor;

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;
    meanNumber = W1*meanN1 + W2*meanN2;

    //G4cout<<"meanN= " <<  meanNumber << " meanN2= " << meanN2 
    //	  << " dNdxCut2= " << dNdxCut2 << G4endl;
  }
  if( meanNumber <= 0.0) return 0.0; 

  G4int numOfCollisions = (G4int)G4Poisson(meanNumber);

  //G4cout << "N= " << numOfCollisions << G4endl;

  if( 0 == numOfCollisions) return 0.0; 

  for(G4int i=0; i< numOfCollisions; ++i) 
  {
    G4double rand = G4UniformRand();
    position = dNdxCut1 + ((*v1)[0]/e1 - dNdxCut1)*rand;
    omega = GetEnergyPlasmonTransfer(coupleIndex, iPlace, position);

    //G4cout << "omega(keV)= " << omega/keV << G4endl;

    if(!one) 
    {
      position = dNdxCut2 + ((*v2)[0]/e2 - dNdxCut2)*rand;
      G4double omega2 = GetEnergyPlasmonTransfer(coupleIndex, iPlace+1, position);
      omega = omega*W1 + omega2*W2;
    }
    //G4cout << "omega(keV)= " << omega/keV << G4endl;

    loss += omega;
    if( loss > kinEnergy) { break; }
  }
  
  // G4cout<<"PAIPhotData AlongStepLoss = "<<loss/keV<<" keV, on step = "
  //<<step/mm<<" mm"<<G4endl; 

  if     ( loss > kinEnergy) loss = kinEnergy; 
  else if( loss < 0.)        loss = 0.;
 
  return loss;
}

///////////////////////////////////////////////////////////////////////
//
// Returns post step PAI energy transfer > cut electron energy 
// according to passed scaled kinetic energy of particle

G4double G4PAIPhotData::SamplePostStepTransfer(G4int coupleIndex, 
						G4double scaledTkin) const
{  
  //G4cout<<"G4PAIPhotData::GetPostStepTransfer"<<G4endl;
  G4double transfer = 0.0;
  G4double rand = G4UniformRand();

  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

  //  std::size_t iTransfer, iTr1, iTr2;
  G4double position, dNdxCut1, dNdxCut2, E1, E2, W1, W2, W;

  G4PhysicsVector* cutv = fdNdxCutTable[coupleIndex];

  // Fermi plato, try from left
  if( scaledTkin >= fParticleEnergyVector->GetMaxEnergy()) 
  {
    position = (*cutv)[nPlace]*rand;
    transfer = GetEnergyTransfer(coupleIndex, nPlace, position);
  }
  else if( scaledTkin <= fParticleEnergyVector->Energy(0) )
  {
    position = (*cutv)[0]*rand;
    transfer = GetEnergyTransfer(coupleIndex, 0, position);
  }
  else 
  {  
    std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);

    dNdxCut1 = (*cutv)[iPlace];  
    dNdxCut2 = (*cutv)[iPlace+1];  

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W  = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;

    //G4cout<<"iPlace= " << "  dNdxCut1 = "<<dNdxCut1 
    //	  <<" dNdxCut2 = "<<dNdxCut2<< " W1= " << W1 << " W2= " << W2 <<G4endl;

    position = dNdxCut1*rand;
    G4double tr1 = GetEnergyTransfer(coupleIndex, iPlace, position);

    position = dNdxCut2*rand;
    G4double tr2 = GetEnergyTransfer(coupleIndex, iPlace+1, position);

    transfer = tr1*W1 + tr2*W2;
  }
  //G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"<<G4endl; 
  if(transfer < 0.0 ) { transfer = 0.0; }
  return transfer;
}

////////////////////////////////////////////////////////////////////////

G4double G4PAIPhotData::SamplePostStepPhotonTransfer(G4int coupleIndex, 
						G4double scaledTkin) const
{  
  //G4cout<<"G4PAIPhotData::GetPostStepTransfer"<<G4endl;
  G4double transfer = 0.0;
  G4double rand = G4UniformRand();

  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

  //  std::size_t iTransfer, iTr1, iTr2;
  G4double position, dNdxCut1, dNdxCut2, E1, E2, W1, W2, W;

  G4PhysicsVector* cutv = fdNdxCutPhotonTable[coupleIndex];

  // Fermi plato, try from left

  if( scaledTkin >= fParticleEnergyVector->GetMaxEnergy()) 
  {
    position = (*cutv)[nPlace]*rand;
    transfer = GetEnergyPhotonTransfer(coupleIndex, nPlace, position);
  }
  else if( scaledTkin <= fParticleEnergyVector->Energy(0) )
  {
    position = (*cutv)[0]*rand;
    transfer = GetEnergyPhotonTransfer(coupleIndex, 0, position);
  }
  else 
  {  
    std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);

    dNdxCut1 = (*cutv)[iPlace];  
    dNdxCut2 = (*cutv)[iPlace+1];  

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W  = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;

    //G4cout<<"iPlace= " << "  dNdxCut1 = "<<dNdxCut1 
    //	  <<" dNdxCut2 = "<<dNdxCut2<< " W1= " << W1 << " W2= " << W2 <<G4endl;

    position = dNdxCut1*rand;

    G4double tr1 = GetEnergyPhotonTransfer(coupleIndex, iPlace, position);

    position = dNdxCut2*rand;
    G4double tr2 = GetEnergyPhotonTransfer(coupleIndex, iPlace+1, position);

    transfer = tr1*W1 + tr2*W2;
  }
  //G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"<<G4endl; 
  if(transfer < 0.0 ) { transfer = 0.0; }
  return transfer;
}

//////////////////////////////////////////////////////////////////////////

G4double G4PAIPhotData::SamplePostStepPlasmonTransfer(G4int coupleIndex, 
						G4double scaledTkin) const
{  
  //G4cout<<"G4PAIPhotData::GetPostStepTransfer"<<G4endl;
  G4double transfer = 0.0;
  G4double rand = G4UniformRand();

  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

  //  std::size_t iTransfer, iTr1, iTr2;
  G4double position, dNdxCut1, dNdxCut2, E1, E2, W1, W2, W;

  G4PhysicsVector* cutv = fdNdxCutPlasmonTable[coupleIndex];

  // Fermi plato, try from left
  if( scaledTkin >= fParticleEnergyVector->GetMaxEnergy()) 
  {
    position = (*cutv)[nPlace]*rand;
    transfer = GetEnergyPlasmonTransfer(coupleIndex, nPlace, position);
  }
  else if( scaledTkin <= fParticleEnergyVector->Energy(0) )
  {
    position = (*cutv)[0]*rand;
    transfer = GetEnergyPlasmonTransfer(coupleIndex, 0, position);
  }
  else 
  {  
    std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);

    dNdxCut1 = (*cutv)[iPlace];  
    dNdxCut2 = (*cutv)[iPlace+1];  

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W  = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;

    //G4cout<<"iPlace= " << "  dNdxCut1 = "<<dNdxCut1 
    //	  <<" dNdxCut2 = "<<dNdxCut2<< " W1= " << W1 << " W2= " << W2 <<G4endl;

    position = dNdxCut1*rand;
    G4double tr1 = GetEnergyPlasmonTransfer(coupleIndex, iPlace, position);

    position = dNdxCut2*rand;
    G4double tr2 = GetEnergyPlasmonTransfer(coupleIndex, iPlace+1, position);

    transfer = tr1*W1 + tr2*W2;
  }
  //G4cout<<"PAImodel PostStepPlasmonTransfer = "<<transfer/keV<<" keV"<<G4endl; 

  if(transfer < 0.0 )  transfer = 0.0;
 
  return transfer;
}

///////////////////////////////////////////////////////////////////////
//
// Returns PAI energy transfer according to passed 
// indexes of particle kinetic enegry and random x-section

G4double G4PAIPhotData::GetEnergyTransfer(G4int coupleIndex, 
					 std::size_t iPlace, 
					 G4double position) const
{ 
  G4PhysicsVector* v = (*(fPAIxscBank[coupleIndex]))(iPlace); 
  if(position*v->Energy(0) >= (*v)[0]) { return v->Energy(0); }

  std::size_t iTransferMax = v->GetVectorLength() - 1;

  std::size_t iTransfer;
  G4double x1(0.0), x2(0.0), y1(0.0), y2(0.0), energyTransfer;

  for(iTransfer=1; iTransfer<=iTransferMax; ++iTransfer) {
    x2 = v->Energy(iTransfer);
    y2 = (*v)[iTransfer]/x2;
    if(position >= y2) { break; }
  }

  x1 = v->Energy(iTransfer-1);
  y1 = (*v)[iTransfer-1]/x1;
  //G4cout << "i= " << iTransfer << " imax= " << iTransferMax
  //	 << " x1= " << x1 << " x2= " << x2 << G4endl;

  energyTransfer = x1;
  if ( x1 != x2 ) {
    if ( y1 == y2  ) {
      energyTransfer += (x2 - x1)*G4UniformRand();
    } else {
      if(x1*1.1 < x2) {
	const G4int nbins = 5;
        G4double del = (x2 - x1)/G4int(nbins);
        x2  = x1;
        for(G4int i=1; i<=nbins; ++i) {
          x2 += del;
          y2 = v->Value(x2)/x2;
          if(position >= y2) { break; }
          x1 = x2;
          y1 = y2;
	}
      }
      energyTransfer = (y2 - y1)*x1*x2/(position*(x1 - x2) - y1*x1 + y2*x2);
    }
  }
  //  G4cout << "x1(keV)= " << x1/keV << " x2(keV)= " << x2/keV
  //	 << " y1= " << y1 << " y2= " << y2 << " pos= " << position
  //	 << " E(keV)= " << energyTransfer/keV << G4endl; 
  return energyTransfer;
}

/////////////////////////////////////////////////////////////////

G4double G4PAIPhotData::GetEnergyPhotonTransfer(G4int coupleIndex, 
					        std::size_t iPlace, 
					        G4double position) const
{ 
  G4PhysicsVector* v = (*(fPAIphotonBank[coupleIndex]))(iPlace); 
  if(position*v->Energy(0) >= (*v)[0])  return v->Energy(0); 

  std::size_t iTransferMax = v->GetVectorLength() - 1;

  std::size_t iTransfer;
  G4double x1(0.0), x2(0.0), y1(0.0), y2(0.0), energyTransfer;

  for(iTransfer=1; iTransfer<=iTransferMax; ++iTransfer) 
  {
    x2 = v->Energy(iTransfer);
    y2 = (*v)[iTransfer]/x2;
    if(position >= y2)  break; 
  }
  x1 = v->Energy(iTransfer-1);
  y1 = (*v)[iTransfer-1]/x1;

  //G4cout << "i= " << iTransfer << " imax= " << iTransferMax
  //	 << " x1= " << x1 << " x2= " << x2 << G4endl;

  energyTransfer = x1;

  if ( x1 != x2 ) 
  {
    if ( y1 == y2  ) 
    {
      energyTransfer += (x2 - x1)*G4UniformRand();
    } 
    else 
    {
      if( x1*1.1 < x2 ) 
      {
	const G4int nbins = 5;
        G4double del = (x2 - x1)/G4int(nbins);
        x2  = x1;

        for(G4int i=1; i<=nbins; ++i) 
        {
          x2 += del;
          y2 = v->Value(x2)/x2;
          if(position >= y2) { break; }
          x1 = x2;
          y1 = y2;
	}
      }
      energyTransfer = (y2 - y1)*x1*x2/(position*(x1 - x2) - y1*x1 + y2*x2);
    }
  }
  //  G4cout << "x1(keV)= " << x1/keV << " x2(keV)= " << x2/keV
  //	 << " y1= " << y1 << " y2= " << y2 << " pos= " << position
  //	 << " E(keV)= " << energyTransfer/keV << G4endl;
 
  return energyTransfer;
}

/////////////////////////////////////////////////////////////////////////

G4double G4PAIPhotData::GetEnergyPlasmonTransfer(G4int coupleIndex, 
					         std::size_t iPlace, 
					         G4double position) const
{ 
  G4PhysicsVector* v = (*(fPAIplasmonBank[coupleIndex]))(iPlace); 

  if( position*v->Energy(0) >= (*v)[0] )  return v->Energy(0); 

  std::size_t iTransferMax = v->GetVectorLength() - 1;

  std::size_t iTransfer;
  G4double x1(0.0), x2(0.0), y1(0.0), y2(0.0), energyTransfer;

  for(iTransfer = 1; iTransfer <= iTransferMax; ++iTransfer) 
  {
    x2 = v->Energy(iTransfer);
    y2 = (*v)[iTransfer]/x2;
    if(position >= y2)  break; 
  }
  x1 = v->Energy(iTransfer-1);
  y1 = (*v)[iTransfer-1]/x1;

  //G4cout << "i= " << iTransfer << " imax= " << iTransferMax
  //	 << " x1= " << x1 << " x2= " << x2 << G4endl;

  energyTransfer = x1;

  if ( x1 != x2 ) 
  {
    if ( y1 == y2  ) 
    {
      energyTransfer += (x2 - x1)*G4UniformRand();
    } 
    else 
    {
      if(x1*1.1 < x2) 
      {
	const G4int nbins = 5;
        G4double del = (x2 - x1)/G4int(nbins);
        x2  = x1;

        for( G4int i = 1; i <= nbins; ++i ) 
        {
          x2 += del;
          y2 = v->Value(x2)/x2;

          if(position >= y2)  break; 

          x1 = x2;
          y1 = y2;
	}
      }
      energyTransfer = (y2 - y1)*x1*x2/(position*(x1 - x2) - y1*x1 + y2*x2);
    }
  }
  //  G4cout << "x1(keV)= " << x1/keV << " x2(keV)= " << x2/keV
  //	 << " y1= " << y1 << " y2= " << y2 << " pos= " << position
  //	 << " E(keV)= " << energyTransfer/keV << G4endl; 

  return energyTransfer;
}

//
//
//////////////////////////////////////////////////////////////////////


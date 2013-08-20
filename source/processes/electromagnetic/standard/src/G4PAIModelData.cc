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
// $Id: G4PAIModelData.cc 72008 2013-07-03 08:46:39Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class
// File name:     G4PAIModelData.cc
//
// Author:        V. Ivanchenko based on V.Grichine code of G4PAIModel
//
// Creation date: 16.08.2013
//
// Modifications:
//

#include "G4PAIModelData.hh"
#include "G4PAIModel.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4SandiaTable.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"

////////////////////////////////////////////////////////////////////////

using namespace std;

G4PAIModelData::G4PAIModelData(G4double tmin, G4double tmax, G4int ver)
{ 
  const G4int nPerDecade = 20; 
  const G4double lowestTkin = 0.1*MeV;
  const G4double highestTkin = 10*TeV;

  fPAIySection.SetVerbose(ver);

  fLowestKineticEnergy  = std::max(tmin, lowestTkin);
  fHighestKineticEnergy = tmax;
  if(tmax < 10*fLowestKineticEnergy) { 
    fHighestKineticEnergy = 10*fLowestKineticEnergy;
  } else if(tmax > highestTkin) {
    fHighestKineticEnergy = std::max(highestTkin, 10*fLowestKineticEnergy);
  }
  fTotBin = (G4int)(nPerDecade*
		    std::log10(fHighestKineticEnergy/fLowestKineticEnergy));

  fParticleEnergyVector = new G4PhysicsLogVector(fLowestKineticEnergy,
						 fHighestKineticEnergy,
						 fTotBin);
  if(0 < ver) {
    G4cout << "G4PAIModelData: Nbins= " << fTotBin
	   << " Tmin(MeV)= " << fLowestKineticEnergy/MeV
	   << " Tmax(GeV)= " << fHighestKineticEnergy/GeV << G4endl;
  }
}

////////////////////////////////////////////////////////////////////////////

G4PAIModelData::~G4PAIModelData()
{
  delete fParticleEnergyVector;

  size_t n = fPAIxscBank.size();
  if(0 < n) {
    for(size_t i=0; i<n; ++i) {
      if(fPAIxscBank[i]) {
	  fPAIxscBank[i]->clearAndDestroy();
	  delete fPAIxscBank[i];
      }
      if(fPAIdEdxBank[i]) {
	fPAIdEdxBank[i]->clearAndDestroy();
	delete fPAIdEdxBank[i];
      }
      delete fdEdxTable[i];
      delete fdNdxCutTable[i];
    }
  }
}

///////////////////////////////////////////////////////////////////////////////

void G4PAIModelData::Initialise(const G4MaterialCutsCouple* couple,
                                G4double cut, G4PAIModel* model)
{
  const G4Material* mat = couple->GetMaterial();     
  fSandia.Initialize(const_cast<G4Material*>(mat));

  G4PhysicsTable* PAItransferTable = new G4PhysicsTable(fTotBin+1);
  G4PhysicsTable* PAIdEdxTable = new G4PhysicsTable(fTotBin+1);

  G4PhysicsLogVector* dEdxCutVector =
    new G4PhysicsLogVector(fLowestKineticEnergy,
			   fHighestKineticEnergy,
			   fTotBin);

  G4PhysicsLogVector* dNdxCutVector = 
    new G4PhysicsLogVector(fLowestKineticEnergy,
			   fHighestKineticEnergy,
			   fTotBin);

  // low energy Sandia interval
  G4double Tmin = fSandia.GetSandiaMatTablePAI(0,0); 

  // energy safety
  const G4double deltaLow = 100.*eV; 

  for (G4int i = 0; i <= fTotBin; ++i) {

    G4double kinEnergy = fParticleEnergyVector->Energy(i);
    G4double Tmax = model->ComputeMaxEnergy(kinEnergy);
    G4double tau = kinEnergy/proton_mass_c2;
    G4double bg2 = tau*( tau + 2. );

    G4double Tkin = Tmax;
    if (Tmax < Tmin + deltaLow ) { Tkin = Tmin + deltaLow; }

    fPAIySection.Initialize(mat, Tkin, bg2, &fSandia);

    // G4cout<<"ionloss = "<<ionloss*cm/keV<<" keV/cm"<<endl;

    G4int n = fPAIySection.GetSplineSize();
    G4PhysicsFreeVector* transferVector = new G4PhysicsFreeVector(n);
    G4PhysicsFreeVector* dEdxVector = new G4PhysicsFreeVector(n);

    for( G4int k = 0; k < n; k++ )
    {
      transferVector->PutValue( k ,
                                fPAIySection.GetSplineEnergy(k+1),
                                fPAIySection.GetIntegralPAIySection(k+1) );
      dEdxVector->PutValue( k ,
                                fPAIySection.GetSplineEnergy(k+1),
                                fPAIySection.GetIntegralPAIdEdx(k+1) );
    }
    //G4cout << *transferVector << G4endl;

    G4double ionloss = fPAIySection.GetMeanEnergyLoss();//  total <dE/dx>

    if(ionloss < 0.0) { ionloss = 0.0; }
    dEdxCutVector->PutValue(i,ionloss);

    G4double dNdxCut = transferVector->Value(cut);
    //G4cout << "i= " << i << " x= " << dNdxCut << G4endl;
    if(dNdxCut < 0.0) { dNdxCut = 0.0; }
    dNdxCutVector->PutValue(i, dNdxCut);

    PAItransferTable->insertAt(i,transferVector);
    PAIdEdxTable->insertAt(i,dEdxVector);

  } // end of Tkin loop

  fPAIxscBank.push_back(PAItransferTable);
  fPAIdEdxBank.push_back(PAIdEdxTable);
  fdEdxTable.push_back(dEdxCutVector);
  fdNdxCutTable.push_back(dNdxCutVector);

}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModelData::DEDXPerVolume(G4int coupleIndex, G4double scaledTkin,
			 G4double cut) const
{
  // VI: iPlace is the low edge index of the bin
  // iPlace is in interval from 0 to (N-1)
  size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);

  // VI: apply interpolation of the vector
  G4double dEdx = fdEdxTable[coupleIndex]->Value(scaledTkin) - 
    (*(fPAIdEdxBank[coupleIndex]))(iPlace)->Value(cut);
  if( dEdx < 0.) { dEdx = 0.; }
  return dEdx;
}

/////////////////////////////////////////////////////////////////////////

G4double 
G4PAIModelData::CrossSectionPerVolume(G4int coupleIndex, 
				      G4double scaledTkin,
				      G4double tcut, G4double tmax) const
{
  G4double cross, cross1, cross2;

  // iPlace is in interval from 0 to (N-1)
  size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);

  //G4cout<<"iPlace = "<<iPlace<<"; tmax = "
  // <<tmax<<"; cutEnergy = "<<cutEnergy<<G4endl;  
  cross1 = (*(fPAIxscBank[coupleIndex]))(iPlace)->Value(tmax);
  //G4cout<<"cross1 = "<<cross1<<G4endl;  
  cross2 = (*(fPAIxscBank[coupleIndex]))(iPlace)->Value(tcut);
  //G4cout<<"cross2 = "<<cross2<<G4endl;  
  cross  = (cross2-cross1);
  //G4cout<<"cross = "<<cross<<G4endl;  
  if( cross < 0.0) { cross = 0.0; }
  return cross;
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIModelData::SampleAlongStepTransfer(G4int coupleIndex, 
                                                 G4double kinEnergy,
						 G4double scaledTkin,
						 G4double stepFactor) const
{
  G4double loss = 0.0;
  G4double omega; 
  G4double position, E1, E2, W1, W2, W, dNdxCut1, dNdxCut2, meanNumber;

  size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;
 
  G4bool one = true;
  if(scaledTkin >= fParticleEnergyVector->Energy(nPlace)) { iPlace = nPlace; }
  else if(scaledTkin > fParticleEnergyVector->Energy(0)) { 
    one = false; 
  }
  dNdxCut1 = (*(fdNdxCutTable[coupleIndex]))(iPlace);

  G4double meanN1 =
    ((*(*(fPAIxscBank[coupleIndex]))(iPlace))(0) - dNdxCut1)*stepFactor;
  meanNumber = meanN1;

  //G4cout<<"iPlace = "<<iPlace<< " meanN1= " << meanN1 
  //	<< " dNdxCut1= " << dNdxCut1 << G4endl;

  dNdxCut2 = dNdxCut1;
  W1 = 1.0;
  W2 = 0.0;
  if(!one) {
    dNdxCut2 = (*(fdNdxCutTable[coupleIndex]))(iPlace+1);
    G4double meanN2 =
      ((*(*(fPAIxscBank[coupleIndex]))(iPlace+1))(0) - dNdxCut2)*stepFactor;
    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;
    meanNumber = W1*meanN1 + W2*meanN2;
    //G4cout<<"meanN= " <<  meanNumber << " meanN2= " << meanN2 
    //	  << " dNdxCut2= " << dNdxCut2 << G4endl;
  }
  if(meanNumber < 0.0) { return 0.0; }

  G4int numOfCollisions = G4Poisson(meanNumber);

  //G4cout << "N= " << numOfCollisions << G4endl;

  if(0 == numOfCollisions) { return 0.0; }

  for(G4int i=0; i< numOfCollisions; ++i) {
    G4double rand = G4UniformRand();
    position = dNdxCut1+
      ((*(*(fPAIxscBank[coupleIndex]))(iPlace))(0) - dNdxCut1)*rand;
    omega = GetEnergyTransfer(coupleIndex, iPlace, position);
    //G4cout << "omega(keV)= " << omega/keV << G4endl;
    if(!one) {
      position = dNdxCut2+
	((*(*(fPAIxscBank[coupleIndex]))(iPlace+1))(0) - dNdxCut2)*rand;
      G4double omega2 = GetEnergyTransfer(coupleIndex,iPlace+1, position);
      omega = omega*W1 + omega2*W2;
    }
    //G4cout << "omega(keV)= " << omega/keV << G4endl;

    loss += omega;
    if(loss > kinEnergy) { break; }
  }
  
  // G4cout<<"PAIModelData AlongStepLoss = "<<loss/keV<<" keV, on step = "
  //<<step/mm<<" mm"<<G4endl; 
  if(loss > kinEnergy) { loss = kinEnergy; }
  else if(loss < 0.)   { loss = 0.; }
  return loss;
}

///////////////////////////////////////////////////////////////////////
//
// Returns post step PAI energy transfer > cut electron energy 
// according to passed scaled kinetic energy of particle

G4double G4PAIModelData::SamplePostStepTransfer(G4int coupleIndex, 
						G4double scaledTkin) const
{  
  // G4cout<<"G4PAIModelData::GetPostStepTransfer"<<G4endl;
  G4double transfer = 0.0;
  G4double rand = G4UniformRand();

  size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

  //  size_t iTransfer, iTr1, iTr2;
  G4double position, dNdxCut1, dNdxCut2, E1, E2, W1, W2, W;

  // Fermi plato, try from left
  if(scaledTkin >= fParticleEnergyVector->GetMaxEnergy()) 
  {
    position = (*(fdNdxCutTable[coupleIndex]))(nPlace)*rand;
    transfer = GetEnergyTransfer(coupleIndex, iPlace, position);
  }
  else if(scaledTkin <= fParticleEnergyVector->Energy(0))
  {
    position = (*(fdNdxCutTable[coupleIndex]))(0)*rand;
    transfer = GetEnergyTransfer(coupleIndex, iPlace, position);
  }
  else 
  {  
    dNdxCut1 = (*(fdNdxCutTable[coupleIndex]))(iPlace);  
    dNdxCut2 = (*(fdNdxCutTable[coupleIndex]))(iPlace+1);  

    // G4cout<<"dNdxCut1 = "<<dNdxCut1<<G4endl;
    // G4cout<<"dNdxCut2 = "<<dNdxCut2<<G4endl;

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W  = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;

    position = dNdxCut1*rand;
    G4double tr1 = GetEnergyTransfer(coupleIndex, iPlace, position);

    position = dNdxCut2*rand;
    G4double tr2 = GetEnergyTransfer(coupleIndex, iPlace+1, position);

    transfer = tr1*W1 + tr2*W2;
  }
  // G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"<<G4endl; 
  if(transfer < 0.0 ) { transfer = 0.0; }
  return transfer;
}

///////////////////////////////////////////////////////////////////////
//
// Returns PAI energy transfer according to passed 
// indexes of particle kinetic enegry and random x-section

G4double G4PAIModelData::GetEnergyTransfer(G4int coupleIndex, 
					   size_t iPlace, 
					   G4double position) const
{ 
  if(position > (*(*(fPAIxscBank[coupleIndex]))(iPlace))[0]) {
    return (*(*(fPAIxscBank[coupleIndex]))(iPlace))[0];
  }

  size_t iTransfer;
  size_t iTransferMax = 
    (*(fPAIxscBank[coupleIndex]))(iPlace)->GetVectorLength();
  G4double x1, x2, y1, y2(0.0), energyTransfer;

  for(iTransfer=1; iTransfer < iTransferMax; ++iTransfer) {
    y2 = (*(*(fPAIxscBank[coupleIndex]))(iPlace))[iTransfer];
    if(position > y2) { break; }
  }

  y1 = (*(*(fPAIxscBank[coupleIndex]))(iPlace))[iTransfer-1];

  x1 = (*(fPAIxscBank[coupleIndex]))(iPlace)->Energy(iTransfer-1);
  x2 = (*(fPAIxscBank[coupleIndex]))(iPlace)->Energy(iTransfer);

  energyTransfer = x1;
  if ( x1 != x2 ) {
    if ( y1 == y2  ) {
      energyTransfer += (x2 - x1)*G4UniformRand();
    } else {
      energyTransfer += (position - y1)*(x2 - x1)/(y2 - y1);
    }
  }
  //G4cout << "x1(keV)= " << x1/keV << " x2(keV)= " << x2/keV
  //	 << " y1= " << y1 << " y2= " << y2 << " pos= " << position
  //	 << " E(keV)= " << energyTransfer/keV << G4endl; 
  return energyTransfer;
}

//////////////////////////////////////////////////////////////////////


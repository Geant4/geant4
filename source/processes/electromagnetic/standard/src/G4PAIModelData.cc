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
  const G4int nPerDecade = 10; 
  const G4double lowestTkin = 50*keV;
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
    G4cout << "### G4PAIModelData: Nbins= " << fTotBin
	   << " Tlowest(MeV)= " << fLowestKineticEnergy/MeV
	   << " Tmin(keV)= " << tmin/keV 
	   << " Tmax(GeV)= " << fHighestKineticEnergy/GeV 
	   << G4endl;
  }
}

////////////////////////////////////////////////////////////////////////////

G4PAIModelData::~G4PAIModelData()
{
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
      delete fdEdxCutTable[i];
    }
  }
  delete fParticleEnergyVector;
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

  G4PhysicsLogVector* dEdxMeanVector =
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

    if (Tmax < Tmin + deltaLow ) { Tmax = Tmin + deltaLow; }

    fPAIySection.Initialize(mat, Tmax, bg2, &fSandia);
    /*
    G4cout << i << ". TransferMax(keV)= "<< Tmax/keV << "  cut(keV)= " 
    	   << cut/keV << "  E(MeV)= " << kinEnergy/MeV << G4endl;
    */
    G4int n = fPAIySection.GetSplineSize();
    G4int kmin = 0;
    for(G4int k = 0; k < n; ++k) {
      if(fPAIySection.GetIntegralPAIySection(k+1) <= 0.0) { 
	kmin = k;
      } else {
	break;
      }
    }
    n -= kmin;

    G4PhysicsFreeVector* transferVector = new G4PhysicsFreeVector(n);
    G4PhysicsFreeVector* dEdxVector = new G4PhysicsFreeVector(n);

    //G4double tr0 = 0.0;
    G4double tr = 0.0;
    for(G4int k = kmin; k < n; ++k)
    {
      G4double t  = fPAIySection.GetSplineEnergy(k+1);
      tr = fPAIySection.GetIntegralPAIySection(k+1);
      //if(tr >= tr0) { tr0 = tr; }
      //else { G4cout << "G4PAIModelData::Initialise Warning: Ekin(MeV)= "
      //		    << t/MeV << " IntegralTransfer= " << tr 
      //		    << " < " << tr0 << G4endl; }
      transferVector->PutValue(k , t, t*tr);
      dEdxVector->PutValue(k, t, fPAIySection.GetIntegralPAIdEdx(k+1));
    }
    //G4cout << *transferVector << G4endl;

    G4double ionloss = fPAIySection.GetMeanEnergyLoss();//  total <dE/dx>

    if(ionloss < 0.0) ionloss = 0.0; 

    dEdxMeanVector->PutValue(i,ionloss);

    G4double dNdxCut = transferVector->Value(cut)/cut;
    G4double dEdxCut = dEdxVector->Value(cut);
    //G4cout << "i= " << i << " x= " << dNdxCut << G4endl;
    if(dNdxCut < 0.0) { dNdxCut = 0.0; }
    if(dEdxCut < 0.0) { dEdxCut = 0.0; }
    dNdxCutVector->PutValue(i, dNdxCut);
    dEdxCutVector->PutValue(i, dEdxCut);

    PAItransferTable->insertAt(i,transferVector);
    PAIdEdxTable->insertAt(i,dEdxVector);

    //transferVector->SetSpline(true);
    //transferVector->FillSecondDerivatives();
    //dEdxVector->SetSpline(true);
    //dEdxVector->FillSecondDerivatives();

  } // end of Tkin loop`
  fPAIxscBank.push_back(PAItransferTable);
  fPAIdEdxBank.push_back(PAIdEdxTable);
  /*
  dEdxMeanVector->SetSpline(true);
  dEdxMeanVector->FillSecondDerivatives();
  dNdxCutVector->SetSpline(true);
  dNdxCutVector->FillSecondDerivatives();
  dEdxCutVector->SetSpline(true);
  dEdxCutVector->FillSecondDerivatives();
  */
  fdEdxTable.push_back(dEdxMeanVector);
  fdNdxCutTable.push_back(dNdxCutVector);
  fdEdxCutTable.push_back(dEdxCutVector);
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModelData::DEDXPerVolume(G4int coupleIndex, G4double scaledTkin,
			 G4double cut) const
{
  // VI: iPlace is the low edge index of the bin
  // iPlace is in interval from 0 to (N-1)
  size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

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
G4PAIModelData::CrossSectionPerVolume(G4int coupleIndex, 
				      G4double scaledTkin,
				      G4double tcut, G4double tmax) const
{
  G4double cross, cross1, cross2;

  // iPlace is in interval from 0 to (N-1)
  size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

  G4bool one = true;
  if(scaledTkin >= fParticleEnergyVector->Energy(nPlace)) { iPlace = nPlace; }
  else if(scaledTkin > fParticleEnergyVector->Energy(0)) { 
    one = false; 
  }
  G4PhysicsTable* table = fPAIxscBank[coupleIndex];

  //G4cout<<"iPlace = "<<iPlace<<"; tmax = "
  // <<tmax<<"; cutEnergy = "<<cutEnergy<<G4endl;  
  cross1 = (*table)(iPlace)->Value(tmax)/tmax;
  //G4cout<<"cross1 = "<<cross1<<G4endl;  
  cross2 = (*table)(iPlace)->Value(tcut)/tcut;
  //G4cout<<"cross2 = "<<cross2<<G4endl;  
  cross  = (cross2-cross1);
  //G4cout<<"cross = "<<cross<<G4endl;  
  if(!one) {
    cross2 = (*table)(iPlace+1)->Value(tcut)/tcut 
      - (*table)(iPlace+1)->Value(tmax)/tmax;

    G4double E1 = fParticleEnergyVector->Energy(iPlace); 
    G4double E2 = fParticleEnergyVector->Energy(iPlace+1);
    G4double W  = 1.0/(E2 - E1);
    G4double W1 = (E2 - scaledTkin)*W;
    G4double W2 = (scaledTkin - E1)*W;
    cross *= W1;
    cross += W2*cross2;
  }

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
  G4PhysicsLogVector* vcut = fdNdxCutTable[coupleIndex];
  G4PhysicsVector* v1 = (*(fPAIxscBank[coupleIndex]))(iPlace);
  G4PhysicsVector* v2 = 0;

  dNdxCut1 = (*vcut)[iPlace];
  G4double e1 = v1->Energy(0);
  G4double e2 = e1;

  G4double meanN1 = ((*v1)[0]/e1 - dNdxCut1)*stepFactor;
  meanNumber = meanN1;

  //G4cout<<"iPlace = "<<iPlace<< " meanN1= " << meanN1 
  //	<< " dNdxCut1= " << dNdxCut1 << G4endl;

  dNdxCut2 = dNdxCut1;
  W1 = 1.0;
  W2 = 0.0;
  if(!one) {
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
  if(meanNumber < 0.0) { return 0.0; }

  G4int numOfCollisions = G4Poisson(meanNumber);

  //G4cout << "N= " << numOfCollisions << G4endl;

  if(0 == numOfCollisions) { return 0.0; }

  for(G4int i=0; i< numOfCollisions; ++i) {
    G4double rand = G4UniformRand();
    position = dNdxCut1 + ((*v1)[0]/e1 - dNdxCut1)*rand;
    omega = GetEnergyTransfer(coupleIndex, iPlace, position);
    //G4cout << "omega(keV)= " << omega/keV << G4endl;
    if(!one) {
      position = dNdxCut2 + ((*v2)[0]/e2 - dNdxCut2)*rand;
      G4double omega2 = GetEnergyTransfer(coupleIndex, iPlace+1, position);
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
						G4double scaledTkin,
						G4double tmax) const
{  
  //G4cout<<"G4PAIModelData::GetPostStepTransfer idx= "<< coupleIndex 
  //	<< " Tkin= " << scaledTkin << "  Tmax= " << tmax << G4endl;
  G4double transfer = 0.0;
  G4double rand = G4UniformRand();

  size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;
  size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);

  G4bool one = true;
  if(scaledTkin >= fParticleEnergyVector->Energy(nPlace)) { iPlace = nPlace; }
  else if(scaledTkin > fParticleEnergyVector->Energy(0)) { 
    one = false; 
  }
  G4PhysicsTable* table = fPAIxscBank[coupleIndex];
  G4PhysicsLogVector* vcut = fdNdxCutTable[coupleIndex];
  G4PhysicsVector* v1 = (*table)[iPlace];

  G4double position;

  //G4cout << *v1 << G4endl;

  G4double dNdxCut1 = (*vcut)[iPlace];  
  G4double emax = std::min(tmax, v1->GetMaxEnergy());
  G4double dNdxCutM = v1->Value(emax)/emax;  
  /*
  G4cout << "iPlace= " << iPlace << " nPlace= " << nPlace << "  emax= " << emax 
	 << " dNdxCut1= " << dNdxCut1 << " dNdxCutM= " << dNdxCutM 
	 << " one= " << one << G4endl;
  */
  if(one) {
    position = dNdxCutM + (dNdxCut1 - dNdxCutM)*rand;
    transfer = GetEnergyTransfer(coupleIndex, iPlace, position);

    //G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"
    //	  << " position= " << position << G4endl; 

  } else {

    G4double E1, E2, W1, W2, W;
    G4double dNdxCut2 = (*vcut)[iPlace+1];  
    G4PhysicsVector* v2 = (*table)[iPlace+1];
    emax = std::min(tmax, v2->GetMaxEnergy());
    G4double dNdxCutM2 = v2->Value(emax)/emax;  

    //G4cout << "  emax2= " << emax 
    //	   << " dNdxCut2= " << dNdxCut2 << " dNdxCutM2= " << dNdxCutM2 << G4endl;

    E1 = fParticleEnergyVector->Energy(iPlace); 
    E2 = fParticleEnergyVector->Energy(iPlace+1);
    W  = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;
    /*
    G4cout<< "E1= " << E1 << " E2= " << E2 <<" iPlace= " << iPlace 
	  << " dNdxCut1 = "<<dNdxCut1 <<" dNdxCut2 = "<<dNdxCut2
	  << " W1= " << W1 << " W2= " << W2 <<G4endl;
    */
    position = dNdxCutM + (dNdxCut1 - dNdxCutM)*rand;
    G4double tr1 = 0.0;
    if(position > 0.0) {
      tr1 = GetEnergyTransfer(coupleIndex, iPlace, position);
    }
    //G4cout<<"PAImodel PostStepTransfer1 = "<<tr1/keV<<" keV"
    //	  << " position= " << position << G4endl; 

    position = dNdxCutM2 + (dNdxCut2 - dNdxCutM2)*rand;
    G4double tr2 = GetEnergyTransfer(coupleIndex, iPlace+1, position);
    //G4cout<<"PAImodel PostStepTransfer2 = "<<tr2/keV<<" keV"
    //	  << " position= " << position << G4endl; 

    transfer = tr1*W1 + tr2*W2;
  }
  //G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"
  //	<< " position= " << position << G4endl; 
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
  G4PhysicsVector* v = (*(fPAIxscBank[coupleIndex]))(iPlace); 
  if(position*v->Energy(0) >= (*v)[0]) { return v->Energy(0); }

  size_t iTransferMax = v->GetVectorLength() - 1;

  size_t iTransfer;
  G4double x1(0.0), x2(0.0), y1(0.0), y2(0.0), energyTransfer;

  //G4cout << "iPlace= " << iPlace << " iTransferMax= " << iTransferMax << G4endl;
  for(iTransfer=1; iTransfer<=iTransferMax; ++iTransfer) {
    x2 = v->Energy(iTransfer);
    y2 = (*v)[iTransfer]/x2;
    if(position >= y2) { break; }
    if(iTransfer == iTransferMax) { return v->GetMaxEnergy(); }
  }

  x1 = v->Energy(iTransfer-1);
  y1 = (*v)[iTransfer-1]/x1;
  /*
  G4cout << "i= " << iTransfer << " imax= " << iTransferMax
  	 << " x1= " << x1 << " x2= " << x2 << G4endl;
  */
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
      //G4cout << "x1(keV)= " << x1/keV << " x2(keV)= " << x2/keV
      //     << " y1= " << y1 << " y2= " << y2 << " pos= " << position << G4endl;
      energyTransfer = (y2 - y1)*x1*x2/(position*(x1 - x2) - y1*x1 + y2*x2);
    }
  }
  //  G4cout << "x1(keV)= " << x1/keV << " x2(keV)= " << x2/keV
  //	 << " y1= " << y1 << " y2= " << y2 << " pos= " << position
  //	 << " E(keV)= " << energyTransfer/keV << G4endl; 
  return energyTransfer;
}

//////////////////////////////////////////////////////////////////////


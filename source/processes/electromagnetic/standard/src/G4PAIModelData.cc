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
	   << " Tlowest(keV)= " << lowestTkin/keV
	   << " Tmin(keV)= " << fLowestKineticEnergy/keV 
	   << " Tmax(GeV)= " << fHighestKineticEnergy/GeV 
	   << G4endl;
  }
}

////////////////////////////////////////////////////////////////////////////

G4PAIModelData::~G4PAIModelData()
{
  std::size_t n = fPAIxscBank.size();
  if(0 < n) {
    for(std::size_t i=0; i<n; ++i) {
      if(fPAIxscBank[i]) {
        fPAIxscBank[i]->clearAndDestroy();
	delete fPAIxscBank[i];
      }
      if(fPAIdEdxBank[i]) {
        fPAIdEdxBank[i]->clearAndDestroy();
	delete fPAIdEdxBank[i];
      }
      delete fdEdxTable[i];
    }
  }
  delete fParticleEnergyVector;
}

///////////////////////////////////////////////////////////////////////////////

void G4PAIModelData::Initialise(const G4MaterialCutsCouple* couple,
                                G4PAIModel* model)
{
  const G4Material* mat = couple->GetMaterial();     
  fSandia.Initialize(const_cast<G4Material*>(mat));

  auto PAItransferTable = new G4PhysicsTable(fTotBin+1);
  auto PAIdEdxTable = new G4PhysicsTable(fTotBin+1);
  auto dEdxMeanVector =
    new G4PhysicsLogVector(fLowestKineticEnergy,
			   fHighestKineticEnergy,
			   fTotBin);
  // low energy Sandia interval
  G4double Tmin = fSandia.GetSandiaMatTablePAI(0,0); 

  // energy safety
  static const G4double deltaLow = 100.*eV; 

  for (G4int i = 0; i <= fTotBin; ++i) {

    G4double kinEnergy = fParticleEnergyVector->Energy(i);
    G4double Tmax = model->ComputeMaxEnergy(kinEnergy);
    G4double tau = kinEnergy/proton_mass_c2;
    G4double bg2 = tau*( tau + 2. );

    if (Tmax < Tmin + deltaLow ) { Tmax = Tmin + deltaLow; }

    fPAIySection.Initialize(mat, Tmax, bg2, &fSandia);
    
    //G4cout << i << ". TransferMax(keV)= "<< Tmax/keV  
    //	   << "  E(MeV)= " << kinEnergy/MeV << G4endl;
    
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

    auto transferVector = new G4PhysicsFreeVector(n);
    auto dEdxVector = new G4PhysicsFreeVector(n);

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
      transferVector->PutValue(k, t, t*tr);
      dEdxVector->PutValue(k, t, fPAIySection.GetIntegralPAIdEdx(k+1));
    }
    //G4cout << "TransferVector:" << G4endl;
    //G4cout << *transferVector << G4endl;
    //G4cout << "DEDXVector:" << G4endl;
    //G4cout << *dEdxVector << G4endl;

    G4double ionloss = fPAIySection.GetMeanEnergyLoss();//  total <dE/dx>

    if(ionloss < 0.0) ionloss = 0.0; 

    dEdxMeanVector->PutValue(i,ionloss);

    PAItransferTable->insertAt(i,transferVector);
    PAIdEdxTable->insertAt(i,dEdxVector);

  } // end of Tkin loop`
  fPAIxscBank.push_back(PAItransferTable);
  fPAIdEdxBank.push_back(PAIdEdxTable);
  //G4cout << "dEdxMeanVector: " << G4endl;
  //G4cout << *dEdxMeanVector << G4endl;
  fdEdxTable.push_back(dEdxMeanVector);
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModelData::DEDXPerVolume(G4int coupleIndex, G4double scaledTkin,
			 G4double cut) const
{
  // VI: iPlace is the low edge index of the bin
  // iPlace is in interval from 0 to (N-1)
  std::size_t iPlace(0);
  G4double dEdx = fdEdxTable[coupleIndex]->Value(scaledTkin, iPlace);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;
  /*
  G4cout << "G4PAIModelData::DEDXPerVolume: coupleIdx= " << coupleIndex
	 << " Tscaled= " << scaledTkin << " cut= " << cut 
	 << " iPlace= " << iPlace << " nPlace= " << nPlace << G4endl;
  */
  G4bool one = true;
  if(scaledTkin >= fParticleEnergyVector->Energy(nPlace)) { iPlace = nPlace; }
  else if(scaledTkin > fParticleEnergyVector->Energy(0)) { 
    one = false; 
  }

  // VI: apply interpolation of the vector
  G4double del  = (*(fPAIdEdxBank[coupleIndex]))(iPlace)->Value(cut);
  //G4cout << "dEdx= " << dEdx << " del= " << del << G4endl;
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
  //G4cout << "dEdx= " << dEdx << " del= " << del << G4endl;

  dEdx = std::max(dEdx, 0.); 
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
  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;

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

  cross = std::max(cross, 0.0); 
  return cross;
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIModelData::SampleAlongStepTransfer(G4int coupleIndex, 
                                                 G4double kinEnergy,
						 G4double scaledTkin,
						 G4double tmax,
						 G4double stepFactor) const
{
  //G4cout << "=== G4PAIModelData::SampleAlongStepTransfer" << G4endl;
  G4double loss = 0.0;

  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);
  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;
 
  G4bool one = true;
  if(scaledTkin >= fParticleEnergyVector->Energy(nPlace)) { iPlace = nPlace; }
  else if(scaledTkin > fParticleEnergyVector->Energy(0)) { 
    one = false; 
  }

  G4double meanNumber = 0.0;
  G4double meanN11 = 0.0;
  G4double meanN12 = 0.0;
  G4double meanN21 = 0.0;
  G4double meanN22 = 0.0;

  G4PhysicsVector* v1 = (*(fPAIxscBank[coupleIndex]))(iPlace);
  G4PhysicsVector* v2 = nullptr;

  G4double e1 = v1->Energy(0);
  G4double e2 = std::min(tmax, v1->GetMaxEnergy());

  if(e2 >= e1) {
    meanN11 = (*v1)[0]/e1;
    meanN12 = v1->Value(e2)/e2;
    meanNumber = (meanN11 - meanN12)*stepFactor;
  }
  //G4cout<<"iPlace = "<<iPlace<< " meanN11= " << meanN11
  //	<< " meanN12= " << meanN12 << G4endl;

  G4double W1 = 1.0;
  G4double W2 = 0.0;
  if(!one) {
    v2 = (*(fPAIxscBank[coupleIndex]))(iPlace+1);

    e1 = v2->Energy(0);
    e2 = std::min(tmax, v2->GetMaxEnergy());
    if(e2 >= e1) {
      meanN21 = (*v2)[0]/e1;
      meanN22 = v2->Value(e2)/e2;
      G4double E1 = fParticleEnergyVector->Energy(iPlace); 
      G4double E2 = fParticleEnergyVector->Energy(iPlace+1);
      G4double W = 1.0/(E2 - E1);
      W1 = (E2 - scaledTkin)*W;
      W2 = (scaledTkin - E1)*W;
      meanNumber *= W1;
      meanNumber += (meanN21 - meanN22)*stepFactor*W2;
    }
  }

  if(meanNumber < 0.0) { return loss; }
  G4int numOfCollisions = (G4int)G4Poisson(meanNumber);

  //G4cout << "meanNumber= " <<  meanNumber << " N= " << numOfCollisions << G4endl;

  if(0 == numOfCollisions) { return loss; }

  G4double position, omega, omega2;
  for(G4int i=0; i< numOfCollisions; ++i) {
    G4double rand = G4UniformRand();
    position = meanN12 + (meanN11 - meanN12)*rand;
    omega = GetEnergyTransfer(coupleIndex, iPlace, position);
    //G4cout << "omega(keV)= " << omega/keV << G4endl;
    if(!one) {
      position = meanN22 + (meanN21 - meanN22)*rand;
      omega2 = GetEnergyTransfer(coupleIndex, iPlace+1, position);
      omega *= W1;
      omega += omega2*W2;
    }
    //G4cout << "omega(keV)= " << omega/keV << G4endl;

    loss += omega;
    if(loss > kinEnergy) { break; }
  }
  
  //G4cout<<"PAIModelData AlongStepLoss = "<<loss/keV<<" keV"<<G4endl; 
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
						G4double tmin,
						G4double tmax) const
{  
  //G4cout<<"=== G4PAIModelData::SamplePostStepTransfer idx= "<< coupleIndex 
  //	<< " Tkin= " << scaledTkin << "  Tmax= " << tmax << G4endl;
  G4double transfer = 0.0;
  G4double rand = G4UniformRand();

  std::size_t nPlace = fParticleEnergyVector->GetVectorLength() - 1;
  std::size_t iPlace = fParticleEnergyVector->FindBin(scaledTkin, 0);

  G4bool one = true;
  if(scaledTkin >= fParticleEnergyVector->Energy(nPlace)) { iPlace = nPlace; }
  else if(scaledTkin > fParticleEnergyVector->Energy(0)) { 
    one = false; 
  }
  G4PhysicsTable* table = fPAIxscBank[coupleIndex];
  G4PhysicsVector* v1 = (*table)[iPlace];

  G4double emin = std::max(tmin, v1->Energy(0));
  G4double emax = std::min(tmax, v1->GetMaxEnergy());
  if(emax < emin) { return transfer; }

  G4double dNdx1 = v1->Value(emin)/emin;  
  G4double dNdx2 = v1->Value(emax)/emax;  
  /*
  G4cout << "iPlace= " << iPlace << " nPlace= " << nPlace 
	 << "  emin= " << emin << "  emax= " << emax 
	 << " dNdx1= " << dNdx1 << " dNdx2= " << dNdx2
	 << " one: " << one << G4endl;
  */
  G4double position = dNdx2 + (dNdx1 - dNdx2)*rand;
  transfer = GetEnergyTransfer(coupleIndex, iPlace, position);

  //G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"
  //	<< " position= " << position << G4endl; 

  if(!one) {

    G4PhysicsVector* v2 = (*table)[iPlace+1];
    emin = std::max(tmin, v2->Energy(0));
    emax = std::min(tmax, v2->GetMaxEnergy());
    if(emin <= emax) {
      dNdx1 = v2->Value(emin)/emin;  
      dNdx2 = v2->Value(emax)/emax;  

      //G4cout << "  emax2= " << emax 
      //     << " dNdx2= " << dNdx2 << " dNdx1= " << dNdx1 << G4endl;

      G4double E1 = fParticleEnergyVector->Energy(iPlace); 
      G4double E2 = fParticleEnergyVector->Energy(iPlace+1);
      G4double W  = 1.0/(E2 - E1);
      G4double W1 = (E2 - scaledTkin)*W;
      G4double W2 = (scaledTkin - E1)*W;
    
      //G4cout<< "E1= " << E1 << " E2= " << E2 <<" iPlace= " << iPlace 
      //    << " W1= " << W1 << " W2= " << W2 <<G4endl;
    
      position = dNdx2 + (dNdx1 - dNdx2)*rand;
      G4double tr2 = GetEnergyTransfer(coupleIndex, iPlace+1, position);

      //G4cout<<"PAImodel PostStepTransfer1 = "<<tr2/keV<<" keV"
      //    << " position= " << position << G4endl; 
      transfer *= W1;
      transfer += tr2*W2;
    }
  }
  //G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"
  //	<< " position= " << position << G4endl; 
  transfer = std::max(transfer, 0.0);
  return transfer;
}

///////////////////////////////////////////////////////////////////////
//
// Returns PAI energy transfer according to passed 
// indexes of particle kinetic enegry and random x-section

G4double G4PAIModelData::GetEnergyTransfer(G4int coupleIndex, 
					   std::size_t iPlace, 
					   G4double position) const
{ 
  G4PhysicsVector* v = (*(fPAIxscBank[coupleIndex]))(iPlace); 
  if(position*v->Energy(0) >= (*v)[0]) { return v->Energy(0); }

  std::size_t iTransferMax = v->GetVectorLength() - 1;

  std::size_t iTransfer;
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
  	 << " x1= " << x1 << " x2= " << x2 
  	 << " y1= " << y1 << " y2= " << y2 << G4endl;
  */
  energyTransfer = x1;
  if ( x1 != x2 ) {
    if ( y1 == y2  ) {
      energyTransfer += (x2 - x1)*G4UniformRand();
    } else {
      if(x1*1.1 < x2) {
	const G4int nbins = 5;
        G4double del = (x2 - x1)/G4int(nbins);
        x2 = x1;
        for(G4int i=1; i<=nbins; ++i) {
          x2 += del;
          y2 = v->Value(x2)/x2;
          if(position >= y2) { 
	    break; 
	  }
          x1 = x2;
          y1 = y2;
	}
      }
      //G4cout << "x1(keV)= " << x1/keV << " x2(keV)= " << x2/keV
      //   << " y1= " << y1 << " y2= " << y2 << " pos= " << position << G4endl;
      energyTransfer = (y2 - y1)*x1*x2/(position*(x1 - x2) - y1*x1 + y2*x2);
    }
  }
  //G4cout << "x1(keV)= " << x1/keV << " x2(keV)= " << x2/keV
  //	 << " y1= " << y1 << " y2= " << y2 << " pos= " << position
  //	 << " E(keV)= " << energyTransfer/keV << G4endl; 
  return energyTransfer;
}

//////////////////////////////////////////////////////////////////////


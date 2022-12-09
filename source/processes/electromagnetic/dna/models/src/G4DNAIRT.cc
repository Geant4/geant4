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
/*
 * G4DNAIRT.cc
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
*/


#include "G4DNAIRT.hh"
#include "G4ErrorFunction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4ITReactionChange.hh"
#include "G4ITTrackHolder.hh"
#include "G4ITReaction.hh"
#include "G4Scheduler.hh"

using namespace std;

G4DNAIRT::G4DNAIRT()  :
G4VITReactionProcess(),
fMolReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable)),
fpReactionModel(nullptr),
fTrackHolder(G4ITTrackHolder::Instance()),
fReactionSet(nullptr)
{
  timeMin = G4Scheduler::Instance()->GetStartTime();
  timeMax = G4Scheduler::Instance()->GetEndTime();

  fXMin = 1e9*nm;
  fYMin = 1e9*nm;
  fZMin = 1e9*nm;

  fXMax = 0e0*nm;
  fYMax = 0e0*nm;
  fZMax = 0e0*nm;

  fNx = 0;
  fNy = 0;
  fNz = 0;

  xiniIndex = 0, yiniIndex = 0, ziniIndex = 0;
  xendIndex = 0, yendIndex = 0, zendIndex = 0;

  fRCutOff =
    1.45 * nm + 2 * std::sqrt(8*9.46e9*nm*nm/s * timeMax); // 95% confidence level

  erfc = new G4ErrorFunction();
}


G4DNAIRT::G4DNAIRT(G4VDNAReactionModel* pReactionModel)
  : G4DNAIRT()
{
  fpReactionModel = pReactionModel;
}

G4DNAIRT::~G4DNAIRT()
{
  delete erfc;
}

void G4DNAIRT::Initialize(){

  fTrackHolder = G4ITTrackHolder::Instance();

  fReactionSet = G4ITReactionSet::Instance();
  fReactionSet->CleanAllReaction();
  fReactionSet->SortByTime();

  spaceBinned.clear();

  timeMin = G4Scheduler::Instance()->GetStartTime();
  timeMax = G4Scheduler::Instance()->GetEndTime();

  xiniIndex = 0;
  yiniIndex = 0;
  ziniIndex = 0;
  xendIndex = 0;
  yendIndex = 0;
  zendIndex = 0;

  fXMin = 1e9*nm;
  fYMin = 1e9*nm;
  fZMin = 1e9*nm;

  fXMax = 0e0*nm;
  fYMax = 0e0*nm;
  fZMax = 0e0*nm;

  fNx = 0;
  fNy = 0;
  fNz = 0;

  SpaceBinning();		// 1. binning the space
  IRTSampling();		// 2. Sampling of the IRT

  //hoang : if the first IRTSampling won't give any reactions, end the simu.
  if(fReactionSet->Empty())
  {
    for (auto pTrack : *fTrackHolder->GetMainList())
    {
      pTrack->SetGlobalTime(G4Scheduler::Instance()->GetEndTime());
    }
  }
}

void G4DNAIRT::SpaceBinning(){
  auto it_begin = fTrackHolder->GetMainList()->begin();
  while(it_begin != fTrackHolder->GetMainList()->end()){

    G4ThreeVector position = it_begin->GetPosition();

    if ( fXMin > position.x() ) fXMin = position.x(); 
    if ( fYMin > position.y() ) fYMin = position.y();
    if ( fZMin > position.z() ) fZMin = position.z();

    if ( fXMax < position.x() ) fXMax = position.x();
    if ( fYMax < position.y() ) fYMax = position.y();
    if ( fZMax < position.z() ) fZMax = position.z();

    ++it_begin;
  }

  fNx = G4int((fXMax-fXMin)/fRCutOff) == 0 ? 1 : G4int((fXMax-fXMin)/fRCutOff);
  fNy = G4int((fYMax-fYMin)/fRCutOff) == 0 ? 1 : G4int((fYMax-fYMin)/fRCutOff);
  fNz = G4int((fZMax-fZMin)/fRCutOff) == 0 ? 1 : G4int((fZMax-fZMin)/fRCutOff);

}

void G4DNAIRT::IRTSampling(){

  auto it_begin = fTrackHolder->GetMainList()->begin();
  while(it_begin != fTrackHolder->GetMainList()->end()){
    G4int I = FindBin(fNx, fXMin, fXMax, it_begin->GetPosition().x());
    G4int J = FindBin(fNy, fYMin, fYMax, it_begin->GetPosition().y());
    G4int K = FindBin(fNz, fZMin, fZMax, it_begin->GetPosition().z());

    spaceBinned[I][J][K].push_back(*it_begin);

    Sampling(*it_begin);
    ++it_begin;
  }
}

void G4DNAIRT::Sampling(G4Track* track){
  G4Molecule* molA = G4Molecule::GetMolecule(track);
  const G4MolecularConfiguration* molConfA = molA->GetMolecularConfiguration();
  if(molConfA->GetDiffusionCoefficient() == 0) return;

  const vector<const G4MolecularConfiguration*>* reactivesVector =
    fMolReactionTable->CanReactWith(molConfA);

  if(reactivesVector == nullptr) return;

  G4double globalTime = G4Scheduler::Instance()->GetGlobalTime();
  G4double minTime = timeMax;

  xiniIndex = FindBin(fNx, fXMin, fXMax, track->GetPosition().x()-fRCutOff);
  xendIndex = FindBin(fNx, fXMin, fXMax, track->GetPosition().x()+fRCutOff);
  yiniIndex = FindBin(fNy, fYMin, fYMax, track->GetPosition().y()-fRCutOff);
  yendIndex = FindBin(fNy, fYMin, fYMax, track->GetPosition().y()+fRCutOff);
  ziniIndex = FindBin(fNz, fZMin, fZMax, track->GetPosition().z()-fRCutOff);
  zendIndex = FindBin(fNz, fZMin, fZMax, track->GetPosition().z()+fRCutOff);

  for ( int ii = xiniIndex; ii <= xendIndex; ii++ ) {
    for ( int jj = yiniIndex; jj <= yendIndex; jj++ ) {
      for ( int kk = ziniIndex; kk <= zendIndex; kk++ ) {

        std::vector<G4Track*> spaceBin = spaceBinned[ii][jj][kk];
        for ( int n = 0; n < (int)spaceBinned[ii][jj][kk].size(); n++ ) {
          if(!spaceBin[n] || track == spaceBin[n]) continue;
          if(spaceBin[n]->GetTrackStatus() == fStopButAlive) continue;

          G4Molecule* molB = G4Molecule::GetMolecule(spaceBin[n]);
          if(!molB) continue;

          const G4MolecularConfiguration* molConfB = molB->GetMolecularConfiguration();
          if(molConfB->GetDiffusionCoefficient() == 0) continue;

          auto it = std::find(reactivesVector->begin(), reactivesVector->end(), molConfB);
          if(it == reactivesVector->end()) continue;

          G4ThreeVector orgPosB = spaceBin[n]->GetPosition();
          G4double dt = track->GetGlobalTime() - spaceBin[n]->GetGlobalTime();
          G4ThreeVector newPosB = orgPosB;

          if(dt > 0){
            G4double sigma, x, y, z;
            G4double diffusionCoefficient = G4Molecule::GetMolecule(spaceBin[n])->GetDiffusionCoefficient();

            sigma = std::sqrt(2.0 * diffusionCoefficient * dt);

            x = G4RandGauss::shoot(0., 1.0)*sigma;
            y = G4RandGauss::shoot(0., 1.0)*sigma;
            z = G4RandGauss::shoot(0., 1.0)*sigma;

            newPosB = orgPosB + G4ThreeVector(x,y,z);
          }else if(dt < 0) continue;

          G4double r0 = (newPosB - track->GetPosition()).mag();
          G4double irt = GetIndependentReactionTime(molConfA,
                                                    molConfB,
                                                    r0);
          if(irt>=0 && irt<timeMax - globalTime)
          {
            irt += globalTime;
            if(irt < minTime) minTime = irt;
#ifdef DEBUG
            G4cout<<irt<<'\t'<<molConfA->GetName()<<" "<<track->GetTrackID()<<'\t'<<molConfB->GetName()<<" "<<spaceBin[n]->GetTrackID()<<'\n';
#endif
            fReactionSet->AddReaction(irt,track,spaceBin[n]);
          }
        }
        spaceBin.clear();
      }
    }
  }

// Scavenging & first order reactions

  auto fReactionDatas = fMolReactionTable->GetReactionData(molConfA);
  G4double index = -1;
  //change the scavenging filter of the IRT beyond 1 us proposed by Naoki and Jose
  if(timeMax > 1*us)
  {
    minTime = timeMax;
  }
  //

  for(size_t u=0; u<fReactionDatas->size();u++){
    if((*fReactionDatas)[u]->GetReactant2()->GetDiffusionCoefficient() == 0){
      G4double kObs = (*fReactionDatas)[u]->GetObservedReactionRateConstant();
      if(kObs == 0) continue;
      G4double time = -(std::log(1.0 - G4UniformRand())/kObs) + globalTime;
      if( time < minTime && time >= globalTime && time < timeMax){
        minTime = time;
        index = (G4int)u;
      }
    }
  }

  if(index != -1){
#ifdef DEBUG
    G4cout<<"scavenged: "<<minTime<<'\t'<<molConfA->GetName()<<it_begin->GetTrackID()<<'\n';
#endif
    G4Molecule* fakeMol = new G4Molecule((*fReactionDatas)[index]->GetReactant2());
    G4Track* fakeTrack = fakeMol->BuildTrack(globalTime,track->GetPosition());
    fTrackHolder->Push(fakeTrack);
    fReactionSet->AddReaction(minTime, track, fakeTrack);
  }
}


G4double G4DNAIRT::GetIndependentReactionTime(const G4MolecularConfiguration* molA, const G4MolecularConfiguration* molB, G4double distance) {
  const auto pMoleculeA = molA;
  const auto pMoleculeB = molB;
  auto fReactionData = fMolReactionTable->GetReactionData(pMoleculeA, pMoleculeB);
  G4int reactionType = fReactionData->GetReactionType();
  G4double r0 = distance;
  if(r0 == 0) r0 += 1e-3*nm;
  G4double irt = -1 * ps;
  G4double D = molA->GetDiffusionCoefficient() +
               molB->GetDiffusionCoefficient();
  if(D == 0) D += 1e-20*(m2/s);
  G4double rc = fReactionData->GetOnsagerRadius();

  if ( reactionType == 0){
    G4double sigma = fReactionData->GetEffectiveReactionRadius();

    if(sigma > r0) return 0; // contact reaction
    if( rc != 0) r0 = -rc / (1-std::exp(rc/r0));

    G4double Winf = sigma/r0;
    G4double W = G4UniformRand();

    if ( W > 0 && W < Winf ) irt = (0.25/D) * std::pow( (r0-sigma)/erfc->erfcInv(r0*W/sigma), 2 );

    return irt;
  }
  else if ( reactionType == 1 ){
    G4double sigma = fReactionData->GetReactionRadius();
    G4double kact = fReactionData->GetActivationRateConstant();
    G4double kdif = fReactionData->GetDiffusionRateConstant();
    G4double kobs = fReactionData->GetObservedReactionRateConstant();

    G4double a, b, Winf;

    if ( rc == 0 ) {
      a = 1/sigma * kact / kobs;
      b = (r0 - sigma) / 2;
    } else {
      G4double v = kact/Avogadro/(4*CLHEP::pi*pow(sigma,2) * exp(-rc / sigma));
      G4double alpha = v+rc*D/(pow(sigma,2)*(1-exp(-rc/sigma)));
      a = 4*pow(sigma,2)*alpha/(D*pow(rc,2))*pow(sinh(rc/(2*sigma)),2);
      b = rc/4*(cosh(rc/(2*r0))/sinh(rc/(2*r0))-cosh(rc/(2*sigma))/sinh(rc/(2*sigma)));
      r0 = -rc/(1-std::exp(rc/r0));
      sigma = fReactionData->GetEffectiveReactionRadius();
    }

    if(sigma > r0){
      if(fReactionData->GetProbability() > G4UniformRand()) return 0;
      else return irt;
    }
    Winf = sigma / r0 * kobs / kdif;

    if(Winf > G4UniformRand()) irt = SamplePDC(a,b)/D;
    return irt;
  }

  return -1 * ps;
}

G4int G4DNAIRT::FindBin(G4int n, G4double xmin, G4double xmax, G4double value) {

  G4int bin = -1;
  if ( value <= xmin )
    bin = 0; //1;
  else if ( value >= xmax)  //!(xmax < value) ) //value >= xmax )
    bin = n-1; //n;
  else
    bin = G4int( n * ( value - xmin )/( xmax - xmin ) ); //bin = 1 + G4int( n * ( value - xmin )/( xmax - xmin ) );

  if ( bin < 0 ) bin = 0;
  if ( bin >= n ) bin = n-1;

  return bin;
}

G4double G4DNAIRT::SamplePDC(G4double a, G4double b) {

  G4double p = 2.0 * std::sqrt(2.0*b/a);
  G4double q = 2.0 / std::sqrt(2.0*b/a);
  G4double M = max(1.0/(a*a),3.0*b/a);

  G4double X, U, lambdax;

  G4int ntrials = 0;
  while(1) {

    // Generate X
    U = G4UniformRand();
    if ( U < p/(p + q * M) ) X = pow(U * (p + q * M) / 2, 2);
    else X = pow(2/((1-U)*(p+q*M)/M),2);

    U = G4UniformRand();

    lambdax = std::exp(-b*b/X) * ( 1.0 - a * std::sqrt(CLHEP::pi * X) * erfc->erfcx(b/std::sqrt(X) + a*std::sqrt(X)));

    if ((X <= 2.0*b/a && U <= lambdax) ||
        (X >= 2.0*b/a && U*M/X <= lambdax)) break;

    ntrials++;

    if ( ntrials > 10000 ){
      G4cout<<"Totally rejected"<<'\n';
      return -1.0;
    }
  }
  return X;
}

std::unique_ptr<G4ITReactionChange> G4DNAIRT::MakeReaction(const G4Track& trackA,
                                                         const G4Track& trackB)
{

  std::unique_ptr<G4ITReactionChange> pChanges(new G4ITReactionChange());
  pChanges->Initialize(trackA, trackB);

  const auto pMoleculeA = GetMolecule(trackA)->GetMolecularConfiguration();
  const auto pMoleculeB = GetMolecule(trackB)->GetMolecularConfiguration();
  const auto pReactionData = fMolReactionTable->GetReactionData(pMoleculeA, pMoleculeB);

  G4double globalTime = G4Scheduler::Instance()->GetGlobalTime();
  G4double effectiveReactionRadius = pReactionData->GetEffectiveReactionRadius();

  const G4double D1 = pMoleculeA->GetDiffusionCoefficient();
  const G4double D2 = pMoleculeB->GetDiffusionCoefficient();

  G4ThreeVector r1 = trackA.GetPosition();
  G4ThreeVector r2 = trackB.GetPosition();

  if(r1 == r2) r2 += G4ThreeVector(0,0,1e-3*nm);

  G4ThreeVector S1 = r1 - r2;

  G4double r0 = S1.mag();

  S1.setMag(effectiveReactionRadius);

  G4double dt = globalTime - trackA.GetGlobalTime();

  if(dt != 0 && (D1 + D2) != 0 && r0 != 0){
    G4double s12 = 2.0 * D1 * dt;
    G4double s22 = 2.0 * D2 * dt;
    if(s12 == 0) r2 = r1;
    else if(s22 == 0) r1 = r2;
    else{
      G4double alpha = effectiveReactionRadius * r0 / (2*(D1 + D2)*dt);
      G4ThreeVector S2 = (r1 + (s12 / s22)*r2) + G4ThreeVector(G4RandGauss::shoot(0, s12 + s22 * s22 / s12),
                                                               G4RandGauss::shoot(0, s12 + s22 * s22 / s12),
                                                               G4RandGauss::shoot(0, s12 + s22 * s22 / s12));

      if(alpha == 0){
        return pChanges;
      }
      S1.setPhi(rad * G4UniformRand() * 2.0 * CLHEP::pi);
      S1.setTheta(rad * std::acos(1.0 + 1./alpha * std::log(1.0 - G4UniformRand() * (1 - std::exp(-2.0 * alpha)))));

      r1 = (D1 * S1 + D2 * S2) / (D1 + D2);
      r2 = D2 * (S2 - S1) / (D1 + D2);
    }
  }

  auto pTrackA = const_cast<G4Track*>(pChanges->GetTrackA());
  auto pTrackB = const_cast<G4Track*>(pChanges->GetTrackB());

  pTrackA->SetPosition(r1);
  pTrackB->SetPosition(r2);

  pTrackA->SetGlobalTime(globalTime);
  pTrackB->SetGlobalTime(globalTime);

  pTrackA->SetTrackStatus(fStopButAlive);
  pTrackB->SetTrackStatus(fStopButAlive);

  const G4int nbProducts = pReactionData->GetNbProducts();

  if(nbProducts){

    const G4double sqrD1 = D1 == 0. ? 0. : std::sqrt(D1);
    const G4double sqrD2 = D2 == 0. ? 0. : std::sqrt(D2);
    if((sqrD1 + sqrD2) == 0){
      return pChanges;
    }
    const G4double inv_numerator = 1./(sqrD1 + sqrD2);
    const G4ThreeVector reactionSite = sqrD2 * inv_numerator * trackA.GetPosition()
                                     + sqrD1 * inv_numerator * trackB.GetPosition();

    std::vector<G4ThreeVector> position;

    if(nbProducts == 1){
      position.push_back(reactionSite);
    }else if(nbProducts == 2){
      position.push_back(trackA.GetPosition());
      position.push_back(trackB.GetPosition());
    }else if (nbProducts == 3){
      position.push_back(reactionSite);
      position.push_back(trackA.GetPosition());
      position.push_back(trackB.GetPosition());
    }

    for(G4int u = 0; u < nbProducts; u++){

      auto product = new G4Molecule(pReactionData->GetProduct(u));
      auto productTrack = product->BuildTrack(globalTime,
                                              position[u]);

      productTrack->SetTrackStatus(fAlive);
      fTrackHolder->Push(productTrack);

      pChanges->AddSecondary(productTrack);

      G4int I = FindBin(fNx, fXMin, fXMax, position[u].x());
      G4int J = FindBin(fNy, fYMin, fYMax, position[u].y());
      G4int K = FindBin(fNz, fZMin, fZMax, position[u].z());

      spaceBinned[I][J][K].push_back(productTrack);

      Sampling(productTrack);
    }
  }

  fTrackHolder->MergeSecondariesWithMainList();
  pChanges->KillParents(true);
  return pChanges;
}


std::vector<std::unique_ptr<G4ITReactionChange>> G4DNAIRT::FindReaction(
  G4ITReactionSet* pReactionSet,
  const G4double /*currentStepTime*/,
  const G4double fGlobalTime,
  const G4bool /*reachedUserStepTimeLimit*/)
{
  std::vector<std::unique_ptr<G4ITReactionChange>> fReactionInfo;
  fReactionInfo.clear();

  if (pReactionSet == nullptr)
  {
    return fReactionInfo;
  }

  auto fReactionsetInTime = pReactionSet->GetReactionsPerTime();
  assert(fReactionsetInTime.begin() != fReactionsetInTime.end());

  auto it_begin = fReactionsetInTime.begin();
  while(it_begin != fReactionsetInTime.end())
  {
    G4double irt = it_begin->get()->GetTime();

    if(fGlobalTime < irt) break;

    pReactionSet->SelectThisReaction(*it_begin);

    G4Track* pTrackA = it_begin->get()->GetReactants().first;
    G4Track* pTrackB = it_begin->get()->GetReactants().second;
    auto pReactionChange = MakeReaction(*pTrackA, *pTrackB);

    if(pReactionChange){
      fReactionInfo.push_back(std::move(pReactionChange));
    }

    fReactionsetInTime = pReactionSet->GetReactionsPerTime();
    it_begin = fReactionsetInTime.begin();
  }

  return fReactionInfo;
}

G4bool G4DNAIRT::TestReactibility(const G4Track& /*trackA*/,
  const G4Track& /*trackB*/,
  G4double /*currentStepTime*/,
  G4bool /*userStepTimeLimit*/) /*const*/
{
  return true;
}

void G4DNAIRT::SetReactionModel(G4VDNAReactionModel* model)
{
  fpReactionModel = model;
}

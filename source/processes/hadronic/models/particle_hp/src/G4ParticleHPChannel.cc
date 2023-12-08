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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 070523 bug fix for G4FPE_DEBUG on by A. Howard ( and T. Koi)
// 071031 bug fix T. Koi on behalf of A. Howard
// 081203 bug fix in Register method by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
// June-2019 - E. Mendoza --> Modification to allow using an incomplete
//   data library if the G4NEUTRONHP_SKIP_MISSING_ISOTOPES environmental
//   flag is defined. The missing XS are set to 0.

#include "G4ParticleHPChannel.hh"

#include "G4HadTmpUtil.hh"
#include "G4ParticleHPElasticFS.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPReactionWhiteBoard.hh"
#include "G4ParticleHPThermalBoost.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include <cstdlib>

G4ParticleHPChannel::G4ParticleHPChannel(G4ParticleDefinition* p)
{
  fManager = G4ParticleHPManager::GetInstance();
  if (fManager->GetUseWendtFissionModel()) {
    wendtFissionGenerator = G4WendtFissionFragmentGenerator::GetInstance();
    // Make sure both fission fragment models are not active at same time
    fManager->SetProduceFissionFragments(false);
  }
  theProjectile = (nullptr == p) ? G4Neutron::Neutron() : p;
  theChannelData = new G4ParticleHPVector;
}

G4ParticleHPChannel::~G4ParticleHPChannel()
{
  delete theChannelData;
  // Following statement disabled to avoid SEGV
  // theBuffer is also deleted as "theChannelData" in
  delete[] theIsotopeWiseData;
  if (theFinalStates != nullptr) {
    for (G4int i = 0; i < niso; i++) {
      delete theFinalStates[i];
    }
    delete[] theFinalStates;
  }
  delete[] active;
}

G4double G4ParticleHPChannel::GetXsec(G4double energy)
{
  return std::max(0., theChannelData->GetXsec(energy));
}

G4double G4ParticleHPChannel::GetWeightedXsec(G4double energy,
					      G4int isoNumber)
{
  return theIsotopeWiseData[isoNumber].GetXsec(energy);
}

G4double G4ParticleHPChannel::GetFSCrossSection(G4double energy,
						G4int isoNumber)
{
  return theFinalStates[isoNumber]->GetXsec(energy);
}

void G4ParticleHPChannel::Init(G4Element* anElement, 
			       const G4String& dirName, const G4String& aFSType)
{
  theFSType = aFSType;
  Init(anElement, dirName);
}

void G4ParticleHPChannel::Init(G4Element* anElement, const G4String& dirName)
{
  theDir = dirName;
  theElement = anElement;
}

G4bool G4ParticleHPChannel::Register(G4ParticleHPFinalState* theFS)
{
  ++registerCount;
  G4int Z = theElement->GetZasInt();

  niso = (G4int)theElement->GetNumberOfIsotopes();

  delete[] theIsotopeWiseData;
  theIsotopeWiseData = new G4ParticleHPIsoData[niso];
  delete[] active;
  active = new G4bool[niso];

  delete[] theFinalStates;
  theFinalStates = new G4ParticleHPFinalState*[niso];
  delete theChannelData;
  theChannelData = new G4ParticleHPVector;
  for (G4int i = 0; i < niso; ++i) {
    theFinalStates[i] = theFS->New();
    theFinalStates[i]->SetProjectile(theProjectile);
  }
  if (niso != 0 && registerCount == 0) {
    for (G4int i1 = 0; i1 < niso; ++i1) {
      G4int A = theElement->GetIsotope(i1)->GetN();
      G4int M = theElement->GetIsotope(i1)->Getm();
      //G4cout <<" Init: normal case i=" << i1 
      //     << " Z=" << Z << " A=" << A << G4endl;
      G4double frac = theElement->GetRelativeAbundanceVector()[i1] / perCent;
      theFinalStates[i1]->SetA_Z(A, Z, M);
      UpdateData(A, Z, M, i1, frac, theProjectile);
    }
  }
  G4bool result = HasDataInAnyFinalState();

  // To avoid issuing hash by worker threads
  if (result) theChannelData->Hash();

  return result;
}

void G4ParticleHPChannel::UpdateData(G4int A, G4int Z, G4int M, G4int index,
                                     G4double abundance,
                                     G4ParticleDefinition* projectile)
{
  // Initialze the G4FissionFragment generator for this isomer if needed
  if (wendtFissionGenerator != nullptr) {
    wendtFissionGenerator->InitializeANucleus(A, Z, M, theDir);
  }

  theFinalStates[index]->Init(A, Z, M, theDir, theFSType, projectile);
  if (!theFinalStates[index]->HasAnyData()) return;
  // nothing there for exactly this isotope.

  // the above has put the X-sec into the FS
  theBuffer = nullptr;
  if (theFinalStates[index]->HasXsec()) {
    theBuffer = theFinalStates[index]->GetXsec();
    theBuffer->Times(abundance / 100.);
    theIsotopeWiseData[index].FillChannelData(theBuffer);
  }
  else  // get data from CrossSection directory
  {
    G4String tString = "/CrossSection";
    active[index] = theIsotopeWiseData[index].Init(A, Z, M, abundance,
                                                   theDir, tString);
    if (active[index]) theBuffer = theIsotopeWiseData[index].MakeChannelData();
  }
  if (theBuffer != nullptr) Harmonise(theChannelData, theBuffer);
}

void G4ParticleHPChannel::Harmonise(G4ParticleHPVector*& theStore,
                                    G4ParticleHPVector* theNew)
{
  G4int s_tmp = 0, n = 0, m_tmp = 0;
  auto theMerge = new G4ParticleHPVector;
  G4ParticleHPVector* anActive = theStore;
  G4ParticleHPVector* aPassive = theNew;
  G4ParticleHPVector* tmp;
  G4int a = s_tmp, p = n, t;
  while (a < anActive->GetVectorLength() && p < aPassive->GetVectorLength())
    // Loop checking, 11.05.2015, T. Koi
  {
    if (anActive->GetEnergy(a) <= aPassive->GetEnergy(p)) {
      G4double xa = anActive->GetEnergy(a);
      theMerge->SetData(m_tmp, xa, anActive->GetXsec(a) + std::max(0., aPassive->GetXsec(xa)));
      m_tmp++;
      a++;
      G4double xp = aPassive->GetEnergy(p);
      if (std::abs(std::abs(xp - xa) / xa) < 0.001) {
        ++p;
      }
    }
    else {
      tmp = anActive;
      t = a;
      anActive = aPassive;
      a = p;
      aPassive = tmp;
      p = t;
    }
  }
  while (a != anActive->GetVectorLength())  // Loop checking, 11.05.2015, T. Koi
  {
    theMerge->SetData(m_tmp++, anActive->GetEnergy(a), anActive->GetXsec(a));
    ++a;
  }
  while (p != aPassive->GetVectorLength())  // Loop checking, 11.05.2015, T. Koi
  {
    if (std::abs(theMerge->GetEnergy(std::max(0, m_tmp - 1)) -
		 aPassive->GetEnergy(p)) / aPassive->GetEnergy(p) > 0.001)
      theMerge->SetData(m_tmp++, aPassive->GetEnergy(p), aPassive->GetXsec(p));
    ++p;
  }
  delete theStore;
  theStore = theMerge;
}

G4HadFinalState*
G4ParticleHPChannel::ApplyYourself(const G4HadProjectile& theTrack,
				   G4int anIsotope, G4bool isElastic)
{
  //G4cout << "G4ParticleHPChannel::ApplyYourself niso=" << niso
  //	 << " ni=" << anIsotope << " isElastic=" << isElastic <<G4endl;
  if (anIsotope != -1 && anIsotope != -2) {
    // Inelastic Case
    //G4cout << "G4ParticleHPChannel Inelastic Case"
    //<< " Z= " << GetZ(anIsotope) << " A = " << GetN(anIsotope) << G4endl;
    fManager->GetReactionWhiteBoard()->SetTargA((G4int)GetN(anIsotope));
    fManager->GetReactionWhiteBoard()->SetTargZ((G4int)GetZ(anIsotope));
    return theFinalStates[anIsotope]->ApplyYourself(theTrack);
  }
  G4double sum = 0;
  G4int it = 0;
  auto xsec = new G4double[niso];
  G4ParticleHPThermalBoost aThermalE;
  for (G4int i = 0; i < niso; i++) {
    if (theFinalStates[i]->HasAnyData()) {
      /*
      G4cout << "FS: " << i << theTrack.GetDefinition()->GetParticleName()
	     << " Z=" << theFinalStates[i]->GetZ() 
	     << " A=" << theFinalStates[i]->GetN() 
	     << G4endl;
      */
      xsec[i] = theIsotopeWiseData[i].GetXsec(
        aThermalE.GetThermalEnergy(theTrack, theFinalStates[i]->GetN(),
                                   theFinalStates[i]->GetZ(),
                                   theTrack.GetMaterial()->GetTemperature()));
      sum += xsec[i];
    }
    else {
      xsec[i] = 0;
    }
  }
  if (sum == 0) {
    it = G4lrint(niso * G4UniformRand());
  }
  else {
    G4double random = G4UniformRand();
    G4double running = 0;
    for (G4int ix = 0; ix < niso; ix++) {
      running += xsec[ix];
      if (sum == 0 || random <= running / sum) {
        it = ix;
        break;
      }
    }
    if (it == niso) it--;
  }
  delete[] xsec;
  G4HadFinalState* theFinalState = nullptr;
  const auto A = (G4int)this->GetN(it);
  const auto Z = (G4int)this->GetZ(it);
  const auto M = (G4int)this->GetM(it);

  //-2:Marker for Fission
  if ((wendtFissionGenerator != nullptr) && anIsotope == -2) {
    theFinalState = wendtFissionGenerator->ApplyYourself(theTrack, Z, A);
  }

  // Use the standard procedure if the G4FissionFragmentGenerator model fails
  if (theFinalState == nullptr) {
    G4int icounter = 0;
    G4int icounter_max = 1024;
    while (theFinalState == nullptr)  // Loop checking, 11.05.2015, T. Koi
    {
      icounter++;
      if (icounter > icounter_max) {
        G4cout << "Loop-counter exceeded the threshold value at " 
               << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
        break;
      }
      if (isElastic) {
        // Register 0 K cross-section for DBRC for Doppler broadened elastic scattering kernel
        G4ParticleHPVector* xsforFS = theIsotopeWiseData[it].MakeChannelData();
        // Only G4ParticleHPElasticFS has the RegisterCrossSection method
        static_cast<G4ParticleHPElasticFS*>(theFinalStates[it])->RegisterCrossSection(xsforFS);
      }
      theFinalState = theFinalStates[it]->ApplyYourself(theTrack);
    }
  }

  // G4cout <<"THE IMPORTANT RETURN"<<G4endl;
  // G4cout << "TK G4ParticleHPChannel Elastic, Capture and Fission Cases "
  //<< " Z= " << this->GetZ(it) << " A = " << this->GetN(it) << G4endl;
  fManager->GetReactionWhiteBoard()->SetTargA(A);
  fManager->GetReactionWhiteBoard()->SetTargZ(Z);
  fManager->GetReactionWhiteBoard()->SetTargM(M);

  return theFinalState;
}

void G4ParticleHPChannel::DumpInfo()
{
  G4cout << " Element: " << theElement->GetName() << G4endl;
  G4cout << " Directory name: " << theDir << G4endl;
  G4cout << " FS name: " << theFSType << G4endl;
  G4cout << " Number of Isotopes: " << niso << G4endl;
  G4cout << " Have cross sections: " << G4endl;
  for (int i = 0; i < niso; i++) {
    G4cout << theFinalStates[i]->HasXsec() << "  ";
  }
  G4cout << G4endl;
  if (theChannelData != nullptr) {
    G4cout << " Cross Section (total for this channel):" << G4endl;
    int np = theChannelData->GetVectorLength();
    G4cout << np << G4endl;
    for (int i = 0; i < np; i++) {
      G4cout << theChannelData->GetEnergy(i) / eV << "  " << theChannelData->GetXsec(i) << G4endl;
    }
  }
}

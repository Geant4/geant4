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
// $Id: G4MuBremsstrahlungModel.cc 97392 2016-06-02 10:10:32Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4MuBremsstrahlungModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 24.06.2002
//
// Modifications:
//
// 04-12-02 Change G4DynamicParticle constructor in PostStepDoIt (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 24-01-03 Fix for compounds (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 10-02-04 Add lowestKinEnergy (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 03-08-05 Angular correlations according to PRM (V.Ivanchenko)
// 13-02-06 add ComputeCrossSectionPerAtom (mma)
// 21-03-06 Fix problem of initialisation in case when cuts are not defined (VI)
// 07-11-07 Improve sampling of final state (A.Bogdanov)
// 28-02-08 Use precomputed Z^1/3 and Log(A) (V.Ivanchenko)
// 31-05-13 Use element selectors instead of local data structure (V.Ivanchenko)
//

//
// Class Description:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuBremsstrahlungModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

const G4double G4MuBremsstrahlungModel::xgi[] = 
  {0.03377,0.16940,0.38069,0.61931,0.83060,0.96623};
const G4double G4MuBremsstrahlungModel::wgi[] = 
  {0.08566,0.18038,0.23396,0.23396,0.18038,0.08566};
G4double G4MuBremsstrahlungModel::fDN[] = {0.0};

G4MuBremsstrahlungModel::G4MuBremsstrahlungModel(const G4ParticleDefinition* p,
                                                 const G4String& nam)
  : G4VEmModel(nam),
    particle(nullptr),
    sqrte(sqrt(G4Exp(1.))),
    bh(202.4),
    bh1(446.),
    btf(183.),
    btf1(1429.),
    fParticleChange(nullptr),
    lowestKinEnergy(1.0*GeV),
    minThreshold(0.9*keV)
{
  theGamma = G4Gamma::Gamma();
  nist = G4NistManager::Instance();

  lowestKinEnergy = 1.*GeV;  

  mass = rmass = cc = coeff = 1.0;

  if(0.0 == fDN[1]) {
    for(G4int i=1; i<93; ++i) {
      G4double dn = 1.54*nist->GetA27(i);
      fDN[i] = dn;
      if(1 < i) {
        fDN[i] /= std::pow(dn, 1./G4double(i));
      }
    }
  }

  if(p) { SetParticle(p); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuBremsstrahlungModel::~G4MuBremsstrahlungModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBremsstrahlungModel::MinEnergyCut(const G4ParticleDefinition*,
                                               const G4MaterialCutsCouple*)
{
  return minThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBremsstrahlungModel::MinPrimaryEnergy(const G4Material*,
                                                   const G4ParticleDefinition*,
                                                   G4double cut)
{
  return std::max(lowestKinEnergy,cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuBremsstrahlungModel::Initialise(const G4ParticleDefinition* p,
                                         const G4DataVector& cuts)
{
  if(p) { SetParticle(p); }

  // define pointer to G4ParticleChange
  if(!fParticleChange) { fParticleChange = GetParticleChangeForLoss(); }

  if(IsMaster() && p == particle && lowestKinEnergy < HighEnergyLimit()) { 
    InitialiseElementSelectors(p, cuts); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuBremsstrahlungModel::InitialiseLocal(const G4ParticleDefinition* p,
                                              G4VEmModel* masterModel)
{
  if(p == particle && lowestKinEnergy < HighEnergyLimit()) {
    SetElementSelectors(masterModel->GetElementSelectors());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBremsstrahlungModel::ComputeDEDXPerVolume(
                                              const G4Material* material,
                                              const G4ParticleDefinition*,
                                                    G4double kineticEnergy,
                                                    G4double cutEnergy)
{
  G4double dedx = 0.0;
  if (kineticEnergy <= lowestKinEnergy) { return dedx; }

  G4double tmax = kineticEnergy;
  G4double cut  = std::min(cutEnergy,tmax);
  if(cut < minThreshold) { cut = minThreshold; }

  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector =
    material->GetAtomicNumDensityVector();

  //  loop for elements in the material
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4double loss = 
      ComputMuBremLoss((*theElementVector)[i]->GetZ(), kineticEnergy, cut);

    dedx += loss*theAtomicNumDensityVector[i];
  }
  //  G4cout << "BR e= " << kineticEnergy << "  dedx= " << dedx << G4endl;
  if(dedx < 0.) dedx = 0.;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBremsstrahlungModel::ComputMuBremLoss(G4double Z,
                                                   G4double tkin, G4double cut)
{
  G4double totalEnergy = mass + tkin;
  static const G4double ak1 = 0.05;
  static const G4int    k2=5;
  G4double loss = 0.;

  G4double vcut = cut/totalEnergy;
  G4double vmax = tkin/totalEnergy;

  G4double aaa = 0.;
  G4double bbb = vcut;
  if(vcut>vmax) { bbb = vmax; }
  G4int kkk = (G4int)((bbb-aaa)/ak1)+k2;
  if(kkk < 1) { kkk = 1; }

  G4double hhh=(bbb-aaa)/G4double(kkk);

  G4double aa = aaa;
  for(G4int l=0; l<kkk; l++)
  {
    for(G4int i=0; i<6; i++)
    {
      G4double ep = (aa + xgi[i]*hhh)*totalEnergy;
      loss += ep*wgi[i]*ComputeDMicroscopicCrossSection(tkin, Z, ep);
    }
    aa += hhh;
  }

  loss *=hhh*totalEnergy ;

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBremsstrahlungModel::ComputeMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double cut)
{
  G4double totalEnergy = tkin + mass;
  static const G4double ak1 = 2.3;
  static const G4int    k2  = 4;
  G4double cross = 0.;

  if(cut >= tkin) return cross;

  G4double vcut = cut/totalEnergy;
  G4double vmax = tkin/totalEnergy;

  G4double aaa = G4Log(vcut);
  G4double bbb = G4Log(vmax);
  G4int    kkk = (G4int)((bbb-aaa)/ak1)+k2 ;
  if(kkk < 1) { kkk = 1; }

  G4double hhh = (bbb-aaa)/G4double(kkk);

  G4double aa = aaa;

  for(G4int l=0; l<kkk; l++)
  {
    for(G4int i=0; i<6; i++)
    {
      G4double ep = G4Exp(aa + xgi[i]*hhh)*totalEnergy;
      cross += ep*wgi[i]*ComputeDMicroscopicCrossSection(tkin, Z, ep);
    }
    aa += hhh;
  }

  cross *=hhh;

  //G4cout << "BR e= " << tkin<< "  cross= " << cross/barn << G4endl;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBremsstrahlungModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double gammaEnergy)
//  differential cross section
{
  G4double dxsection = 0.;

  if(gammaEnergy > tkin) { return dxsection; }

  G4double E = tkin + mass ;
  G4double v = gammaEnergy/E ;
  G4double delta = 0.5*mass*mass*v/(E-gammaEnergy) ;
  G4double rab0  = delta*sqrte ;

  G4int iz = G4lrint(Z);
  if(iz < 1) { iz = 1; }
  else if(iz > 92) { iz = 92; }

  G4double z13 = 1.0/nist->GetZ13(iz);
  G4double dnstar = fDN[iz];

  G4double b,b1;

  if(1 == iz) {
    b  = bh;
    b1 = bh1;
  } else {
    b  = btf;
    b1 = btf1;
  }

  // nucleus contribution logarithm
  G4double rab1=b*z13;
  G4double fn=G4Log(rab1/(dnstar*(electron_mass_c2+rab0*rab1))*
              (mass+delta*(dnstar*sqrte-2.))) ;
  if(fn <0.) { fn = 0.; }
  // electron contribution logarithm
  G4double epmax1=E/(1.+0.5*mass*rmass/E) ;
  G4double fe=0.;
  if(gammaEnergy<epmax1)
  {
    G4double rab2=b1*z13*z13 ;
    fe=G4Log(rab2*mass/((1.+delta*rmass/(electron_mass_c2*sqrte))*
                              (electron_mass_c2+rab0*rab2))) ;
    if(fe<0.) { fe=0.; }
  }

  dxsection = coeff*(1.-v*(1. - 0.75*v))*Z*(fn*Z + fe)/gammaEnergy;

  return dxsection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBremsstrahlungModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition*,
                                                 G4double kineticEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double cross = 0.0;
  if (kineticEnergy <= lowestKinEnergy) return cross;
  G4double tmax = std::min(maxEnergy, kineticEnergy);
  G4double cut  = std::min(cutEnergy, kineticEnergy);
  if(cut < minThreshold) cut = minThreshold;
  if (cut >= tmax) return cross;

  cross = ComputeMicroscopicCrossSection (kineticEnergy, Z, cut);
  if(tmax < kineticEnergy) {
    cross -= ComputeMicroscopicCrossSection(kineticEnergy, Z, tmax);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuBremsstrahlungModel::SampleSecondaries(
                              std::vector<G4DynamicParticle*>* vdp,
                              const G4MaterialCutsCouple* couple,
                              const G4DynamicParticle* dp,
                              G4double minEnergy,
                              G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  // check against insufficient energy
  G4double tmax = std::min(kineticEnergy, maxEnergy);
  G4double tmin = std::min(kineticEnergy, minEnergy);
  if(tmin < minThreshold) tmin = minThreshold;
  if(tmin >= tmax) return;

  // ===== sampling of energy transfer ======

  G4ParticleMomentum partDirection = dp->GetMomentumDirection();

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(couple,particle,kineticEnergy);
  G4double Z = anElement->GetZ();

  G4double totalEnergy   = kineticEnergy + mass;
  G4double totalMomentum = sqrt(kineticEnergy*(kineticEnergy + 2.0*mass));

  G4double func1 = tmin*
    ComputeDMicroscopicCrossSection(kineticEnergy,Z,tmin);

  G4double lnepksi, epksi;
  G4double func2;

  G4double xmin = G4Log(tmin/MeV);
  G4double xmax = G4Log(kineticEnergy/tmin);

  do {
    lnepksi = xmin + G4UniformRand()*xmax;
    epksi   = MeV*G4Exp(lnepksi);
    func2   = epksi*ComputeDMicroscopicCrossSection(kineticEnergy,Z,epksi);

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while(func2 < func1*G4UniformRand());

  G4double gEnergy = epksi;

  // ===== sample angle =====

  G4double gam  = totalEnergy/mass;
  G4double rmax = gam*std::min(1.0, totalEnergy/gEnergy - 1.0);
  G4double rmax2= rmax*rmax;
  G4double x = G4UniformRand()*rmax2/(1.0 + rmax2);

  G4double theta = sqrt(x/(1.0 - x))/gam;
  G4double sint  = sin(theta);
  G4double phi   = twopi * G4UniformRand() ;
  G4double dirx  = sint*cos(phi), diry = sint*sin(phi), dirz = cos(theta) ;

  G4ThreeVector gDirection(dirx, diry, dirz);
  gDirection.rotateUz(partDirection);

  partDirection *= totalMomentum;
  partDirection -= gEnergy*gDirection;
  partDirection = partDirection.unit();

  // primary change
  kineticEnergy -= gEnergy;
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(partDirection);

  // save secondary
  G4DynamicParticle* aGamma = 
    new G4DynamicParticle(theGamma,gDirection,gEnergy);
  vdp->push_back(aGamma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
// $Id: G4hLDMBremModel.cc 74020 2013-09-19 13:38:38Z gcosmo $
//
// -------------------------------------------------------------------
//
// 21.03.17 V. Grichine based on G4hBremsstrahlungModel
//
// Class Description:
//
// Implementation of energy loss for LDMPhoton emission by hadrons
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4LDMBremModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4LDMPhoton.hh"
#include "G4ParticleChangeForLoss.hh"
#include "TestParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4LDMBremModel::G4LDMBremModel(const G4ParticleDefinition* p,
                               const G4String& nam)
  : G4MuBremsstrahlungModel(p, nam)
{
  fEpsilon = TestParameters::GetPointer()->GetAlphaFactor();
  theLDMPhoton = G4LDMPhoton::LDMPhoton();
  fLDMPhotonMass = theLDMPhoton->GetPDGMass();
  minThreshold = 1.2*fLDMPhotonMass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LDMBremModel::~G4LDMBremModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4LDMBremModel::ComputeDEDXPerVolume(const G4Material*,
                                              const G4ParticleDefinition*,
                                              G4double, G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4LDMBremModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double gammaEnergy)
//  differential cross section
{
  G4double dxsection = 0.;

  if( gammaEnergy > tkin || tkin < minThreshold) return dxsection;
  /*
  G4cout << "G4LDMBremModel m= " << mass 
         << "  " << particle->GetParticleName() 
         << "  Egamma(GeV)= " << gammaEnergy/GeV 
         << "  Ekin(GeV)= " << tkin/GeV << G4endl;
  */
  G4double E = tkin + mass ;
  G4double v = gammaEnergy/E ;
  G4double delta = 0.5*mass*mass*v/(E-gammaEnergy) ;
  G4double rab0=delta*sqrte ;

  G4int iz = std::max(1,std::min(G4lrint(Z),99));

  G4double z13 = 1.0/nist->GetZ13(iz);
  G4double dn  = mass*nist->GetA27(iz)/(70.*MeV);

  G4double    b = btf;
  if(1 == iz) b = bh;

  // nucleus contribution logarithm
  G4double rab1=b*z13;
  G4double fn=G4Log(rab1/(dn*(electron_mass_c2+rab0*rab1))*
              (mass+delta*(dn*sqrte-2.))) ;
  if(fn <0.) fn = 0. ;

  G4double x = 1.0 - v;

  if(particle->GetPDGSpin() != 0) { x += 0.75*v*v; }

  dxsection  = coeff*x*Z*Z*fn/gammaEnergy;
  return dxsection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4LDMBremModel::ComputeCrossSectionPerAtom(
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
  
  cut = std::max(cut, minThreshold);
  if (cut >= tmax) return cross;

  cross = ComputeMicroscopicCrossSection (kineticEnergy, Z, cut);

  if(tmax < kineticEnergy) 
  {
    cross -= ComputeMicroscopicCrossSection(kineticEnergy, Z, tmax);
  }
  cross *= fEpsilon*fEpsilon;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LDMBremModel::SampleSecondaries(
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
  tmin = std::max(tmin, minThreshold);
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
  G4double xmax = G4Log(tmax/tmin);

  do 
  {
    lnepksi = xmin + G4UniformRand()*xmax;
    epksi   = MeV*G4Exp(lnepksi);
    func2   = epksi*ComputeDMicroscopicCrossSection(kineticEnergy,Z,epksi);

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while(func2 < func1*G4UniformRand());

  G4double gEnergy = std::max(epksi, fLDMPhotonMass);
  G4double gMomentum = 
    std::sqrt((epksi - fLDMPhotonMass)*(epksi + fLDMPhotonMass));

  // ===== sample angle =====

  G4double gam  = totalEnergy/mass;
  G4double rmax = gam*std::min(1.0, totalEnergy/gEnergy - 1.0);
  G4double rmax2= rmax*rmax;
  G4double x = G4UniformRand()*rmax2/(1.0 + rmax2);

  G4double theta = std::sqrt(x/(1.0 - x))/gam;
  G4double sint  = std::sin(theta);
  G4double phi   = twopi * G4UniformRand() ;
  G4double dirx  = sint*cos(phi), diry = sint*sin(phi), dirz = cos(theta) ;

  G4ThreeVector gDirection(dirx, diry, dirz);
  gDirection.rotateUz(partDirection);

  partDirection *= totalMomentum;
  partDirection -= gMomentum*gDirection;
  partDirection = partDirection.unit();

  // primary change

  kineticEnergy -= gEnergy;

  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(partDirection);

  // save secondary
  G4DynamicParticle* aLDMPhoton = 
    new G4DynamicParticle(theLDMPhoton,gDirection,gEnergy);
  vdp->push_back(aLDMPhoton);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

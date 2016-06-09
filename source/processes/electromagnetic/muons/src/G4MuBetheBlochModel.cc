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
// $Id: G4MuBetheBlochModel.cc,v 1.14 2004/12/03 17:32:03 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4MuBetheBlochModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 09.08.2002
//
// Modifications:
//
// 04-12-02 Fix problem of G4DynamicParticle constructor (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 10-02-04 Calculation of radiative corrections using R.Kokoulin model (V.Ivanchenko)
//

//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MuBetheBlochModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

G4double G4MuBetheBlochModel::xgi[]={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801 };
G4double G4MuBetheBlochModel::wgi[]={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4MuBetheBlochModel::G4MuBetheBlochModel(const G4ParticleDefinition* p,
                                         const G4String& nam)
  : G4VEmModel(nam),
  particle(0),
  limitKinEnergy(100.*keV),
  logLimitKinEnergy(log(limitKinEnergy)),
  highKinEnergy(100.*TeV),
  lowKinEnergy(1.0*GeV),
  twoln10(2.0*log(10.0)),
  bg2lim(0.0169),
  taulim(8.4146e-3),
  alphaprime(fine_structure_const/twopi)
{
  if(p) SetParticle(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuBetheBlochModel::~G4MuBetheBlochModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuBetheBlochModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  massSquare = mass*mass;
  ratio = electron_mass_c2/mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::HighEnergyLimit(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::LowEnergyLimit(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::MinEnergyCut(const G4ParticleDefinition*,
                                           const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuBetheBlochModel::IsInCharge(const G4ParticleDefinition* p)
{
  if(!particle) SetParticle(p);
  return (p->GetPDGCharge() != 0.0 && p->GetPDGMass() > 10.*MeV
                                   && p->GetPDGSpin() == 0.5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuBetheBlochModel::Initialise(const G4ParticleDefinition* p,
                                     const G4DataVector&)
{
  if(!particle) SetParticle(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::ComputeDEDX(const G4MaterialCutsCouple* couple,
                                               const G4ParticleDefinition* p,
                                                     G4double kineticEnergy,
                                                     G4double cut)
{
  G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tau   = kineticEnergy/mass;
  G4double cutEnergy = min(cut,tmax);
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau+2.0);
  G4double beta2 = bg2/(gam*gam);

  const G4Material* material = couple->GetMaterial();
  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc2 = eexc*eexc;
  G4double taul  = material->GetIonisation()->GetTaul();
  G4double cden  = material->GetIonisation()->GetCdensity();
  G4double mden  = material->GetIonisation()->GetMdensity();
  G4double aden  = material->GetIonisation()->GetAdensity();
  G4double x0den = material->GetIonisation()->GetX0density();
  G4double x1den = material->GetIonisation()->GetX1density();
  G4double* shellCorrectionVector =
            material->GetIonisation()->GetShellCorrectionVector();
  G4double eDensity = material->GetElectronDensity();

  G4double dedx = log(2.0*electron_mass_c2*bg2*cutEnergy/eexc2)-(1.0 + cutEnergy/tmax)*beta2;

  G4double totEnergy = kineticEnergy + mass;
  G4double del = 0.5*cutEnergy/totEnergy;
  dedx += del*del;

  // density correction
  G4double x = log(bg2)/twoln10;
  if ( x >= x0den ) {
    dedx -= twoln10*x - cden ;
    if ( x < x1den ) dedx -= aden*pow((x1den-x),mden) ;
  }

  // shell correction
  G4double sh = 0.0;
  x  = 1.0;

  if ( bg2 > bg2lim ) {
    for (G4int k=0; k<3; k++) {
	x *= bg2 ;
	sh += shellCorrectionVector[k]/x;
    }

  } else {
    for (G4int k=0; k<3; k++) {
	x *= bg2lim ;
	sh += shellCorrectionVector[k]/x;
    }
    sh *= log(tau/taul)/log(taulim/taul);
  }
  dedx -= sh;

  // now compute the total ionization loss

  if (dedx < 0.0) dedx = 0.0 ;

  // radiative corrections of R. Kokoulin
  if (cutEnergy > limitKinEnergy) {

    G4double logtmax = log(cutEnergy);
    G4double logstep = logtmax - logLimitKinEnergy;
    G4double dloss = 0.0;
    G4double ftot2= 0.5/(totEnergy*totEnergy);

    for (G4int ll=0; ll<8; ll++)
    {
      G4double ep = exp(logLimitKinEnergy + xgi[ll]*logstep);
      G4double a1 = log(1.0 + 2.0*ep/electron_mass_c2);
      G4double a3 = log(4.0*totEnergy*(totEnergy - ep)/massSquare);
      dloss += wgi[ll]*(1.0 - beta2*ep/tmax + ep*ep*ftot2)*a1*(a3 - a1);
    }
    dedx += dloss*logstep*alphaprime;
  }

  dedx *= twopi_mc2_rcl2*eDensity/beta2;

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::CrossSection(const G4MaterialCutsCouple* couple,
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxKinEnergy)
{
  G4double cross = 0.0;
  G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  G4double maxEnergy = min(tmax,maxKinEnergy);
  if(cutEnergy < maxEnergy) {

    G4double totEnergy = kineticEnergy + mass;
    G4double energy2   = totEnergy*totEnergy;
    G4double beta2     = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;

    cross = 1.0/cutEnergy - 1.0/maxEnergy - beta2*log(maxEnergy/cutEnergy)/tmax
          + 0.5*(maxEnergy - cutEnergy)/energy2;

    // radiative corrections of R. Kokoulin
    if (maxEnergy > limitKinEnergy) {

      G4double logtmax = log(maxEnergy);
      G4double logtmin = log(max(cutEnergy,limitKinEnergy));
      G4double logstep = logtmax - logtmin;
      G4double dcross  = 0.0;

      for (G4int ll=0; ll<8; ll++)
      {
        G4double ep = exp(logtmin + xgi[ll]*logstep);
        G4double a1 = log(1.0 + 2.0*ep/electron_mass_c2);
        G4double a3 = log(4.0*totEnergy*(totEnergy - ep)/massSquare);
        dcross += wgi[ll]*(1.0/ep - beta2/tmax + 0.5*ep/energy2)*a1*(a3 - a1);
      }

      cross += dcross*logstep*alphaprime;
    }

    cross *= twopi_mc2_rcl2*(couple->GetMaterial()->GetElectronDensity())/beta2;

  }

  //  G4cout << "tmin= " << cutEnergy << " tmax= " << tmax
  //       << " cross= " << cross << G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* G4MuBetheBlochModel::SampleSecondary(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle* dp,
                                   G4double minEnergy,
                                   G4double maxEnergy)
{
  G4double tmax = MaxSecondaryEnergy(dp);
  G4double maxKinEnergy = min(maxEnergy,tmax);
  G4double minKinEnergy = min(minEnergy,maxKinEnergy);

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double totEnergy     = kineticEnergy + mass;
  G4double etot2         = totEnergy*totEnergy;
  G4double beta2         = kineticEnergy*(kineticEnergy + 2.0*mass)/etot2;
 
  G4double grej  = 1.;
  if(tmax > limitKinEnergy) {
    G4double a0    = log(2.*totEnergy/mass);
    grej  += alphaprime*a0*a0;
  }

  G4double deltaKinEnergy, f;

  // sampling follows ...
  do {
    G4double q = G4UniformRand();
    deltaKinEnergy = minKinEnergy*maxKinEnergy/(minKinEnergy*(1.0 - q) + maxKinEnergy*q);


    f = 1.0 - beta2*deltaKinEnergy/tmax + 0.5*deltaKinEnergy*deltaKinEnergy/etot2;

    if(deltaKinEnergy > limitKinEnergy) {
      G4double a1 = log(1.0 + 2.0*deltaKinEnergy/electron_mass_c2);
      G4double a3 = log(4.0*totEnergy*(totEnergy - deltaKinEnergy)/massSquare);
      f *= (1. + alphaprime*a1*(a3 - a1));
    }

    if(f > grej) {
        G4cout << "G4MuBetheBlochModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for edelta= " << deltaKinEnergy
               << " tmin= " << minKinEnergy << " max= " << maxKinEnergy
               << G4endl;
    }


  } while( grej*G4UniformRand() > f );

  G4double deltaMomentum =
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double totMomentum = totEnergy*sqrt(beta2);
  G4double cost = deltaKinEnergy * (totEnergy + electron_mass_c2) /
                                   (deltaMomentum * totMomentum);

  G4double sint = sqrt(1.0 - cost*cost);

  G4double phi = twopi * G4UniformRand() ;

  G4ThreeVector deltaDirection(sint*cos(phi),sint*sin(phi), cost) ;
  G4ThreeVector direction = dp->GetMomentumDirection();
  deltaDirection.rotateUz(direction);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = new G4DynamicParticle(G4Electron::Electron(),
                                                   deltaDirection,deltaKinEnergy);

  return delta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

vector<G4DynamicParticle*>* G4MuBetheBlochModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{
  vector<G4DynamicParticle*>* vdp = new vector<G4DynamicParticle*>;
  G4DynamicParticle* delta = SampleSecondary(couple,dp,tmin,maxEnergy);
  vdp->push_back(delta);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

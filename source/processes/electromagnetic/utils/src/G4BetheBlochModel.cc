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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BetheBlochModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 03.01.2002
//
// Modifications: 04.12.2002 VI Fix problem of G4DynamicParticle constructor
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4BetheBlochModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BetheBlochModel::G4BetheBlochModel(const G4ParticleDefinition* p) 
  : G4VEmModel(),
  particle(0),
  highKinEnergy(100.*TeV),
  lowKinEnergy(2.0*MeV),
  twoln10(2.0*log(10.0)),
  bg2lim(0.0169), 
  taulim(8.4146e-3)
{
  if(p) SetParticle(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BetheBlochModel::~G4BetheBlochModel() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BetheBlochModel::SetParticle(const G4ParticleDefinition* p) 
{
  particle = p;
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  G4double q = particle->GetPDGCharge()/eplus;
  chargeSquare = q*q;
  lowKinEnergy *= mass/proton_mass_c2;
  ratio = electron_mass_c2/mass;
  qc = mass/ratio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BetheBlochModel::HighEnergyLimit(const G4ParticleDefinition* p,
                                            const G4Material*) 
{
  if(!particle) SetParticle(p);
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4double G4BetheBlochModel::LowEnergyLimit(const G4ParticleDefinition* p,
                                           const G4Material*) 
{
  if(!particle) SetParticle(p);
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BetheBlochModel::MinEnergyCut(const G4ParticleDefinition* p,
                                         const G4Material* material) 
{
  return material->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4bool G4BetheBlochModel::IsInCharge(const G4ParticleDefinition* p,
			             const G4Material*) 
{
  return (p->GetPDGCharge() != 0.0 && p->GetPDGMass() > 10.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BetheBlochModel::ComputeDEDX(const G4Material* material,
                                        const G4ParticleDefinition* p,
                                              G4double kineticEnergy,
                                              G4double cutEnergy) 
{
  if(!particle) SetParticle(p);
  G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tau   = kineticEnergy/mass;
  G4double x     = 1.0;
  if(cutEnergy < tmax) x = cutEnergy/tmax;
  G4double gam   = tau + 1.0;    
  G4double beta2 = 1. - 1./(gam*gam);
  G4double bg2   = tau * (tau+2.0);

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
        
  G4double dedx = log(2.0*electron_mass_c2*bg2*tmax*x/eexc2)-(1.0 + x)*beta2;
    
  if(0.5 == spin) {
    G4double del = 0.5*x*tmax/(kineticEnergy + mass);
    dedx += del*del;
  }

  // density correction     
  x = log(bg2)/twoln10;
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
  
  dedx *= twopi_mc2_rcl2*chargeSquare*eDensity/beta2;

  return dedx; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BetheBlochModel::CrossSection(const G4Material* material,
                                         const G4ParticleDefinition* p,
                                               G4double kineticEnergy,
                                               G4double cutEnergy,
                                               G4double maxEnergy) 
{
  if(!particle) SetParticle(p);
  G4double cross = 0.0;
  G4double tmax = G4std::min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);
  if(cutEnergy < tmax) {
    
    G4double x      = cutEnergy/tmax;
    G4double energy = kineticEnergy + mass;
    G4double gam    = energy/mass;
    G4double beta2  = 1. - 1./(gam*gam);
    cross = (1.0 - x*(1.0 - beta2*log(x)))/cutEnergy;
    
    // +term for spin=1/2 particle
    if( 0.5 == spin ) {
      cross +=  0.5 * (tmax - cutEnergy) / (energy*energy);
    
    // +term for spin=1 particle
    } else if( 0.9 < spin ) {
      
      cross += -log(x)/(3.0*qc) +
	(tmax - cutEnergy) * ((1.0+ 0.25*tmax*(1.0 + x)/qc)/(energy*energy)
	- beta2 / (tmax * qc) )/3.0;
    }
    
    cross *= twopi_mc2_rcl2*chargeSquare*material->GetElectronDensity()/beta2;
  }
  //  G4cout << "tmin= " << cutEnergy << " tmax= " << tmax 
  //       << " cross= " << cross << G4endl; 
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4DynamicParticle*>* G4BetheBlochModel::SampleSecondary(
                             const G4Material* material,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy) 
{
  G4double tmax = MaxSecondaryEnergy(dp);
  G4double xmin = tmin/tmax;
  G4double xmax = G4std::min(tmax, maxEnergy)/tmax;
  if(xmin >= xmax) return 0;

  G4ThreeVector momentum = dp->GetMomentumDirection();

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double energy = kineticEnergy + mass;
  G4double beta2 = 1. - mass*mass/(energy*energy);
  G4double totMomentum = sqrt(energy*energy - mass*mass);
  G4double x = 0.0;
  G4double y = 0.0;
  G4double grej = 1.0 - beta2*xmin;
  if(spin > 0.0) {
    x = 0.5*tmax*tmax/(energy*energy);
    if(spin < 0.9) {
      grej += xmax*xmax*x;
    } else {
      y = tmax/(3.0*qc);
      grej += (1.0 - beta2*xmin)*xmin*y + x*xmin*xmin*(2.0/3.0 + xmin*y);
    }
  }
  G4double z, f;

  // sampling follows ...      
  do {
    G4double q = G4UniformRand();
    z = xmin*xmax/(xmin*(1.0 - q) + xmax*q);

    f = 1.0 - z * (beta2 - z*x);

    if (spin > 0.9) {
      f += (1.0 - beta2*z)*z*y + x*z*z*(2.0/3.0 + z*y);
    }
    if(f > grej) {
        G4cout << "G4BetheBlochModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for x= " << z
               << G4endl; 
    }

  } while( grej*G4UniformRand() > f );
  
  G4double deltaKinEnergy = z * tmax;
    
  G4double deltaMomentum = 
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
                                   (deltaMomentum * totMomentum);
  G4double sint = sqrt(1.0 - cost*cost);
 
  G4double phi = twopi * G4UniformRand() ; 


  G4ThreeVector deltaDirection(sint*cos(phi),sint*sin(phi), cost) ;
  deltaDirection.rotateUz(momentum);
  
  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = new G4DynamicParticle();
  delta->SetDefinition(G4Electron::Electron());
  delta->SetKineticEnergy(deltaKinEnergy);
  delta->SetMomentumDirection(deltaDirection);
  G4std::vector<G4DynamicParticle*>* vdp = new G4std::vector<G4DynamicParticle*>;
  vdp->push_back(delta);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

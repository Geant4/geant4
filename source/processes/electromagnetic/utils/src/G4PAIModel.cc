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
// File name:     G4PAIModel.cc
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 03.01.2002
//
// Modifications:
//


#include "G4PAIModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"


////////////////////////////////////////////////////////////////////////

G4PAIModel::G4PAIModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),G4VEmFluctuationModel(nam),
  fParticle(0),
  fHighKinEnergy(100.*TeV),
  fLowKinEnergy(2.0*MeV),
  fTwoln10(2.0*log(10.0)),
  fBg2lim(0.0169),
  fTaulim(8.4146e-3)
{
  if(p) SetParticle(p);
}

////////////////////////////////////////////////////////////////////////////

G4PAIModel::~G4PAIModel()
{}

///////////////////////////////////////////////////////////////////////////////

void G4PAIModel::SetParticle(const G4ParticleDefinition* p)
{
  fParticle = p;
  fMass = fParticle->GetPDGMass();
  fSpin = fParticle->GetPDGSpin();
  G4double q = fParticle->GetPDGCharge()/eplus;
  fChargeSquare = q*q;
  fLowKinEnergy *= fMass/proton_mass_c2;
  fRatio = electron_mass_c2/fMass;
  fQc = fMass/fRatio;
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::HighEnergyLimit(const G4ParticleDefinition* p)
{
  if(!fParticle) SetParticle(p);
  return fHighKinEnergy;
}

///////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::LowEnergyLimit(const G4ParticleDefinition* p)
{
  if(!fParticle) SetParticle(p);
  return fLowKinEnergy;
}

////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::MinEnergyCut(const G4ParticleDefinition*,
                                         const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

////////////////////////////////////////////////////////////////////////////

G4bool G4PAIModel::IsInCharge(const G4ParticleDefinition* p)
{
  if(!fParticle) SetParticle(p);
  return (p->GetPDGCharge() != 0.0 && p->GetPDGMass() > 10.*MeV);
}

////////////////////////////////////////////////////////////////////////////

void G4PAIModel::Initialise(const G4ParticleDefinition* p,
                                   const G4DataVector&)
{
  if(!fParticle) SetParticle(p);
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::ComputeDEDX(const G4Material* material,
                                        const G4ParticleDefinition* p,
                                              G4double kineticEnergy,
                                              G4double cutEnergy)
{
  G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tau   = kineticEnergy/fMass;
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
    
  if(0.5 == fSpin) 
  {
    G4double del = 0.5*x*tmax/(kineticEnergy + fMass);
    dedx += del*del;
  }

  // density correction     
  x = log(bg2)/fTwoln10;
  if ( x >= x0den ) 
  {
    dedx -= fTwoln10*x - cden ;
    if ( x < x1den ) dedx -= aden*pow((x1den-x),mden) ;
  }
    
  // shell correction 
  G4double sh = 0.0;      
  x  = 1.0;

  if ( bg2 > fBg2lim ) 
  {
    for (G4int k=0; k<3; k++) 
    {
	x *= bg2 ;
	sh += shellCorrectionVector[k]/x;
    }

  } 
  else 
  {
    for (G4int k=0; k<3; k++) 
    {
	x *= fBg2lim ;
	sh += shellCorrectionVector[k]/x;
    }
    sh *= log(tau/taul)/log(fTaulim/taul);
  }
  dedx -= sh;
    
  // now compute the total ionization loss

  if (dedx < 0.0) dedx = 0.0 ;
  
  dedx *= twopi_mc2_rcl2*fChargeSquare*eDensity/beta2;

  return dedx; 
}

/////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::CrossSection(const G4Material* material,
                                         const G4ParticleDefinition* p,
                                               G4double kineticEnergy,
                                               G4double cutEnergy,
                                               G4double maxEnergy) 
{
  G4double cross = 0.0;
  G4double tmax = std::min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);
  if(cutEnergy < tmax) 
  {

    G4double x      = cutEnergy/tmax;
    G4double energy = kineticEnergy + fMass;
    G4double gam    = energy/fMass;
    G4double beta2  = 1. - 1./(gam*gam);
    cross = (1.0 - x*(1.0 - beta2*log(x)))/cutEnergy;

    // +term for spin=1/2 particle
    if( 0.5 == fSpin ) 
    {
      cross +=  0.5 * (tmax - cutEnergy) / (energy*energy);

    // +term for spin=1 particle
    } 
    else if( 0.9 < fSpin ) 
    {

      cross += -log(x)/(3.0*fQc) +
	(tmax - cutEnergy) * ((1.0+ 0.25*tmax*(1.0 + x)/fQc)/(energy*energy)
	- beta2 / (tmax * fQc) )/3.0;
    }
    
    cross *= twopi_mc2_rcl2*fChargeSquare*material->GetElectronDensity()/beta2;
  }
  //  G4cout << "tmin= " << cutEnergy << " tmax= " << tmax 
  //       << " cross= " << cross << G4endl; 
  return cross;
}

/////////////////////////////////////////////////////////////////////////////

G4DynamicParticle* G4PAIModel::SampleSecondary(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{
  G4double tmax = MaxSecondaryEnergy(dp);
  G4double xmin = tmin/tmax;
  G4double xmax = std::min(tmax, maxEnergy)/tmax;
  if(xmin >= xmax) return 0;

  G4ThreeVector momentum = dp->GetMomentumDirection();

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double energy = kineticEnergy + fMass;
  G4double beta2 = 1. - fMass*fMass/(energy*energy);
  G4double totMomentum = sqrt(energy*energy - fMass*fMass);
  G4double x = 0.0;
  G4double y = 0.0;
  G4double grej = 1.0 - beta2*xmin;
  if(fSpin > 0.0) 
  {
    x = 0.5*tmax*tmax/(energy*energy);
    if(fSpin < 0.9) 
    {
      grej += xmax*xmax*x;
    } 
    else 
    {
      y = tmax/(3.0*fQc);
      grej += (1.0 - beta2*xmin)*xmin*y + x*xmin*xmin*(2.0/3.0 + xmin*y);
    }
  }
  G4double z, f;

  // sampling follows ...
  do 
  {
    G4double q = G4UniformRand();
    z = xmin*xmax/(xmin*(1.0 - q) + xmax*q);

    f = 1.0 - z * (beta2 - z*x);

    if (fSpin > 0.9) 
    {
      f += (1.0 - beta2*z)*z*y + x*z*z*(2.0/3.0 + z*y);
    }
    if(f > grej) 
    {
        G4cout << "G4PAIModel::SampleSecondary Warning! "
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

  return delta;
}

////////////////////////////////////////////////////////////////////////////

std::vector<G4DynamicParticle*>* G4PAIModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{
  std::vector<G4DynamicParticle*>* vdp = new std::vector<G4DynamicParticle*>;
  G4DynamicParticle* delta = SampleSecondary(couple, dp, tmin, maxEnergy);
  vdp->push_back(delta);

  return vdp;
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIModel::SampleFluctuations(const G4Material* material,
                                                const G4DynamicParticle* dp,
				                      G4double& tmax,
					              G4double& length,
                                                      G4double& meanLoss)
{
  G4double loss=0.,minLoss=0.;

  if(meanLoss < minLoss) return meanLoss;

  return loss;
}

////////////////////////////////////////////////////////////////////////..


G4double G4PAIModel::Dispersion( const G4Material* material, 
                                 const G4DynamicParticle* dp,
 				       G4double& tmax, 
			               G4double& length       )
{
  G4double electronDensity = material->GetElectronDensity();
  G4double particleMass = dp->GetDefinition()->GetPDGMass();
  G4double gam   = (dp->GetKineticEnergy())/particleMass + 1.0; 
  G4double beta2 = 1.0 - 1.0/(gam*gam);
 
  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                 * electronDensity * fChargeSquare;

  return siga;
}


//
//
/////////////////////////////////////////////////







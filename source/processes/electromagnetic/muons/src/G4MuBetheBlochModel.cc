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
// File name:     G4MuBetheBlochModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 09.08.2002
//
// Modifications: 
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MuBetheBlochModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuBetheBlochModel::G4MuBetheBlochModel(const G4ParticleDefinition* p) 
  : G4VEmModel(),
  particle(0),
  mass(proton_mass_c2),
  chargeSquare(1.0),
  highKinEnergy(100.*TeV),
  lowKinEnergy(2.0*MeV),
  twoln10(2.0*log(10.0)),
  bg2lim(0.0169), 
  taulim(8.4146e-3)
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
  G4double q = particle->GetPDGCharge()/eplus;
  chargeSquare = q*q;
  lowKinEnergy *= mass/proton_mass_c2;
  ratio = electron_mass_c2/mass;
  qc = mass/ratio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::HighEnergyLimit(const G4ParticleDefinition* p,
                                            const G4Material*) 
{
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4double G4MuBetheBlochModel::LowEnergyLimit(const G4ParticleDefinition* p,
                                           const G4Material*) 
{
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::MinEnergyCut(const G4ParticleDefinition* p,
                                         const G4Material* material) 
{
  return material->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4bool G4MuBetheBlochModel::IsInCharge(const G4ParticleDefinition* p,
			             const G4Material*) 
{
  return (p->GetPDGCharge() != 0.0 && p->GetPDGMass() > 10.*MeV 
                                   && p->GetPDGSpin() == 0.5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::ComputeDEDX(const G4Material* material,
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

  G4double totEnergy = kineticEnergy + mass;
  G4double del = 0.5*x*tmax/totEnergy;
  dedx += del*del;
    
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

     // correction of R. Kokoulin  // has been taken out *************** 
     //  G4double apar = log(2.*(tmax/electron_mass_c2 + 1.0));
     //  dedx += fine_structure_const*(log(2.*totEnergy/mass)-apar/3.)*apar*apar/twopi; 

  
  dedx *= twopi_mc2_rcl2*chargeSquare*eDensity/beta2;

  return dedx; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBetheBlochModel::CrossSection(const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy) 
{
  if(!particle) SetParticle(p);
  G4double cross = 0.0;
  G4double tmaxSecondary = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tmax = G4std::min(tmaxSecondary, maxEnergy);
  if(cutEnergy < tmax) {
    
    const G4ElementVector* theElementVector = material->GetElementVector();
    const G4double* NbOfAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
    G4int NumberOfElements = material->GetNumberOfElements();
 
    for (G4int iel=0; iel<NumberOfElements; iel++ ) {

      cross +=  NbOfAtomsPerVolume[iel]* CrossSectionPerAtom( 
                                    (*theElementVector)[iel]->GetZ(), 
                                    kineticEnergy, cutEnergy, tmax, tmaxSecondary);

    }

  }
  //  G4cout << "tmin= " << cutEnergy << " tmax= " << tmax 
  //       << " cross= " << cross << G4endl; 
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBetheBlochModel::CrossSectionPerAtom(
                                 G4double Z,
                                 G4double kineticEnergy,
                                 G4double tmin,
                                 G4double tmax,
                                 G4double tmaxSecondary)
{
     
  G4double cross = 0.;
    
  const G4double xgi[] = {0.06943,0.33001,0.66999,0.93057};
  const G4double wgi[] = {0.17393,0.32607,0.32607,0.17393};
  const G4double ak1 = 4.6;
  const G4int k2 = 2;

  G4double aaa = log(tmin);
  G4double bbb = log(tmax);
  G4int    kkk = int((bbb-aaa)/ak1)+k2;
  G4double hhh = (bbb-aaa)/(float)kkk;
  G4double step = exp(hhh);
  G4double ymax = 1./tmax;

     
  for (G4int k=0; k<kkk; k++) {

    G4double ymin = ymax;
    ymax = ymin*step;
    G4double hhy = ymax-ymin;

    for (G4int i=0; i<4; i++) {

      G4double y = ymin+hhy*xgi[i];
      G4double ep = 1./y ;
      cross += ep*ep*wgi[i]*hhy*DifCrossSectionPerAtom(kineticEnergy, ep, tmaxSecondary);
    }
  }
  cross *= twopi_mc2_rcl2*Z*chargeSquare;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBetheBlochModel::DifCrossSectionPerAtom(G4double kineticEnergy,
                                                     G4double knockonEnergy,
                                                     G4double tmaxSecondary)

 // Calculates the differential cross section per atom
 //   using the cross section formula of R.P. Kokoulin (10/98)
{
  const G4double alphaprime = fine_structure_const/twopi;
  G4double totalEnergy = kineticEnergy + mass;
  G4double gam         = totalEnergy/mass;
  G4double beta2       = 1. - 1./(gam*gam);

  G4double v = knockonEnergy/totalEnergy;
  G4double cross = (1.-beta2*knockonEnergy/tmaxSecondary + 0.5*v*v)/
                   (beta2*knockonEnergy*knockonEnergy);
  G4double a1 = log(1.+2.*knockonEnergy/electron_mass_c2);
  G4double a3 = log(4.*totalEnergy*(totalEnergy-knockonEnergy)/(mass*mass));
  cross  *= (1.+alphaprime*a1*(a3-a1)); 

  return cross;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4DynamicParticle*>* G4MuBetheBlochModel::SampleSecondary(
                             const G4Material* material,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy) 
{
  G4double tmax = MaxSecondaryEnergy(dp);
  G4double xmin = tmin/tmax;
  G4double xmax = G4std::min(tmax, maxEnergy)/tmax;
  if(xmin >= xmax) return 0;

  G4ThreeVector momentum = dp->GetMomentum();

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double totEnergy = kineticEnergy + mass;
  G4double beta2 = 1. - mass*mass/(totEnergy*totEnergy);
  G4double x = 0.5*tmax*tmax/(totEnergy*totEnergy);

  const G4double alphaprime = fine_structure_const/twopi; 
  G4double a0=log(2.*totEnergy/mass); 

  G4double grej = (1.0 - beta2*xmin + xmax*xmax*x)*(1. + alphaprime*a0*a0);
  
  G4double z, f, a1, twoep;
      
  // sampling follows ...      
  do {
    G4double q = G4UniformRand();
    z = xmin*xmax/(xmin*(1.0 - q) + xmax*q);

    twoep = 2.0*z*tmax;
    a1= log(1.0 + twoep/electron_mass_c2); 
    
    f = (1.0 - beta2 * z + z*z*x)
      * (1. + alphaprime*a1*(a0 + log((2.0*totEnergy-twoep)/mass) - a1));

    if(f > grej) {
        G4cout << "G4MuBetheBlochModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for x= " << z
               << G4endl; 
    }

  } while( grej*G4UniformRand() >= f );
  
  G4double deltaKinEnergy = z * tmax;
    
  G4double deltaMomentum = 
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double cost = deltaKinEnergy * (totEnergy + electron_mass_c2) /
                                   (deltaMomentum * momentum.mag());
  G4double sint = sqrt(1.0 - cost*cost);
 
  G4double phi = twopi * G4UniformRand() ; 

  G4ThreeVector deltaDirection(sint*cos(phi),sint*sin(phi), cost) ;
  deltaDirection.rotateUz(momentum);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = new G4DynamicParticle(G4Electron::Electron(),
                                                   deltaKinEnergy,
                                                   deltaDirection);

  G4std::vector<G4DynamicParticle*>* vdp = new G4std::vector<G4DynamicParticle*>;
  vdp->push_back(delta);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

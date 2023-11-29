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
// GEANT4 Class file
//
//
// File name:     G4eBremParametrizedModel
//
// Author:        Andreas Schaelicke 
//
// Creation date: 06.04.2011
//
// Modifications:
//
// Main References:
//  - based on G4eBremsstrahlungModel and G4eBremsstrahlungRelModel
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eBremParametrizedModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4LossTableManager.hh"
#include "G4ModifiedTsai.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4double G4eBremParametrizedModel::xgi[]={ 0.0199, 0.1017, 0.2372, 0.4083,
						  0.5917, 0.7628, 0.8983, 0.9801 };
const G4double G4eBremParametrizedModel::wgi[]={ 0.0506, 0.1112, 0.1569, 0.1813,
					    0.1813, 0.1569, 0.1112, 0.0506 };

static const G4double tlow = 1.*CLHEP::MeV;

//
// GEANT4 internal units.
//
static const G4double
     ah10 = 4.67733E+00, ah11 =-6.19012E-01, ah12 = 2.02225E-02,
     ah20 =-7.34101E+00, ah21 = 1.00462E+00, ah22 =-3.20985E-02,
     ah30 = 2.93119E+00, ah31 =-4.03761E-01, ah32 = 1.25153E-02;

static const G4double
     bh10 = 4.23071E+00, bh11 =-6.10995E-01, bh12 = 1.95531E-02,
     bh20 =-7.12527E+00, bh21 = 9.69160E-01, bh22 =-2.74255E-02,
     bh30 = 2.69925E+00, bh31 =-3.63283E-01, bh32 = 9.55316E-03;

static const G4double
     al00 =-2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
     al10 =-7.69748E-02, al11 =-6.91499E-02, al12 = 2.22453E-03,
     al20 = 4.06463E-02, al21 =-1.01281E-02, al22 = 3.40919E-04;

static const G4double
     bl00 = 1.04133E+00, bl01 =-9.43291E-03, bl02 =-4.54758E-04,
     bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 =-1.30718E-03,
     bl20 =-1.59391E-02, bl21 = 7.27752E-03, bl22 =-1.94405E-04;

using namespace std;

G4eBremParametrizedModel::G4eBremParametrizedModel(const G4ParticleDefinition* p,
						   const G4String& nam)
  : G4VEmModel(nam),
    particle(nullptr),
    fMigdalConstant(classic_electr_radius*electron_Compton_length*electron_Compton_length*4.0*pi),
    bremFactor(fine_structure_const*classic_electr_radius*classic_electr_radius*16./3.),
    isInitialised(false),
    isElectron(true)
{
  theGamma = G4Gamma::Gamma();

  minThreshold = 0.1*keV;
  lowKinEnergy = 10.*MeV;
  SetLowEnergyLimit(lowKinEnergy);  

  nist = G4NistManager::Instance();  

  SetAngularDistribution(new G4ModifiedTsai());

  particleMass = kinEnergy = totalEnergy = currentZ = z13 = z23 = lnZ = Fel = Finel 
    = densityFactor = densityCorr = fMax = fCoulomb = 0.;

  InitialiseConstants();
  if(nullptr != p) { SetParticle(p); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremParametrizedModel::InitialiseConstants()
{
  facFel = G4Log(184.15);
  facFinel = G4Log(1194.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eBremParametrizedModel::~G4eBremParametrizedModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremParametrizedModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  particleMass = p->GetPDGMass();
  if(p == G4Electron::Electron()) { isElectron = true; }
  else                            { isElectron = false;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremParametrizedModel::MinEnergyCut(const G4ParticleDefinition*,
						const G4MaterialCutsCouple*)
{
  return minThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremParametrizedModel::SetupForMaterial(const G4ParticleDefinition*,
						const G4Material* mat, 
						G4double kineticEnergy)
{
  densityFactor = mat->GetElectronDensity()*fMigdalConstant;

  // calculate threshold for density effect
  kinEnergy   = kineticEnergy;
  totalEnergy = kineticEnergy + particleMass;
  densityCorr = densityFactor*totalEnergy*totalEnergy;    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremParametrizedModel::Initialise(const G4ParticleDefinition* p,
					  const G4DataVector& cuts)
{
  if(p) { SetParticle(p); }

  lowKinEnergy  = LowEnergyLimit();

  currentZ = 0.;

  if(IsMaster()) { InitialiseElementSelectors(p, cuts); }

  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForLoss();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremParametrizedModel::InitialiseLocal(const G4ParticleDefinition*,
					       G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremParametrizedModel::ComputeDEDXPerVolume(
					     const G4Material* material,
                                             const G4ParticleDefinition* p,
					     G4double kineticEnergy,
					     G4double cutEnergy)
{
  if(!particle) { SetParticle(p); }
  if(kineticEnergy < lowKinEnergy) { return 0.0; }
  G4double cut = std::min(cutEnergy, kineticEnergy);
  if(cut == 0.0) { return 0.0; }

  SetupForMaterial(particle, material,kineticEnergy);

  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();

  G4double dedx = 0.0;

  //  loop for elements in the material
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4VEmModel::SetCurrentElement((*theElementVector)[i]);
    SetCurrentElement((*theElementVector)[i]->GetZ());

    dedx += theAtomicNumDensityVector[i]*currentZ*currentZ*ComputeBremLoss(cut);
  }
  dedx *= bremFactor;

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremParametrizedModel::ComputeBremLoss(G4double cut)
{
  G4double loss = 0.0;

  // number of intervals and integration step 
  G4double vcut = cut/totalEnergy;
  G4int n = (G4int)(20*vcut) + 3;
  G4double delta = vcut/G4double(n);

  G4double e0 = 0.0;
  G4double xs; 

  // integration
  for(G4int l=0; l<n; l++) {

    for(G4int i=0; i<8; i++) {

      G4double eg = (e0 + xgi[i]*delta)*totalEnergy;

      xs = ComputeDXSectionPerAtom(eg);

      loss += wgi[i]*xs/(1.0 + densityCorr/(eg*eg));
    }
    e0 += delta;
  }

  loss *= delta*totalEnergy;

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremParametrizedModel::ComputeCrossSectionPerAtom(
                                              const G4ParticleDefinition* p,
					      G4double kineticEnergy, 
					      G4double Z,   G4double,
					      G4double cutEnergy, 
					      G4double maxEnergy)
{
  if(!particle) { SetParticle(p); }
  if(kineticEnergy < lowKinEnergy) { return 0.0; }
  G4double cut  = std::min(cutEnergy, kineticEnergy);
  G4double tmax = std::min(maxEnergy, kineticEnergy);

  if(cut >= tmax) { return 0.0; }

  SetCurrentElement(Z);

  G4double cross = ComputeXSectionPerAtom(cut);

  // allow partial integration
  if(tmax < kinEnergy) { cross -= ComputeXSectionPerAtom(tmax); }
  
  cross *= Z*Z*bremFactor;

  return cross;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremParametrizedModel::ComputeXSectionPerAtom(G4double cut)
{
  G4double cross = 0.0;

  // number of intervals and integration step 
  G4double vcut = G4Log(cut/totalEnergy);
  G4double vmax = G4Log(kinEnergy/totalEnergy);
  G4int n = (G4int)(0.45*(vmax - vcut)) + 4;
  //  n=1; //  integration test 
  G4double delta = (vmax - vcut)/G4double(n);

  G4double e0 = vcut;
  G4double xs; 

  // integration
  for(G4int l=0; l<n; l++) {

    for(G4int i=0; i<8; i++) {

      G4double eg = G4Exp(e0 + xgi[i]*delta)*totalEnergy;

      xs = ComputeDXSectionPerAtom(eg);

      cross += wgi[i]*xs/(1.0 + densityCorr/(eg*eg));
    }
    e0 += delta;
  }

  cross *= delta;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// compute the value of the screening function 3*PHI1 - PHI2

G4double G4eBremParametrizedModel::ScreenFunction1(G4double ScreenVariable)
{
  G4double screenVal;

  if (ScreenVariable > 1.)
    screenVal = 42.24 - 8.368*G4Log(ScreenVariable+0.952);
  else
    screenVal = 42.392 - ScreenVariable* (7.796 - 1.961*ScreenVariable);

  return screenVal;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// compute the value of the screening function 1.5*PHI1 - 0.5*PHI2

G4double G4eBremParametrizedModel::ScreenFunction2(G4double ScreenVariable)
{
  G4double screenVal;

  if (ScreenVariable > 1.)
    screenVal = 42.24 - 8.368*G4Log(ScreenVariable+0.952);
  else
    screenVal = 41.734 - ScreenVariable* (6.484 - 1.250*ScreenVariable);

  return screenVal;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Parametrized cross section
G4double G4eBremParametrizedModel::ComputeParametrizedDXSectionPerAtom(
                                   G4double kineticEnergy, 
				   G4double gammaEnergy, G4double Z) 
{
  SetCurrentElement(Z);
  G4double FZ = lnZ* (4.- 0.55*lnZ);
  G4double Z3 = z13; 
  G4double ZZ = z13*nist->GetZ13(G4lrint(Z)+1); 

  totalEnergy = kineticEnergy + electron_mass_c2;

  //  G4double x, epsil, greject, migdal, grejmax, q;
  G4double epsil, greject;
  G4double U  = G4Log(kineticEnergy/electron_mass_c2);
  G4double U2 = U*U;

  // precalculated parameters
  G4double ah, bh;

  if (kineticEnergy > tlow) {
       
    G4double ah1 = ah10 + ZZ* (ah11 + ZZ* ah12);
    G4double ah2 = ah20 + ZZ* (ah21 + ZZ* ah22);
    G4double ah3 = ah30 + ZZ* (ah31 + ZZ* ah32);

    G4double bh1 = bh10 + ZZ* (bh11 + ZZ* bh12);
    G4double bh2 = bh20 + ZZ* (bh21 + ZZ* bh22);
    G4double bh3 = bh30 + ZZ* (bh31 + ZZ* bh32);

    ah = 1.   + (ah1*U2 + ah2*U + ah3) / (U2*U);
    bh = 0.75 + (bh1*U2 + bh2*U + bh3) / (U2*U);

    // limit of the screening variable
    G4double screenfac =
      136.*electron_mass_c2/(Z3*totalEnergy);

    epsil = gammaEnergy/totalEnergy; //         epsil = x*kineticEnergy/totalEnergy;
        G4double screenvar = screenfac*epsil/(1.0-epsil);
        G4double F1 = max(ScreenFunction1(screenvar) - FZ ,0.);
        G4double F2 = max(ScreenFunction2(screenvar) - FZ ,0.);


	greject = (F1 - epsil* (ah*F1 - bh*epsil*F2))/8.; //  1./(42.392 - FZ);

    std::cout << " yy = "<<epsil<<std::endl;
    std::cout << " F1/(...) "<<F1/(42.392 - FZ)<<std::endl;
    std::cout << " F2/(...) "<<F2/(42.392 - FZ)<<std::endl;
    std::cout << " (42.392 - FZ) " << (42.392 - FZ) <<std::endl;

  } else {  

    G4double al0 = al00 + ZZ* (al01 + ZZ* al02);
    G4double al1 = al10 + ZZ* (al11 + ZZ* al12);
    G4double al2 = al20 + ZZ* (al21 + ZZ* al22);
 
    G4double bl0 = bl00 + ZZ* (bl01 + ZZ* bl02);
    G4double bl1 = bl10 + ZZ* (bl11 + ZZ* bl12);
    G4double bl2 = bl20 + ZZ* (bl21 + ZZ* bl22);
 
    ah = al0 + al1*U + al2*U2;
    bh = bl0 + bl1*U + bl2*U2;

    G4double x=gammaEnergy/kineticEnergy;
    greject=(1. + x* (ah + bh*x));

    /*
    // Compute the maximum of the rejection function
    grejmax = max(1. + xmin* (ah + bh*xmin), 1.+ah+bh);
    G4double xm = -ah/(2.*bh);
    if ( xmin < xm && xm < xmax) grejmax = max(grejmax, 1.+ xm* (ah + bh*xm));
    */
  }

  return greject;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremParametrizedModel::ComputeDXSectionPerAtom(G4double gammaEnergy)
{

  if(gammaEnergy < 0.0) { return 0.0; }

  G4double y = gammaEnergy/totalEnergy;

  G4double main=0.;
  //secondTerm=0.;

  // ** form factors complete screening case **
  //  only valid for high energies (and if LPM suppression does not play a role)
  main   = (3./4.*y*y - y + 1.) * ( (Fel-fCoulomb) + Finel/currentZ );
  //  secondTerm = (1.-y)/12.*(1.+1./currentZ);

  std::cout<<" F1(0) "<<ScreenFunction1(0.) <<std::endl;
  std::cout<<" F1(0) "<<ScreenFunction2(0.) <<std::endl;
  std::cout<<"Ekin = "<<kinEnergy<<std::endl;
  std::cout<<"Z = "<<currentZ<<std::endl;
  std::cout<<"main  = "<<main<<std::endl;
  std::cout<<" y = "<<y<<std::endl;
  std::cout<<" Fel-fCoulomb "<< (Fel-fCoulomb) <<std::endl;

  G4double main2 = ComputeParametrizedDXSectionPerAtom(kinEnergy,gammaEnergy,currentZ);
  std::cout<<"main2 = "<<main2<<std::endl;
  std::cout<<"main2tot = "<<main2 * ( (Fel-fCoulomb) + Finel/currentZ )/(Fel-fCoulomb);

  G4double cross =  main2; //main+secondTerm;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremParametrizedModel::SampleSecondaries(
				      std::vector<G4DynamicParticle*>* vdp, 
				      const G4MaterialCutsCouple* couple,
				      const G4DynamicParticle* dp,
				      G4double cutEnergy,
				      G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  if(kineticEnergy < lowKinEnergy) { return; }
  G4double cut  = std::min(cutEnergy, kineticEnergy);
  G4double emax = std::min(maxEnergy, kineticEnergy);
  if(cut >= emax) { return; }

  SetupForMaterial(particle, couple->GetMaterial(),kineticEnergy);

  const G4Element* elm = SelectTargetAtom(couple,particle,kineticEnergy,
                                          dp->GetLogKineticEnergy(),cut,emax);
  SetCurrentElement(elm->GetZ());

  kinEnergy   = kineticEnergy;
  totalEnergy = kineticEnergy + particleMass;
  densityCorr = densityFactor*totalEnergy*totalEnergy;
 
  G4double xmin = G4Log(cut*cut + densityCorr);
  G4double xmax = G4Log(emax*emax  + densityCorr);
  G4double gammaEnergy, f, x; 

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  do {
    x = G4Exp(xmin + rndmEngine->flat()*(xmax - xmin)) - densityCorr;
    if(x < 0.0) x = 0.0;
    gammaEnergy = sqrt(x);
    f = ComputeDXSectionPerAtom(gammaEnergy);

    if ( f > fMax ) {
      G4cout << "### G4eBremParametrizedModel Warning: Majoranta exceeded! "
	     << f << " > " << fMax
	     << " Egamma(MeV)= " << gammaEnergy
	     << " E(mEV)= " << kineticEnergy
	     << G4endl;
    }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (f < fMax*rndmEngine->flat());

  //
  // angles of the emitted gamma. ( Z - axis along the parent particle)
  // use general interface
  //
  G4ThreeVector gammaDirection = 
    GetAngularDistribution()->SampleDirection(dp, totalEnergy-gammaEnergy,
					      G4lrint(currentZ), 
					      couple->GetMaterial());

  // create G4DynamicParticle object for the Gamma
  auto gamma = new G4DynamicParticle(theGamma,gammaDirection, gammaEnergy);
  vdp->push_back(gamma);
  
  G4double totMomentum = sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));
  G4ThreeVector direction = (totMomentum*dp->GetMomentumDirection()
			     - gammaEnergy*gammaDirection).unit();

  // energy of primary
  G4double finalE = kineticEnergy - gammaEnergy;

  // stop tracking and create new secondary instead of primary
  if(gammaEnergy > SecondaryThreshold()) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    auto el =
      new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
			    direction, finalE);
    vdp->push_back(el);

    // continue tracking
  } else {
    fParticleChange->SetProposedMomentumDirection(direction);
    fParticleChange->SetProposedKineticEnergy(finalE);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



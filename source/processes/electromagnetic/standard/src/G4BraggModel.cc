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
// File name:   G4BraggModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications: 
//
// 04-12-02 Fix problem of G4DynamicParticle constructor (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 04-06-03 Fix compilation warnings (V.Ivanchenko)
// 12-09-04 Add lowestKinEnergy and change order of if in DEDX method (VI)
// 11-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 16-06-05 Fix problem of chemical formula (V.Ivantchenko)
// 15-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 25-04-06 Add stopping data from PSTAR (V.Ivanchenko)
// 12-08-08 Added methods GetParticleCharge, GetChargeSquareRatio, 
//          CorrectionsAlongStep needed for ions(V.Ivanchenko)

// Class Description:
//
// Implementation of energy loss and delta-electron production by
// slow charged heavy particles

// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4BraggModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4LossTableManager.hh"
#include "G4EmCorrections.hh"
#include "G4EmParameters.hh"
#include "G4DeltaAngle.hh"
#include "G4ICRU90StoppingData.hh"
#include "G4NistManager.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PSTARStopping* G4BraggModel::fPSTAR = nullptr;

G4BraggModel::G4BraggModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),
    protonMassAMU(1.007276)
{
  SetHighEnergyLimit(2.0*CLHEP::MeV);

  lowestKinEnergy  = 0.25*CLHEP::keV;
  theZieglerFactor = CLHEP::eV*CLHEP::cm2*1.0e-15;
  theElectron = G4Electron::Electron();
  expStopPower125 = 0.0;

  corr = G4LossTableManager::Instance()->EmCorrections();
  if(nullptr != p) { SetParticle(p); }
  else  { SetParticle(theElectron); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BraggModel::~G4BraggModel()
{
  if(IsMaster()) { 
    delete fPSTAR; 
    fPSTAR = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggModel::Initialise(const G4ParticleDefinition* p,
                              const G4DataVector&)
{
  if(p != particle) { SetParticle(p); }

  // always false before the run
  SetDeexcitationFlag(false);

  if(IsMaster()) {
    if(nullptr == fPSTAR)  { fPSTAR = new G4PSTARStopping(); }
    if(particle->GetPDGMass() < CLHEP::GeV) { fPSTAR->Initialise(); }
    if(G4EmParameters::Instance()->UseICRU90Data()) {
      if(!fICRU90) { 
	fICRU90 = G4NistManager::Instance()->GetICRU90StoppingData(); 
      } else if(particle->GetPDGMass() < CLHEP::GeV) { fICRU90->Initialise(); }
    }
  }

  if(nullptr == fParticleChange) {

    if(UseAngularGeneratorFlag() && !GetAngularDistribution()) {
      SetAngularDistribution(new G4DeltaAngle());
    }
    G4String pname = particle->GetParticleName();
    if(particle->GetParticleType() == "nucleus" && 
       pname != "deuteron" && pname != "triton" &&
       pname != "alpha+"   && pname != "helium" &&
       pname != "hydrogen") { isIon = true; }

    fParticleChange = GetParticleChangeForLoss();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::GetChargeSquareRatio(const G4ParticleDefinition* p,
                                            const G4Material* mat,
                                            G4double kineticEnergy)
{
  // this method is called only for ions
  G4double q2 = corr->EffectiveChargeSquareRatio(p,mat,kineticEnergy);
  GetModelOfFluctuations()->SetParticleAndCharge(p, q2);
  return q2*corr->EffectiveChargeCorrection(p,mat,kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::MinEnergyCut(const G4ParticleDefinition*,
                                    const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::GetParticleCharge(const G4ParticleDefinition* p,
                                         const G4Material* mat,
                                         G4double kineticEnergy)
{
  // this method is called only for ions, so no check if it is an ion 
  return corr->GetParticleCharge(p,mat,kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::ComputeCrossSectionPerElectron(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cut,
                                                 G4double maxKinEnergy)
{
  G4double cross = 0.0;
  const G4double tmax      = MaxSecondaryEnergy(p, kineticEnergy);
  const G4double maxEnergy = std::min(tmax, maxKinEnergy);
  const G4double cutEnergy = std::max(cut, lowestKinEnergy*massRate);
  if(cutEnergy < maxEnergy) {

    const G4double energy  = kineticEnergy + mass;
    const G4double energy2 = energy*energy;
    const G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;
    cross = (maxEnergy - cutEnergy)/(cutEnergy*maxEnergy) 
      - beta2*G4Log(maxEnergy/cutEnergy)/tmax;

    if( 0.0 < spin ) { cross += 0.5*(maxEnergy - cutEnergy)/energy2; }

    cross *= CLHEP::twopi_mc2_rcl2*chargeSquare/beta2;
  }
  //   G4cout << "BR: e= " << kineticEnergy << " tmin= " << cutEnergy 
  //          << " tmax= " << tmax << " cross= " << cross << G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4BraggModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition* p,
                                         G4double kineticEnergy,
                                         G4double Z, G4double,
                                         G4double cutEnergy,
                                         G4double maxEnergy)
{
  return 
    Z*ComputeCrossSectionPerElectron(p,kineticEnergy,cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::CrossSectionPerVolume(const G4Material* material,
                                             const G4ParticleDefinition* p,
                                             G4double kineticEnergy,
                                             G4double cutEnergy,
                                             G4double maxEnergy)
{
  return material->GetElectronDensity()
    *ComputeCrossSectionPerElectron(p,kineticEnergy,cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::ComputeDEDXPerVolume(const G4Material* material,
                                            const G4ParticleDefinition* p,
                                            G4double kineticEnergy,
                                            G4double cut)
{
  const G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  const G4double tkin  = kineticEnergy/massRate;
  const G4double cutEnergy = std::max(cut, lowestKinEnergy*massRate);
  G4double dedx  = 0.0;

  if(tkin < lowestKinEnergy) {
    dedx = DEDX(material, lowestKinEnergy)*std::sqrt(tkin/lowestKinEnergy);
  } else {
    dedx = DEDX(material, tkin); 

    if (cutEnergy < tmax) {
      const G4double tau = kineticEnergy/mass;
      const G4double x   = cutEnergy/tmax;

      dedx += (G4Log(x)*(tau + 1.)*(tau + 1.)/(tau * (tau + 2.0)) + 1.0 - x) * 
	CLHEP::twopi_mc2_rcl2 * material->GetElectronDensity();
    }
  }
  dedx = std::max(dedx, 0.0) * chargeSquare;

  //G4cout << "E(MeV)= " << tkin/MeV << " dedx= " << dedx 
  //         << "  " << material->GetName() << G4endl;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
                                     const G4MaterialCutsCouple* couple,
                                     const G4DynamicParticle* dp,
                                     G4double minEnergy,
                                     G4double maxEnergy)
{
  const G4double tmax = MaxSecondaryKinEnergy(dp);
  const G4double xmax = std::min(tmax, maxEnergy);
  const G4double xmin  = std::max(lowestKinEnergy*massRate, minEnergy);
  if(xmin >= xmax) { return; }

  G4double kineticEnergy = dp->GetKineticEnergy();
  const G4double energy  = kineticEnergy + mass;
  const G4double energy2 = energy*energy;
  const G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;
  const G4double grej = 1.0;
  G4double deltaKinEnergy, f;

  CLHEP::HepRandomEngine* rndmEngineMod = G4Random::getTheEngine();
  G4double rndm[2];

  // sampling follows ...
  do {
    rndmEngineMod->flatArray(2, rndm);
    deltaKinEnergy = xmin*xmax/(xmin*(1.0 - rndm[0]) + xmax*rndm[0]);

    f = 1.0 - beta2*deltaKinEnergy/tmax;

    if(f > grej) {
        G4cout << "G4BraggModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for e= " << deltaKinEnergy
               << G4endl;
    }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while( grej*rndm[1] >= f );

  G4ThreeVector deltaDirection;

  if(UseAngularGeneratorFlag()) {
    const G4Material* mat =  couple->GetMaterial();
    G4int Z = SelectRandomAtomNumber(mat);

    deltaDirection = 
      GetAngularDistribution()->SampleDirection(dp, deltaKinEnergy, Z, mat);

  } else {
 
    G4double deltaMomentum =
      std::sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
    G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
      (deltaMomentum * dp->GetTotalMomentum());
    if(cost > 1.0) { cost = 1.0; }
    G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));

    G4double phi = twopi*rndmEngineMod->flat();

    deltaDirection.set(sint*std::cos(phi),sint*std::sin(phi), cost) ;
    deltaDirection.rotateUz(dp->GetMomentumDirection());
  }  

  // create G4DynamicParticle object for delta ray
  auto delta = new G4DynamicParticle(theElectron,deltaDirection,deltaKinEnergy);

  // Change kinematics of primary particle
  kineticEnergy -= deltaKinEnergy;
  G4ThreeVector finalP = dp->GetMomentum() - delta->GetMomentum();
  finalP               = finalP.unit();
  
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(finalP);

  vdp->push_back(delta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::MaxSecondaryEnergy(const G4ParticleDefinition* pd,
                                          G4double kinEnergy)
{
  if(pd != particle) { SetParticle(pd); }
  G4double tau  = kinEnergy/mass;
  G4double tmax = 2.0*electron_mass_c2*tau*(tau + 2.) /
                  (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggModel::HasMaterial(const G4Material* mat)
{
  const G4String& chFormula = mat->GetChemicalFormula();
  if(chFormula.empty()) { return; }

  // ICRU Report N49, 1993. Power's model for H
  static const G4int numberOfMolecula = 11;
  static const G4String molName[numberOfMolecula] = {
    "Al_2O_3",                 "CO_2",                      "CH_4",
    "(C_2H_4)_N-Polyethylene", "(C_2H_4)_N-Polypropylene",  "(C_8H_8)_N",
    "C_3H_8",                  "SiO_2",                     "H_2O",
    "H_2O-Gas",                "Graphite" } ;

  // Search for the material in the table
  for (G4int i=0; i<numberOfMolecula; ++i) {
    if (chFormula == molName[i]) {
      iMolecula = i;  
      return;
    }
  }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::StoppingPower(const G4Material* material,
                                     G4double kineticEnergy) 
{
  G4double ionloss = 0.0 ;

  if (iMolecula >= 0) {
  
    // The data and the fit from: 
    // ICRU Report N49, 1993. Ziegler's model for protons.
    // Proton kinetic energy for parametrisation (keV/amu)  

    G4double T = kineticEnergy/(keV*protonMassAMU) ; 

    static const G4float a[11][5] = {
   {1.187E+1f, 1.343E+1f, 1.069E+4f, 7.723E+2f, 2.153E-2f},
   {7.802E+0f, 8.814E+0f, 8.303E+3f, 7.446E+2f, 7.966E-3f}, 
   {7.294E+0f, 8.284E+0f, 5.010E+3f, 4.544E+2f, 8.153E-3f}, 
   {8.646E+0f, 9.800E+0f, 7.066E+3f, 4.581E+2f, 9.383E-3f}, 
   {1.286E+1f, 1.462E+1f, 5.625E+3f, 2.621E+3f, 3.512E-2f}, 
   {3.229E+1f, 3.696E+1f, 8.918E+3f, 3.244E+3f, 1.273E-1f}, 
   {1.604E+1f, 1.825E+1f, 6.967E+3f, 2.307E+3f, 3.775E-2f}, 
   {8.049E+0f, 9.099E+0f, 9.257E+3f, 3.846E+2f, 1.007E-2f},
   {4.015E+0f, 4.542E+0f, 3.955E+3f, 4.847E+2f, 7.904E-3f}, 
   {4.571E+0f, 5.173E+0f, 4.346E+3f, 4.779E+2f, 8.572E-3f},
   {2.631E+0f, 2.601E+0f, 1.701E+3f, 1.279E+3f, 1.638E-2f} };

    static const G4float atomicWeight[11] = {
    101.96128f, 44.0098f, 16.0426f, 28.0536f, 42.0804f,
    104.1512f, 44.665f, 60.0843f, 18.0152f, 18.0152f, 12.0f};       

    if ( T < 10.0 ) {
      ionloss = ((G4double)(a[iMolecula][0])) * std::sqrt(T) ;
    
    } else if ( T < 10000.0 ) {
      G4double x1 = (G4double)(a[iMolecula][1]);
      G4double x2 = (G4double)(a[iMolecula][2]);
      G4double x3 = (G4double)(a[iMolecula][3]);
      G4double x4 = (G4double)(a[iMolecula][4]);
      G4double slow  = x1 * G4Exp(G4Log(T)* 0.45);
      G4double shigh = G4Log( 1.0 + x3/T  + x4*T ) * x2/T;
      ionloss = slow*shigh / (slow + shigh) ;     
    } 

    ionloss = std::max(ionloss, 0.0);
    if ( 10 == iMolecula ) { 
      static const G4double invLog10 = 1.0/G4Log(10.);

      if (T < 100.0) {
        ionloss *= (1.0+0.023+0.0066*G4Log(T)*invLog10);  
      }
      else if (T < 700.0) {   
        ionloss *=(1.0+0.089-0.0248*G4Log(T-99.)*invLog10);
      } 
      else if (T < 10000.0) {    
        ionloss *=(1.0+0.089-0.0248*G4Log(700.-99.)*invLog10);
      }
    }
    ionloss /= (G4double)atomicWeight[iMolecula];

  // pure material (normally not the case for this function)
  } else if(1 == (material->GetNumberOfElements())) {
    G4double z = material->GetZ() ;
    ionloss = ElectronicStoppingPower( z, kineticEnergy ) ;  
  }
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::ElectronicStoppingPower(G4double z,
                                               G4double kineticEnergy) const
{
  G4double ionloss ;
  G4int i = std::min(std::max(G4lrint(z)-1,0),91);  // index of atom

  // The data and the fit from: 
  // ICRU Report 49, 1993. Ziegler's type of parametrisations.
  // Proton kinetic energy for parametrisation (keV/amu)  

  G4double T = kineticEnergy/(keV*protonMassAMU) ; 
  
  static const G4float a[92][5] = {
   {1.254E+0f, 1.440E+0f, 2.426E+2f, 1.200E+4f, 1.159E-1f},
   {1.229E+0f, 1.397E+0f, 4.845E+2f, 5.873E+3f, 5.225E-2f},
   {1.411E+0f, 1.600E+0f, 7.256E+2f, 3.013E+3f, 4.578E-2f},
   {2.248E+0f, 2.590E+0f, 9.660E+2f, 1.538E+2f, 3.475E-2f},
   {2.474E+0f, 2.815E+0f, 1.206E+3f, 1.060E+3f, 2.855E-2f},
   {2.631E+0f, 2.601E+0f, 1.701E+3f, 1.279E+3f, 1.638E-2f},
   {2.954E+0f, 3.350E+0f, 1.683E+3f, 1.900E+3f, 2.513E-2f},
   {2.652E+0f, 3.000E+0f, 1.920E+3f, 2.000E+3f, 2.230E-2f},
   {2.085E+0f, 2.352E+0f, 2.157E+3f, 2.634E+3f, 1.816E-2f},
   {1.951E+0f, 2.199E+0f, 2.393E+3f, 2.699E+3f, 1.568E-2f},
       // Z= 11-20
   {2.542E+0f, 2.869E+0f, 2.628E+3f, 1.854E+3f, 1.472E-2f},
   {3.791E+0f, 4.293E+0f, 2.862E+3f, 1.009E+3f, 1.397E-2f},
   {4.154E+0f, 4.739E+0f, 2.766E+3f, 1.645E+2f, 2.023E-2f},
   {4.914E+0f, 5.598E+0f, 3.193E+3f, 2.327E+2f, 1.419E-2f},
   {3.232E+0f, 3.647E+0f, 3.561E+3f, 1.560E+3f, 1.267E-2f},
   {3.447E+0f, 3.891E+0f, 3.792E+3f, 1.219E+3f, 1.211E-2f},
   {5.301E+0f, 6.008E+0f, 3.969E+3f, 6.451E+2f, 1.183E-2f},
   {5.731E+0f, 6.500E+0f, 4.253E+3f, 5.300E+2f, 1.123E-2f},
   {5.152E+0f, 5.833E+0f, 4.482E+3f, 5.457E+2f, 1.129E-2f},
   {5.521E+0f, 6.252E+0f, 4.710E+3f, 5.533E+2f, 1.112E-2f},
       // Z= 21-30
   {5.201E+0f, 5.884E+0f, 4.938E+3f, 5.609E+2f, 9.995E-3f},
   {4.858E+0f, 5.489E+0f, 5.260E+3f, 6.511E+2f, 8.930E-3f},
   {4.479E+0f, 5.055E+0f, 5.391E+3f, 9.523E+2f, 9.117E-3f},
   {3.983E+0f, 4.489E+0f, 5.616E+3f, 1.336E+3f, 8.413E-3f},
   {3.469E+0f, 3.907E+0f, 5.725E+3f, 1.461E+3f, 8.829E-3f},
   {3.519E+0f, 3.963E+0f, 6.065E+3f, 1.243E+3f, 7.782E-3f},
   {3.140E+0f, 3.535E+0f, 6.288E+3f, 1.372E+3f, 7.361E-3f},
   {3.553E+0f, 4.004E+0f, 6.205E+3f, 5.551E+2f, 8.763E-3f},
   {3.696E+0f, 4.194E+0f, 4.649E+3f, 8.113E+1f, 2.242E-2f},
   {4.210E+0f, 4.750E+0f, 6.953E+3f, 2.952E+2f, 6.809E-3f},
       // Z= 31-40
   {5.041E+0f, 5.697E+0f, 7.173E+3f, 2.026E+2f, 6.725E-3f},
   {5.554E+0f, 6.300E+0f, 6.496E+3f, 1.100E+2f, 9.689E-3f},
   {5.323E+0f, 6.012E+0f, 7.611E+3f, 2.925E+2f, 6.447E-3f},
   {5.874E+0f, 6.656E+0f, 7.395E+3f, 1.175E+2f, 7.684E-3f},
   {6.658E+0f, 7.536E+0f, 7.694E+3f, 2.223E+2f, 6.509E-3f},
   {6.413E+0f, 7.240E+0f, 1.185E+4f, 1.537E+2f, 2.880E-3f},
   {5.694E+0f, 6.429E+0f, 8.478E+3f, 2.929E+2f, 6.087E-3f},
   {6.339E+0f, 7.159E+0f, 8.693E+3f, 3.303E+2f, 6.003E-3f},
   {6.407E+0f, 7.234E+0f, 8.907E+3f, 3.678E+2f, 5.889E-3f},
   {6.734E+0f, 7.603E+0f, 9.120E+3f, 4.052E+2f, 5.765E-3f},
       // Z= 41-50
   {6.901E+0f, 7.791E+0f, 9.333E+3f, 4.427E+2f, 5.587E-3f},
   {6.424E+0f, 7.248E+0f, 9.545E+3f, 4.802E+2f, 5.376E-3f},
   {6.799E+0f, 7.671E+0f, 9.756E+3f, 5.176E+2f, 5.315E-3f},
   {6.109E+0f, 6.887E+0f, 9.966E+3f, 5.551E+2f, 5.151E-3f},
   {5.924E+0f, 6.677E+0f, 1.018E+4f, 5.925E+2f, 4.919E-3f},
   {5.238E+0f, 5.900E+0f, 1.038E+4f, 6.300E+2f, 4.758E-3f},
   // {5.623f,    6.354f,    7160.0f,   337.6f,    0.013940f}, // Ag Ziegler77
   {5.345E+0f, 6.038E+0f, 6.790E+3f, 3.978E+2f, 1.676E-2f}, // Ag ICRU49
   {5.814E+0f, 6.554E+0f, 1.080E+4f, 3.555E+2f, 4.626E-3f},
   {6.229E+0f, 7.024E+0f, 1.101E+4f, 3.709E+2f, 4.540E-3f},
   {6.409E+0f, 7.227E+0f, 1.121E+4f, 3.864E+2f, 4.474E-3f},
       // Z= 51-60
   {7.500E+0f, 8.480E+0f, 8.608E+3f, 3.480E+2f, 9.074E-3f},
   {6.979E+0f, 7.871E+0f, 1.162E+4f, 3.924E+2f, 4.402E-3f},
   {7.725E+0f, 8.716E+0f, 1.183E+4f, 3.948E+2f, 4.376E-3f},
   {8.337E+0f, 9.425E+0f, 1.051E+4f, 2.696E+2f, 6.206E-3f},
   {7.287E+0f, 8.218E+0f, 1.223E+4f, 3.997E+2f, 4.447E-3f},
   {7.899E+0f, 8.911E+0f, 1.243E+4f, 4.021E+2f, 4.511E-3f},
   {8.041E+0f, 9.071E+0f, 1.263E+4f, 4.045E+2f, 4.540E-3f},
   {7.488E+0f, 8.444E+0f, 1.283E+4f, 4.069E+2f, 4.420E-3f},
   {7.291E+0f, 8.219E+0f, 1.303E+4f, 4.093E+2f, 4.298E-3f},
   {7.098E+0f, 8.000E+0f, 1.323E+4f, 4.118E+2f, 4.182E-3f},
       // Z= 61-70
   {6.909E+0f, 7.786E+0f, 1.343E+4f, 4.142E+2f, 4.058E-3f},
   {6.728E+0f, 7.580E+0f, 1.362E+4f, 4.166E+2f, 3.976E-3f},
   {6.551E+0f, 7.380E+0f, 1.382E+4f, 4.190E+2f, 3.877E-3f},
   {6.739E+0f, 7.592E+0f, 1.402E+4f, 4.214E+2f, 3.863E-3f},
   {6.212E+0f, 6.996E+0f, 1.421E+4f, 4.239E+2f, 3.725E-3f},
   {5.517E+0f, 6.210E+0f, 1.440E+4f, 4.263E+2f, 3.632E-3f},
   {5.220E+0f, 5.874E+0f, 1.460E+4f, 4.287E+2f, 3.498E-3f},
   {5.071E+0f, 5.706E+0f, 1.479E+4f, 4.330E+2f, 3.405E-3f},
   {4.926E+0f, 5.542E+0f, 1.498E+4f, 4.335E+2f, 3.342E-3f},
   {4.788E+0f, 5.386E+0f, 1.517E+4f, 4.359E+2f, 3.292E-3f},
       // Z= 71-80
   {4.893E+0f, 5.505E+0f, 1.536E+4f, 4.384E+2f, 3.243E-3f},
   {5.028E+0f, 5.657E+0f, 1.555E+4f, 4.408E+2f, 3.195E-3f},
   {4.738E+0f, 5.329E+0f, 1.574E+4f, 4.432E+2f, 3.186E-3f},
   {4.587E+0f, 5.160E+0f, 1.541E+4f, 4.153E+2f, 3.406E-3f},
   {5.201E+0f, 5.851E+0f, 1.612E+4f, 4.416E+2f, 3.122E-3f},
   {5.071E+0f, 5.704E+0f, 1.630E+4f, 4.409E+2f, 3.082E-3f},
   {4.946E+0f, 5.563E+0f, 1.649E+4f, 4.401E+2f, 2.965E-3f},
   {4.477E+0f, 5.034E+0f, 1.667E+4f, 4.393E+2f, 2.871E-3f},
   //  {4.856f,    5.460f,    18320.0f,  438.5f,    0.002542f}, //Ziegler77
   {4.844E+0f, 5.458E+0f, 7.852E+3f, 9.758E+2f, 2.077E-2f}, //ICRU49
   {4.307E+0f, 4.843E+0f, 1.704E+4f, 4.878E+2f, 2.882E-3f},
       // Z= 81-90
   {4.723E+0f, 5.311E+0f, 1.722E+4f, 5.370E+2f, 2.913E-3f},
   {5.319E+0f, 5.982E+0f, 1.740E+4f, 5.863E+2f, 2.871E-3f},
   {5.956E+0f, 6.700E+0f, 1.780E+4f, 6.770E+2f, 2.660E-3f},
   {6.158E+0f, 6.928E+0f, 1.777E+4f, 5.863E+2f, 2.812E-3f},
   {6.203E+0f, 6.979E+0f, 1.795E+4f, 5.863E+2f, 2.776E-3f},
   {6.181E+0f, 6.954E+0f, 1.812E+4f, 5.863E+2f, 2.748E-3f},
   {6.949E+0f, 7.820E+0f, 1.830E+4f, 5.863E+2f, 2.737E-3f},
   {7.506E+0f, 8.448E+0f, 1.848E+4f, 5.863E+2f, 2.727E-3f},
   {7.648E+0f, 8.609E+0f, 1.866E+4f, 5.863E+2f, 2.697E-3f},
   {7.711E+0f, 8.679E+0f, 1.883E+4f, 5.863E+2f, 2.641E-3f},
       // Z= 91-92
   {7.407E+0f, 8.336E+0f, 1.901E+4f, 5.863E+2f, 2.603E-3f},
   {7.290E+0f, 8.204E+0f, 1.918E+4f, 5.863E+2f, 2.673E-3f}
  };

  G4double fac = 1.0 ;

    // Carbon specific case for E < 40 keV
  if ( T < 40.0 && 5 == i) {
    fac = std::sqrt(T*0.025);
    T = 40.0;  

    // Free electron gas model
  } else if ( T < 10.0 ) { 
    fac = std::sqrt(T*0.1) ;
    T = 10.0;
  }

  // Main parametrisation
  G4double x1 = (G4double)(a[i][1]);
  G4double x2 = (G4double)(a[i][2]);
  G4double x3 = (G4double)(a[i][3]);
  G4double x4 = (G4double)(a[i][4]);
  G4double slow  = x1 * G4Exp(G4Log(T) * 0.45);
  G4double shigh = G4Log( 1.0 + x3/T + x4*T ) * x2/T;
  ionloss = slow*shigh*fac / (slow + shigh);
  
  ionloss = std::max(ionloss, 0.0);
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::DEDX(const G4Material* material, G4double kineticEnergy) 
{
  G4double eloss = 0.0;

  // check DB
  if(material != currentMaterial) {
    currentMaterial = material;
    baseMaterial = material->GetBaseMaterial() 
      ? material->GetBaseMaterial() : material;
    iPSTAR    = -1;
    iMolecula = -1;
    iICRU90 = fICRU90 ? fICRU90->GetIndex(baseMaterial) : -1;
    
    if(iICRU90 < 0) { 
      iPSTAR = fPSTAR->GetIndex(baseMaterial); 
      if(iPSTAR < 0) { HasMaterial(baseMaterial); }
    }
    //G4cout << "%%% " <<material->GetName() << "  iMolecula= " 
    //       << iMolecula << "  iPSTAR= " << iPSTAR 
    //       << "  iICRU90= " << iICRU90<< G4endl; 
  }

  // ICRU90 parameterisation
  if(iICRU90 >= 0) {
    return fICRU90->GetElectronicDEDXforProton(iICRU90, kineticEnergy)
      *material->GetDensity();
  }
  // PSTAR parameterisation
  if( iPSTAR >= 0 ) {
    return fPSTAR->GetElectronicDEDX(iPSTAR, kineticEnergy)
      *material->GetDensity();

  } 
  const std::size_t numberOfElements = material->GetNumberOfElements();
  const G4double* theAtomicNumDensityVector =
                                 material->GetAtomicNumDensityVector();
  

  if(iMolecula >= 0) {
    eloss = StoppingPower(baseMaterial, kineticEnergy)*
                          material->GetDensity()/amu;

    // Pure material ICRU49 paralmeterisation
  } else if(1 == numberOfElements) {

    G4double z = material->GetZ();
    eloss = ElectronicStoppingPower(z, kineticEnergy)
                               * (material->GetTotNbOfAtomsPerVolume());


  // Experimental data exist only for kinetic energy 125 keV
  } else if( MolecIsInZiegler1988(material) ) { 

    // Loop over elements - calculation based on Bragg's rule 
    G4double eloss125 = 0.0 ;
    const G4ElementVector* theElementVector =
                           material->GetElementVector();
  
    //  Loop for the elements in the material
    for (std::size_t i=0; i<numberOfElements; ++i) {
      const G4Element* element = (*theElementVector)[i] ;
      G4double z = element->GetZ() ;
      eloss    += ElectronicStoppingPower(z,kineticEnergy)
                                    * theAtomicNumDensityVector[i] ;
      eloss125 += ElectronicStoppingPower(z,125.0*keV)
                                    * theAtomicNumDensityVector[i] ;
    }      

    // Chemical factor is taken into account
    eloss *= ChemicalFactor(kineticEnergy, eloss125) ;
 
  // Brugg's rule calculation
  } else {
    const G4ElementVector* theElementVector =
                           material->GetElementVector() ;
  
    //  loop for the elements in the material
    for (std::size_t i=0; i<numberOfElements; ++i)
    {
      const G4Element* element = (*theElementVector)[i] ;
      eloss   += ElectronicStoppingPower(element->GetZ(), kineticEnergy)
                                   * theAtomicNumDensityVector[i];
    }      
  }
  return eloss*theZieglerFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4BraggModel::MolecIsInZiegler1988(const G4Material* material) 
{
  // The list of molecules from
  // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
  // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.
  
  G4String myFormula = G4String(" ") ;
  const G4String chFormula = material->GetChemicalFormula() ;
  if (myFormula == chFormula ) { return false; }
  
  //  There are no evidence for difference of stopping power depended on
  //  phase of the compound except for water. The stopping power of the 
  //  water in gas phase can be predicted using Bragg's rule.
  //  
  //  No chemical factor for water-gas 
   
  myFormula = G4String("H_2O") ;
  const G4State theState = material->GetState() ;
  if( theState == kStateGas && myFormula == chFormula) return false ;
    

  // The coffecient from Table.4 of Ziegler & Manoyan
  static const G4float HeEff = 2.8735f;
  
  static const std::size_t numberOfMolecula = 53;
  static const G4String nameOfMol[numberOfMolecula] = {
    "H_2O",      "C_2H_4O",    "C_3H_6O",  "C_2H_2",             "C_H_3OH",
    "C_2H_5OH",  "C_3H_7OH",   "C_3H_4",   "NH_3",               "C_14H_10",
    "C_6H_6",    "C_4H_10",    "C_4H_6",   "C_4H_8O",            "CCl_4",
    "CF_4",      "C_6H_8",     "C_6H_12",  "C_6H_10O",           "C_6H_10",
    "C_8H_16",   "C_5H_10",    "C_5H_8",   "C_3H_6-Cyclopropane","C_2H_4F_2",
    "C_2H_2F_2", "C_4H_8O_2",  "C_2H_6",   "C_2F_6",             "C_2H_6O",
    "C_3H_6O",   "C_4H_10O",   "C_2H_4",   "C_2H_4O",            "C_2H_4S",
    "SH_2",      "CH_4",       "CCLF_3",   "CCl_2F_2",           "CHCl_2F",
    "(CH_3)_2S", "N_2O",       "C_5H_10O", "C_8H_6",             "(CH_2)_N",
    "(C_3H_6)_N","(C_8H_8)_N", "C_3H_8",   "C_3H_6-Propylene",   "C_3H_6O",
    "C_3H_6S",   "C_4H_4S",    "C_7H_8"
  };

  static const G4float expStopping[numberOfMolecula] = {
     66.1f,  190.4f, 258.7f,  42.2f, 141.5f,
    210.9f,  279.6f, 198.8f,  31.0f, 267.5f,
    122.8f,  311.4f, 260.3f, 328.9f, 391.3f,
    206.6f,  374.0f, 422.0f, 432.0f, 398.0f,
    554.0f,  353.0f, 326.0f,  74.6f, 220.5f,
    197.4f,  362.0f, 170.0f, 330.5f, 211.3f,
    262.3f,  349.6f,  51.3f, 187.0f, 236.9f,
    121.9f,   35.8f, 247.0f, 292.6f, 268.0f,
    262.3f,   49.0f, 398.9f, 444.0f,  22.91f,
     68.0f,  155.0f,  84.0f,  74.2f, 254.7f,
    306.8f,  324.4f, 420.0f
  } ;

  static const G4float expCharge[numberOfMolecula] = {
    HeEff, HeEff, HeEff,  1.0f, HeEff,
    HeEff, HeEff, HeEff,  1.0f,  1.0f,
     1.0f, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff, HeEff,  1.0f, HeEff,
    HeEff, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff,  1.0f, HeEff, HeEff,
    HeEff,  1.0f, HeEff, HeEff, HeEff,
    HeEff,  1.0f, HeEff, HeEff,  1.0f,
     1.0f,  1.0f,  1.0f,  1.0f, HeEff,
    HeEff, HeEff, HeEff
  } ;

  static const G4int numberOfAtomsPerMolecula[numberOfMolecula] = {
    3,  7, 10,  4,  6,
    9, 12,  7,  4, 24,
    12,14, 10, 13,  5,
    5, 14, 18, 17, 17,
    24,15, 13,  9,  8,
    6, 14,  8,  8,  9,
    10,15,  6,  7,  7,
    3,  5,  5,  5,  5,
    9,  3, 16, 14,  3,
    9, 16, 11,  9, 10,
    10, 9,  15};

  // Search for the compaund in the table
  for (std::size_t i=0; i<numberOfMolecula; ++i) {
    if(chFormula == nameOfMol[i]) {
      expStopPower125 = ((G4double)expStopping[i])
        * (material->GetTotNbOfAtomsPerVolume()) /
        ((G4double)(expCharge[i] * numberOfAtomsPerMolecula[i]));
      return true;
    }
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggModel::ChemicalFactor(G4double kineticEnergy, 
                                      G4double eloss125) const
{
  // Approximation of Chemical Factor according to
  // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
  // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.

  static const G4double gamma25  = 1.0 + 25.0*keV /proton_mass_c2;
  static const G4double gamma125 = 1.0 + 125.0*keV/proton_mass_c2;
  static const G4double beta25   = std::sqrt(1.0 - 1.0/(gamma25*gamma25));
  static const G4double beta125  = std::sqrt(1.0 - 1.0/(gamma125*gamma125));
  static const G4double f12525   = 1.0 + G4Exp( 1.48*(beta125/beta25 - 7.0) );
  
  G4double gamma = 1.0 + kineticEnergy/proton_mass_c2;
  G4double beta  = std::sqrt(1.0 - 1.0/(gamma*gamma));
  
  G4double factor = 1.0 + (expStopPower125/eloss125 - 1.0) * f12525/
    (1.0 + G4Exp( 1.48 * ( beta/beta25    - 7.0 ) ) );

  return factor ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4BraggIonModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.10.2004
//
// Modifications:
// 11-05-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 29-11-05 Do not use G4Alpha class (V.Ivantchenko)
// 15-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 25-04-06 Add stopping data from ASTAR (V.Ivanchenko)
// 23-10-06 Reduce lowestKinEnergy to 0.25 keV (V.Ivanchenko)
// 12-08-08 Added methods GetParticleCharge, GetChargeSquareRatio, 
//          CorrectionsAlongStep needed for ions(V.Ivanchenko)
//

// Class Description:
//
// Implementation of energy loss and delta-electron production by
// slow charged heavy particles

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4BraggIonModel.hh"
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

using namespace std;

G4ASTARStopping* G4BraggIonModel::fASTAR = nullptr;

G4BraggIonModel::G4BraggIonModel(const G4ParticleDefinition* p,
                                 const G4String& nam)
  : G4VEmModel(nam),
    theElectron(G4Electron::Electron()),
    HeMass(3.727417*CLHEP::GeV),
    theZieglerFactor(CLHEP::eV*CLHEP::cm2*1.0e-15),
    lowestKinEnergy(0.25*CLHEP::keV)
{
  SetHighEnergyLimit(2.0*CLHEP::MeV);

  rateMassHe2p = HeMass/CLHEP::proton_mass_c2;
  massFactor = 1000.*CLHEP::amu_c2/HeMass;

  if(nullptr != p) { SetParticle(p); }
  else  { SetParticle(theElectron); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BraggIonModel::~G4BraggIonModel()
{
  if(IsMaster()) { delete fASTAR; fASTAR = nullptr; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggIonModel::Initialise(const G4ParticleDefinition* p,
                                 const G4DataVector&)
{
  if(p != particle) { SetParticle(p); }

  // always false before the run
  SetDeexcitationFlag(false);

  // initialise once
  if(nullptr == fParticleChange) {
    const G4String& pname = particle->GetParticleName();
    if(IsMaster()) {
      if(pname == "proton" || pname == "GenericIon" || pname == "alpha") {
	if(nullptr == fASTAR)  { fASTAR = new G4ASTARStopping(); }
	fASTAR->Initialise(); 

	if(G4EmParameters::Instance()->UseICRU90Data()) {
	  fICRU90 = G4NistManager::Instance()->GetICRU90StoppingData();
	  fICRU90->Initialise();
	}
      }
    }
    if(pname == "alpha") { isAlpha = true; }

    if(UseAngularGeneratorFlag() && nullptr == GetAngularDistribution()) {
      SetAngularDistribution(new G4DeltaAngle());
    }
    corr = G4LossTableManager::Instance()->EmCorrections();

    fParticleChange = GetParticleChangeForLoss();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::MinEnergyCut(const G4ParticleDefinition*,
                                       const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::GetChargeSquareRatio(const G4ParticleDefinition* p,
                                               const G4Material* mat,
                                               G4double kineticEnergy)
{
  return corr->EffectiveChargeSquareRatio(p,mat,kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::GetParticleCharge(const G4ParticleDefinition* p,
                                            const G4Material* mat,
                                            G4double kineticEnergy)
{
  return corr->GetParticleCharge(p,mat,kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::ComputeCrossSectionPerElectron(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double minKinEnergy,
                                                 G4double maxKinEnergy)
{
  G4double cross = 0.0;
  const G4double tmax      = MaxSecondaryEnergy(p, kineticEnergy);
  const G4double maxEnergy = std::min(tmax, maxKinEnergy);
  const G4double cutEnergy = std::max(lowestKinEnergy*massRate, minKinEnergy);
  
  if(cutEnergy < tmax) {

    const G4double energy  = kineticEnergy + mass;
    const G4double energy2 = energy*energy;
    const G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;

    cross = (maxEnergy - cutEnergy)/(cutEnergy*maxEnergy) 
      - beta2*G4Log(maxEnergy/cutEnergy)/tmax;
    if( 0.0 < spin ) { cross += 0.5*(maxEnergy - cutEnergy)/energy2; }

    cross *= CLHEP::twopi_mc2_rcl2*chargeSquare/beta2;
    cross = std::max(cross, 0.0);
  }
  //   G4cout << "BR: e= " << kineticEnergy << " tmin= " << cutEnergy 
  //          << " tmax= " << tmax << " cross= " << cross << G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition* p,
                                                 G4double kinEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double sigma = 
    Z*ComputeCrossSectionPerElectron(p,kinEnergy,cutEnergy,maxEnergy);
  if(isAlpha) {
    sigma *= (HeEffChargeSquare(Z, kinEnergy/CLHEP::MeV)/chargeSquare);
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::CrossSectionPerVolume(
                                           const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kinEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double sigma = material->GetElectronDensity()* 
    ComputeCrossSectionPerElectron(p,kinEnergy,cutEnergy,maxEnergy);
  if(isAlpha) {
    const G4double zeff = material->GetTotNbOfElectPerVolume()/
      material->GetTotNbOfAtomsPerVolume();
    sigma *= (HeEffChargeSquare(zeff, kinEnergy/CLHEP::MeV)/chargeSquare);
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::ComputeDEDXPerVolume(const G4Material* material,
                                               const G4ParticleDefinition* p,
                                               G4double kineticEnergy,
                                               G4double minKinEnergy)
{
  const G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  const G4double tmin = std::max(lowestKinEnergy*massRate, minKinEnergy);
  G4double dedx = 0.0;

  // T is alpha energy
  G4double T = kineticEnergy;
  const G4double zeff = material->GetTotNbOfElectPerVolume()/
    material->GetTotNbOfAtomsPerVolume();
  heChargeSquare = HeEffChargeSquare(zeff, T/CLHEP::MeV);
  if(!isAlpha) { T *= rateMassHe2p; }

  if(T < lowestKinEnergy) {
    dedx = DEDX(material, lowestKinEnergy)*std::sqrt(T/lowestKinEnergy);
  } else {
    dedx = DEDX(material, T);
  }
  if(!isAlpha) { dedx /= heChargeSquare; }
  if (tmin < tmax) {
    const G4double tau = kineticEnergy/mass;
    const G4double x   = tmin/tmax;

    G4double del = 
      (G4Log(x)*(tau + 1.)*(tau + 1.)/(tau * (tau + 2.0)) + 1.0 - x) * 
      CLHEP::twopi_mc2_rcl2*material->GetElectronDensity();
    if(isAlpha) { del *= heChargeSquare; }
    dedx += del;
  }
  dedx = std::max(dedx, 0.0);
  /*
  G4cout << "BraggIon: tkin(MeV) = " << tkin/MeV << " dedx(MeV*cm^2/g) = " 
         << dedx*gram/(MeV*cm2*material->GetDensity()) 
         << " q2 = " << chargeSquare <<  G4endl;
  */
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggIonModel::CorrectionsAlongStep(const G4MaterialCutsCouple* couple,
                                           const G4DynamicParticle* dp,
                                           const G4double&,
                                           G4double& eloss)
{
  // no correction for alpha
  if(isAlpha) { return; }

  // no correction at a small step at the last step
  const G4double preKinEnergy = dp->GetKineticEnergy();
  if(eloss >= preKinEnergy || eloss < preKinEnergy*0.05) { return; }

  // corrections only for ions
  const G4ParticleDefinition* p = dp->GetDefinition();
  if(p != particle) { SetParticle(p); }

  // effective energy and charge at a step
  const G4Material* mat = couple->GetMaterial();
  const G4double e = std::max(preKinEnergy - eloss*0.5, preKinEnergy*0.5);
  const G4double q20 = corr->EffectiveChargeSquareRatio(p, mat, preKinEnergy);
  const G4double q2 = corr->EffectiveChargeSquareRatio(p, mat, e);
  const G4double qfactor = q2/q20;
  /*    
    G4cout << "G4BraggIonModel::CorrectionsAlongStep: Epre(MeV)="
    << preKinEnergy << " Eeff(MeV)=" << e
    << " eloss=" << eloss << " elossnew=" << eloss*qfactor 
    << " qfactor=" << qfactor << " Qpre=" << q20 
    << p->GetParticleName() <<G4endl;
  */  
  eloss *= qfactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggIonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
                                        const G4MaterialCutsCouple* couple,
                                        const G4DynamicParticle* dp,
                                        G4double minEnergy,
                                        G4double maxEnergy)
{
  const G4double tmax = MaxSecondaryKinEnergy(dp);
  const G4double xmax = std::min(tmax, maxEnergy);
  const G4double xmin = std::max(lowestKinEnergy*massRate, minEnergy);
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
        G4cout << "G4BraggIonModel::SampleSecondary Warning! "
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
      sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
    G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
      (deltaMomentum * dp->GetTotalMomentum());
    if(cost > 1.0) { cost = 1.0; }
    G4double sint = sqrt((1.0 - cost)*(1.0 + cost));

    G4double phi = twopi*rndmEngineMod->flat();

    deltaDirection.set(sint*cos(phi),sint*sin(phi), cost) ;
    deltaDirection.rotateUz(dp->GetMomentumDirection());
  }  

  // create G4DynamicParticle object for delta ray
  auto delta = new G4DynamicParticle(theElectron,deltaDirection,deltaKinEnergy);

  vdp->push_back(delta);

  // Change kinematics of primary particle
  kineticEnergy -= deltaKinEnergy;
  G4ThreeVector finalP = dp->GetMomentum() - delta->GetMomentum();
  finalP               = finalP.unit();

  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(finalP);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::MaxSecondaryEnergy(const G4ParticleDefinition* pd,
                                             G4double kinEnergy)
{
  if(pd != particle) { SetParticle(pd); }
  G4double tau  = kinEnergy/mass;
  G4double tmax = 2.0*electron_mass_c2*tau*(tau + 2.) /
                  (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggIonModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  G4double q = particle->GetPDGCharge()/CLHEP::eplus;
  chargeSquare = q*q;
  massRate = mass/CLHEP::proton_mass_c2;
  ratio = CLHEP::electron_mass_c2/mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4BraggIonModel::HasMaterial(const G4Material* mat) const
{
  const G4String& chFormula = mat->GetChemicalFormula();
  if(chFormula.empty()) { return -1; }

  // ICRU Report N49, 1993. Ziegler model for He.
  
  static const G4int numberOfMolecula = 11;
  static const G4String molName[numberOfMolecula] = {
    "CaF_2",  "Cellulose_Nitrate",  "LiF", "Policarbonate",  
    "(C_2H_4)_N-Polyethylene",  "(C_2H_4)_N-Polymethly_Methacralate",
    "Polysterene", "SiO_2", "NaI", "H_2O",
    "Graphite" };

  // Search for the material in the table
  for (G4int i=0; i<numberOfMolecula; ++i) {
    if (chFormula == molName[i]) {  
      return i;
    }
  }
  return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::StoppingPower(const G4Material* material,
                                        const G4double kineticEnergy) const
{
  G4double ionloss = 0.0 ;

  if (iMolecula >= 0) {
  
    // The data and the fit from: 
    // ICRU Report N49, 1993. Ziegler's model for alpha
    // He energy in internal units of parametrisation formula (MeV)
    // Input scaled energy of a proton or GenericIon

    //    G4double T = kineticEnergy*rateMassHe2p/CLHEP::MeV;
    G4double T = kineticEnergy/CLHEP::MeV;

    static const G4float a[11][5] = {
       {9.43672f, 0.54398f, 84.341f,  1.3705f, 57.422f},
       {67.1503f, 0.41409f, 404.512f, 148.97f, 20.99f},
       {5.11203f, 0.453f,   36.718f,  50.6f,   28.058f}, 
       {61.793f,  0.48445f, 361.537f, 57.889f, 50.674f},
       {7.83464f, 0.49804f, 160.452f, 3.192f,  0.71922f},
       {19.729f,  0.52153f, 162.341f, 58.35f,  25.668f}, 
       {26.4648f, 0.50112f, 188.913f, 30.079f, 16.509f},
       {7.8655f,  0.5205f,  63.96f,   51.32f,  67.775f},
       {8.8965f,  0.5148f,  339.36f,  1.7205f, 0.70423f},
       {2.959f,   0.53255f, 34.247f,  60.655f, 15.153f}, 
       {3.80133f, 0.41590f, 12.9966f, 117.83f, 242.28f} };   

    static const G4double atomicWeight[11] = {
       101.96128f, 44.0098f, 16.0426f, 28.0536f, 42.0804f,
       104.1512f,  44.665f,  60.0843f, 18.0152f, 18.0152f, 12.0f};       

    G4int i = iMolecula;

    G4double slow = (G4double)(a[i][0]);

    G4double x1 = (G4double)(a[i][1]);
    G4double x2 = (G4double)(a[i][2]);
    G4double x3 = (G4double)(a[i][3]);
    G4double x4 = (G4double)(a[i][4]);

    // Free electron gas model
    if ( T < 0.001 ) {
      G4double shigh = G4Log( 1.0 + x3*1000.0 + x4*0.001 ) *x2*1000.0;
      ionloss  = slow*shigh / (slow + shigh) ;
      ionloss *= sqrt(T*1000.0) ;

      // Main parametrisation
    } else {
      slow  *= G4Exp(G4Log(T*1000.0)*x1) ;
      G4double shigh = G4Log( 1.0 + x3/T + x4*T ) * x2/T ;
      ionloss = slow*shigh / (slow + shigh) ;
       /*
         G4cout << "## " << i << ". T= " << T << " slow= " << slow
         << " a0= " << a[i][0] << " a1= " << a[i][1] 
         << " shigh= " << shigh 
         << " dedx= " << ionloss << " q^2= " <<  HeEffChargeSquare(z, T*MeV)
         << G4endl;
       */
    }
    ionloss = std::max(ionloss, 0.0);

    // He effective charge
    ionloss /= (heChargeSquare*atomicWeight[iMolecula]);

  // pure material (normally not the case for this function)
  } else if(1 == (material->GetNumberOfElements())) {
    const G4double z = material->GetZ() ;
    ionloss = ElectronicStoppingPower( z, kineticEnergy ) ;  
  }
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4BraggIonModel::ElectronicStoppingPower(const G4double z,
                                         const G4double kineticEnergy) const
{
  G4double ionloss ;
  G4int i = std::min(std::max(G4lrint(z)-1,0),91);  // index of atom
  //G4cout << "ElectronicStoppingPower z=" << z << " i=" << i 
  // << " E=" << kineticEnergy << G4endl;
  // The data and the fit from:
  // ICRU Report 49, 1993. Ziegler's type of parametrisations.
  // Proton kinetic energy for parametrisation (keV/amu)
  // He energy in internal units of parametrisation formula (MeV)
  //G4double T = kineticEnergy*rateMassHe2p/CLHEP::MeV;
  G4double T = kineticEnergy/CLHEP::MeV;

  static const G4float a[92][5] = {
    {  0.35485f, 0.6456f, 6.01525f,  20.8933f, 4.3515f
   },{ 0.58f,    0.59f,   6.3f,      130.0f,   44.07f
   },{ 1.42f,    0.49f,   12.25f,    32.0f,    9.161f
   },{ 2.206f,   0.51f,   15.32f,    0.25f,    8.995f //Be Ziegler77
       // },{ 2.1895f,  0.47183,7.2362f,   134.30f,  197.96f //Be from ICRU
   },{ 3.691f,   0.4128f, 18.48f,    50.72f,   9.0f
   },{ 3.83523f, 0.42993f,12.6125f,  227.41f,  188.97f
       // },{ 1.9259f,  0.5550f, 27.15125f, 26.0665f, 6.2768f //too many digits
   },{ 1.9259f,  0.5550f, 27.1513f,  26.0665f, 6.2768f
   },{ 2.81015f, 0.4759f, 50.0253f,  10.556f,  1.0382f
   },{ 1.533f,   0.531f,  40.44f,    18.41f,   2.718f
   },{ 2.303f,   0.4861f, 37.01f,    37.96f,   5.092f
       // Z= 11-20
   },{ 9.894f,   0.3081f, 23.65f,    0.384f,   92.93f
   },{ 4.3f,     0.47f,   34.3f,     3.3f,     12.74f
   },{ 2.5f,     0.625f,  45.7f,     0.1f,     4.359f
   },{ 2.1f,     0.65f,   49.34f,    1.788f,   4.133f
   },{ 1.729f,   0.6562f, 53.41f,    2.405f,   3.845f
   },{ 1.402f,   0.6791f, 58.98f,    3.528f,   3.211f
   },{ 1.117f,   0.7044f, 69.69f,    3.705f,   2.156f
   },{ 2.291f,   0.6284f, 73.88f,    4.478f,   2.066f
   },{ 8.554f,   0.3817f, 83.61f,    11.84f,   1.875f
   },{ 6.297f,   0.4622f, 65.39f,    10.14f,   5.036f
       // Z= 21-30     
   },{ 5.307f,   0.4918f, 61.74f,    12.4f,    6.665f
   },{ 4.71f,    0.5087f, 65.28f,    8.806f,   5.948f
   },{ 6.151f,   0.4524f, 83.0f,     18.31f,   2.71f
   },{ 6.57f,    0.4322f, 84.76f,    15.53f,   2.779f
   },{ 5.738f,   0.4492f, 84.6f,     14.18f,   3.101f
   },{ 5.013f,   0.4707f, 85.8f,     16.55f,   3.211f
   },{ 4.32f,    0.4947f, 76.14f,    10.85f,   5.441f
   },{ 4.652f,   0.4571f, 80.73f,    22.0f,    4.952f
   },{ 3.114f,   0.5236f, 76.67f,    7.62f,    6.385f
   },{ 3.114f,   0.5236f, 76.67f,    7.62f,    7.502f
       // Z= 31-40
   },{ 3.114f,   0.5236f, 76.67f,    7.62f,    8.514f
   },{ 5.746f,   0.4662f, 79.24f,    1.185f,   7.993f
   },{ 2.792f,   0.6346f, 106.1f,    0.2986f,  2.331f
   },{ 4.667f,   0.5095f, 124.3f,    2.102f,   1.667f
   },{ 2.44f,    0.6346f, 105.0f,    0.83f,    2.851f
   },{ 1.413f,   0.7377f, 147.9f,    1.466f,   1.016f
   },{ 11.72f,   0.3826f, 102.8f,    9.231f,   4.371f
   },{ 7.126f,   0.4804f, 119.3f,    5.784f,   2.454f
   },{ 11.61f,   0.3955f, 146.7f,    7.031f,   1.423f
   },{ 10.99f,   0.41f,   163.9f,    7.1f,     1.052f
       // Z= 41-50
   },{ 9.241f,   0.4275f, 163.1f,    7.954f,   1.102f
   },{ 9.276f,   0.418f,  157.1f,    8.038f,   1.29f
   },{ 3.999f,   0.6152f, 97.6f,     1.297f,   5.792f
   },{ 4.306f,   0.5658f, 97.99f,    5.514f,   5.754f
   },{ 3.615f,   0.6197f, 86.26f,    0.333f,   8.689f
   },{ 5.8f,     0.49f,   147.2f,    6.903f,   1.289f
   },{ 5.6f,     0.49f,   130.0f,    10.0f,    2.844f
   },{ 3.55f,    0.6068f, 124.7f,    1.112f,   3.119f
   },{ 3.6f,     0.62f,   105.8f,    0.1692f,  6.026f
   },{ 5.4f,     0.53f,   103.1f,    3.931f,   7.767f
       // Z= 51-60
   },{ 3.97f,    0.6459f, 131.8f,    0.2233f,  2.723f
   },{ 3.65f,    0.64f,   126.8f,    0.6834f,  3.411f
   },{ 3.118f,   0.6519f, 164.9f,    1.208f,   1.51f
   },{ 3.949f,   0.6209f, 200.5f,    1.878f,   0.9126f
   },{ 14.4f,    0.3923f, 152.5f,    8.354f,   2.597f
   },{ 10.99f,   0.4599f, 138.4f,    4.811f,   3.726f
   },{ 16.6f,    0.3773f, 224.1f,    6.28f,    0.9121f
   },{ 10.54f,   0.4533f, 159.3f,    4.832f,   2.529f
   },{ 10.33f,   0.4502f, 162.0f,    5.132f,   2.444f
   },{ 10.15f,   0.4471f, 165.6f,    5.378f,   2.328f
       // Z= 61-70
   },{ 9.976f,   0.4439f, 168.0f,    5.721f,   2.258f
   },{ 9.804f,   0.4408f, 176.2f,    5.675f,   1.997f
   },{ 14.22f,   0.363f,  228.4f,    7.024f,   1.016f
   },{ 9.952f,   0.4318f, 233.5f,    5.065f,   0.9244f
   },{ 9.272f,   0.4345f, 210.0f,    4.911f,   1.258f
   },{ 10.13f,   0.4146f, 225.7f,    5.525f,   1.055f
   },{ 8.949f,   0.4304f, 213.3f,    5.071f,   1.221f
   },{ 11.94f,   0.3783f, 247.2f,    6.655f,   0.849f
   },{ 8.472f,   0.4405f, 195.5f,    4.051f,   1.604f
   },{ 8.301f,   0.4399f, 203.7f,    3.667f,   1.459f
       // Z= 71-80
   },{ 6.567f,   0.4858f, 193.0f,    2.65f,    1.66f
   },{ 5.951f,   0.5016f, 196.1f,    2.662f,   1.589f
   },{ 7.495f,   0.4523f, 251.4f,    3.433f,   0.8619f
   },{ 6.335f,   0.4825f, 255.1f,    2.834f,   0.8228f
   },{ 4.314f,   0.5558f, 214.8f,    2.354f,   1.263f
   },{ 4.02f,    0.5681f, 219.9f,    2.402f,   1.191f
   },{ 3.836f,   0.5765f, 210.2f,    2.742f,   1.305f
   },{ 4.68f,    0.5247f, 244.7f,    2.749f,   0.8962f
   },{ 2.892f,   0.6204f, 208.6f,    2.415f,   1.416f //Au Z77
       // },{ 3.223f,   0.5883f, 232.7f,   2.954f,    1.05  //Au ICRU
   },{ 2.892f,   0.6204f, 208.6f,    2.415f,   1.416f
       // Z= 81-90
   },{ 4.728f,   0.5522f, 217.0f,    3.091f,   1.386f
   },{ 6.18f,    0.52f,   170.0f,    4.0f,     3.224f
   },{ 9.0f,     0.47f,   198.0f,    3.8f,     2.032f
   },{ 2.324f,   0.6997f, 216.0f,    1.599f,   1.399f
   },{ 1.961f,   0.7286f, 223.0f,    1.621f,   1.296f
   },{ 1.75f,    0.7427f, 350.1f,    0.9789f,  0.5507f
   },{ 10.31f,   0.4613f, 261.2f,    4.738f,   0.9899f
   },{ 7.962f,   0.519f,  235.7f,    4.347f,   1.313f
   },{ 6.227f,   0.5645f, 231.9f,    3.961f,   1.379f
   },{ 5.246f,   0.5947f, 228.6f,    4.027f,   1.432f
       // Z= 91-92
   },{ 5.408f,   0.5811f, 235.7f,    3.961f,   1.358f
   },{ 5.218f,   0.5828f, 245.0f,    3.838f,   1.25f}
  };

  G4double slow = (G4double)(a[i][0]);

  G4double x1 = (G4double)(a[i][1]);
  G4double x2 = (G4double)(a[i][2]);
  G4double x3 = (G4double)(a[i][3]);
  G4double x4 = (G4double)(a[i][4]);

  // Free electron gas model
  if ( T < 0.001 ) {
    G4double shigh = G4Log( 1.0 + x3*1000.0 + x4*0.001 )* x2*1000.0;
    ionloss  = slow*shigh*std::sqrt(T*1000.0)  / (slow + shigh) ;

  // Main parametrisation
  } else {
    slow  *= G4Exp(G4Log(T*1000.0)*x1);
    G4double shigh = G4Log( 1.0 + x3/T + x4*T ) * x2/T;
    ionloss = slow*shigh / (slow + shigh) ;
    /*
    G4cout << "## " << i << ". T= " << T << " slow= " << slow
           << " a0= " << a[i][0] << " a1= " << a[i][1] 
           << " shigh= " << shigh 
           << " dedx= " << ionloss << " q^2= " <<  HeEffChargeSquare(z, T) 
           << G4endl;
    */
  }
  ionloss = std::max(ionloss, 0.0);

  // He effective charge
  // ionloss /= heChargeSquare;
  // G4cout << ionloss << G4endl;
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::DEDX(const G4Material* material,
                               const G4double aEnergy)
{
  // aEnergy is energy of alpha
  G4double eloss = 0.0;
  // check DB
  if(material != currentMaterial) {
    currentMaterial = material;
    baseMaterial = material->GetBaseMaterial() 
      ? material->GetBaseMaterial() : material;
    iASTAR    = -1;
    iMolecula = -1;
    iICRU90 = (nullptr != fICRU90) ? fICRU90->GetIndex(baseMaterial) : -1;
    
    if(iICRU90 < 0) { 
      iASTAR = fASTAR->GetIndex(baseMaterial); 
      if(iASTAR < 0) { iMolecula = HasMaterial(baseMaterial); }
    }
    /*    
    G4cout << "%%% " <<material->GetName() << "  iMolecula= " 
           << iMolecula << "  iASTAR= " << iASTAR 
           << "  iICRU90= " << iICRU90<< G4endl; 
    */
  }
  // ICRU90 
  if(iICRU90 >= 0) {
    eloss = fICRU90->GetElectronicDEDXforAlpha(iICRU90, aEnergy);
    if(eloss > 0.0) { return eloss*material->GetDensity(); }
  }
  // ASTAR
  if( iASTAR >= 0 ) {
    eloss = fASTAR->GetElectronicDEDX(iASTAR, aEnergy);
    /*
    G4cout << "ASTAR:  E=" << aEnergy 
	   << " dedx=" << eloss*material->GetDensity() 
	   << "  " << particle->GetParticleName() << G4endl;
    */
    if(eloss > 0.0) { return eloss*material->GetDensity(); }
  }

  const std::size_t numberOfElements = material->GetNumberOfElements();
  const G4double* theAtomicNumDensityVector =
    material->GetAtomicNumDensityVector();

  if(iMolecula >= 0) {

    eloss = StoppingPower(baseMaterial, aEnergy)*material->GetDensity()/amu;

    // pure material
  } else if(1 == numberOfElements) {

    const G4double z = material->GetZ();
    eloss = ElectronicStoppingPower(z, aEnergy)
                               * (material->GetTotNbOfAtomsPerVolume());

  // Brugg's rule calculation
  } else {
    const G4ElementVector* theElmVector = material->GetElementVector();

    //  loop for the elements in the material
    for (std::size_t i=0; i<numberOfElements; ++i) {
      const G4Element* element = (*theElmVector)[i];
      eloss += ElectronicStoppingPower(element->GetZ(), aEnergy)
	* theAtomicNumDensityVector[i];
    }
  }
  return eloss*theZieglerFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4BraggIonModel::HeEffChargeSquare(const G4double z,
                                   const G4double kinEnergyHeInMeV) const
{
  // The aproximation of He effective charge from:
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985

  static const G4double c[6] = {0.2865,  0.1266, -0.001429,
                                0.02402,-0.01135, 0.001475};

  G4double e = std::max(0.0, G4Log(kinEnergyHeInMeV*massFactor));
  G4double x = c[0] ;
  G4double y = 1.0 ;
  for (G4int i=1; i<6; ++i) {
    y *= e;
    x += y * c[i];
  }

  G4double w = 7.6 -  e ;
  w = 1.0 + (0.007 + 0.00005*z) * G4Exp( -w*w ) ;
  w = 4.0 * (1.0 - G4Exp(-x)) * w * w ;

  return w;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


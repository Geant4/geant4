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
// GEANT4 Class header file
//
//
// File name:     G4BetheBlochModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 04-12-02 Fix problem of G4DynamicParticle constructor (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 24-03-05 Add G4EmCorrections (V.Ivanchenko)
// 11-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 11-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 12-02-06 move G4LossTableManager::Instance()->EmCorrections() 
//          in constructor (mma)
// 12-08-08 Added methods GetParticleCharge, GetChargeSquareRatio, 
//          CorrectionsAlongStep needed for ions(V.Ivanchenko)
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4BetheBlochModel.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4EmCorrections.hh"
#include "G4EmParameters.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ICRU90StoppingData.hh"
#include "G4Log.hh"
#include "G4DeltaAngle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4BetheBlochModel::G4BetheBlochModel(const G4ParticleDefinition*, 
                                     const G4String& nam)
  : G4VEmModel(nam),
    twoln10(2.0*G4Log(10.0)),
    fAlphaTlimit(1*CLHEP::GeV),
    fProtonTlimit(10*CLHEP::GeV)
{
  theElectron = G4Electron::Electron();
  corr = G4LossTableManager::Instance()->EmCorrections();  
  nist = G4NistManager::Instance();
  SetLowEnergyLimit(2.0*CLHEP::MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BetheBlochModel::~G4BetheBlochModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BetheBlochModel::Initialise(const G4ParticleDefinition* p,
                                   const G4DataVector&)
{
  if(p != particle) { SetupParameters(p); }

  //G4cout << "G4BetheBlochModel::Initialise for " << p->GetParticleName()
  //         << "  isIon= " << isIon 
  //         << G4endl;

  // always false before the run
  SetDeexcitationFlag(false);

  // initialisation once
  if(nullptr == fParticleChange) {
    const G4String& pname = particle->GetParticleName();
    if(IsMaster() && G4EmParameters::Instance()->UseICRU90Data() &&
       (pname == "proton" || pname == "GenericIon" || pname == "alpha")) {
      fICRU90 = nist->GetICRU90StoppingData();
      fICRU90->Initialise();
    }
    if(particle->GetPDGCharge() > CLHEP::eplus ||
       pname == "GenericIon") { isIon = true; } 
    if(pname == "alpha") { isAlpha = true; } 

    fParticleChange = GetParticleChangeForLoss();
    if(UseAngularGeneratorFlag() && nullptr == GetAngularDistribution()) {
      SetAngularDistribution(new G4DeltaAngle());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::GetChargeSquareRatio(const G4ParticleDefinition* p,
                                                 const G4Material* mat,
                                                 G4double kineticEnergy)
{
  // this method is called only for ions, so no check if it is an ion
  return 
    (!isAlpha) ? corr->EffectiveChargeSquareRatio(p,mat,kineticEnergy) : 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::GetParticleCharge(const G4ParticleDefinition* p,
                                              const G4Material* mat,
                                              G4double kineticEnergy)
{
  // this method is called only for ions, so no check if it is an ion
  return corr->GetParticleCharge(p,mat,kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BetheBlochModel::SetupParameters(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  G4double q = particle->GetPDGCharge()*inveplus;
  isIon = (!isAlpha && q > 1.1); 
  chargeSquare = q*q;
  ratio = electron_mass_c2/mass;
  static const G4double aMag = 1./(0.5*eplus*CLHEP::hbar_Planck*CLHEP::c_squared);
  G4double magmom = particle->GetPDGMagneticMoment()*mass*aMag;
  magMoment2 = magmom*magmom - 1.0;
  formfact = 0.0;
  tlimit = DBL_MAX;
  if(particle->GetLeptonNumber() == 0) {
    G4double x = 0.8426*CLHEP::GeV;
    if(spin == 0.0 && mass < CLHEP::GeV) { x = 0.736*CLHEP::GeV; }
    else if (mass > CLHEP::GeV) {
      G4int iz = G4lrint(std::abs(q));
      if(iz > 1) { x /= nist->GetA27(iz); }  
    }
    formfact = 2.0*CLHEP::electron_mass_c2/(x*x);
    tlimit = 2.0/formfact;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::MinEnergyCut(const G4ParticleDefinition*,
                                         const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4BetheBlochModel::ComputeCrossSectionPerElectron(const G4ParticleDefinition* p,
                                                  G4double kineticEnergy,
                                                  G4double cutEnergy,
                                                  G4double maxKinEnergy)        
{
  G4double cross = 0.0;
  G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  G4double maxEnergy = std::min(tmax, maxKinEnergy);
  if(cutEnergy < maxEnergy) {

    G4double totEnergy = kineticEnergy + mass;
    G4double energy2   = totEnergy*totEnergy;
    G4double beta2     = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;

    cross = (maxEnergy - cutEnergy)/(cutEnergy*maxEnergy) 
      - beta2*G4Log(maxEnergy/cutEnergy)/tmax;

    // +term for spin=1/2 particle
    if( 0.0 < spin ) { cross += 0.5*(maxEnergy - cutEnergy)/energy2; }

    cross *= CLHEP::twopi_mc2_rcl2*chargeSquare/beta2;
  }
  
   // G4cout << "BB: e= " << kineticEnergy << " tmin= " << cutEnergy 
   //        << " tmax= " << tmax << " cross= " << cross << G4endl;
  
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition* p,
                                                 G4double kinEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  return Z*ComputeCrossSectionPerElectron(p,kinEnergy,cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::CrossSectionPerVolume(
                                           const G4Material* mat,
                                           const G4ParticleDefinition* p,
                                                 G4double kinEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double sigma = mat->GetElectronDensity() 
    *ComputeCrossSectionPerElectron(p,kinEnergy,cutEnergy,maxEnergy);
  if(isAlpha) {
    sigma *= corr->EffectiveChargeSquareRatio(p,mat,kinEnergy)/chargeSquare;
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::ComputeDEDXPerVolume(const G4Material* material,
                                                 const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cut)
{
  const G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  // projectile formfactor limit energy loss 
  const G4double cutEnergy = std::min(std::min(cut,tmax), tlimit);

  G4double tau   = kineticEnergy/mass;
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau+2.0);
  G4double beta2 = bg2/(gam*gam);
  G4double xc    = cutEnergy/tmax;

  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc2 = eexc*eexc;

  G4double eDensity = material->GetElectronDensity();

  // added ICRU90 stopping data for limited list of materials
  /*
  G4cout << "### DEDX ICRI90:" << (nullptr != fICRU90) 
	 << " Ekin=" << kineticEnergy 
	 << "  " << p->GetParticleName() 
	 << " q2=" << chargeSquare
	 << " inside  " << material->GetName() << G4endl;
  */
  if(nullptr != fICRU90 && kineticEnergy < fProtonTlimit) {
    if(material != currentMaterial) {
      currentMaterial = material;
      baseMaterial = material->GetBaseMaterial() 
        ? material->GetBaseMaterial() : material;
      iICRU90 = fICRU90->GetIndex(baseMaterial);
    }
    if(iICRU90 >= 0) {
      G4double dedx = 0.0;
      // only for alpha
      if(isAlpha) {
	if(kineticEnergy <= fAlphaTlimit) {
	  dedx = fICRU90->GetElectronicDEDXforAlpha(iICRU90, kineticEnergy);
	} else {
          const G4double e = kineticEnergy*CLHEP::proton_mass_c2/mass;
	  dedx = fICRU90->GetElectronicDEDXforProton(iICRU90, e)*chargeSquare;
	}
      } else {
        dedx = fICRU90->GetElectronicDEDXforProton(iICRU90, kineticEnergy)
	  *chargeSquare;
      }
      dedx *= material->GetDensity();
      if(cutEnergy < tmax) {
        dedx += (G4Log(xc) + (1.0 - xc)*beta2)*CLHEP::twopi_mc2_rcl2
          *(eDensity*chargeSquare/beta2);
      }
      //G4cout << "   iICRU90=" << iICRU90 << "   dedx=" << dedx << G4endl;
      if(dedx > 0.0) { return dedx; }
    }
  }
  // general Bethe-Bloch formula
  G4double dedx = G4Log(2.0*CLHEP::electron_mass_c2*bg2*cutEnergy/eexc2)
                - (1.0 + xc)*beta2;

  if(0.0 < spin) {
    G4double del = 0.5*cutEnergy/(kineticEnergy + mass);
    dedx += del*del;
  }

  // density correction
  G4double x = G4Log(bg2)/twoln10;
  dedx -= material->GetIonisation()->DensityCorrection(x);

  // shell correction
  dedx -= 2.0*corr->ShellCorrection(p,material,kineticEnergy);

  // now compute the total ionization loss
  dedx *= CLHEP::twopi_mc2_rcl2*chargeSquare*eDensity/beta2;

  //High order correction different for hadrons and ions
  if(isIon) {
    dedx += corr->IonBarkasCorrection(p,material,kineticEnergy);
  } else {      
    dedx += corr->HighOrderCorrections(p,material,kineticEnergy,cutEnergy);
  }

  dedx = std::max(dedx, 0.0); 
  /*
  G4cout << "E(MeV)= " << kineticEnergy/CLHEP::MeV << " dedx= " << dedx 
           << "  " << material->GetName() << G4endl;
  */
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BetheBlochModel::CorrectionsAlongStep(const G4MaterialCutsCouple* couple,
                                             const G4DynamicParticle* dp,
                                             const G4double& /*length*/,
                                             G4double& eloss)
{
  // no correction for alpha
  if(isAlpha) { return; }

  // no correction at the last step or at small step
  const G4double preKinEnergy = dp->GetKineticEnergy();
  if(eloss >= preKinEnergy || eloss < preKinEnergy*0.05) { return; }

  // corrections for all charged particles with Q > 1
  const G4ParticleDefinition* p = dp->GetDefinition();
  if(p != particle) { SetupParameters(p); }
  if(!isIon) { return; }

  // effective energy and charge at a step
  const G4double e = std::max(preKinEnergy - eloss*0.5, preKinEnergy*0.5);
  const G4Material* mat = couple->GetMaterial();
  const G4double q20 = corr->EffectiveChargeSquareRatio(p, mat, preKinEnergy);
  const G4double q2 = corr->EffectiveChargeSquareRatio(p, mat, e);
  const G4double qfactor = q2/q20;

  /*    
    G4cout << "G4BetheBlochModel::CorrectionsAlongStep: Epre(MeV)="
    << preKinEnergy << " Eeff(MeV)=" << e
    << " eloss=" << eloss << " elossnew=" << eloss*qfactor 
    << " qfactor=" << qfactor << " Qpre=" << q20 
    << p->GetParticleName() <<G4endl;
  */
  eloss *= qfactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BetheBlochModel::SampleSecondaries(vector<G4DynamicParticle*>* vdp,
                                          const G4MaterialCutsCouple* couple,
                                          const G4DynamicParticle* dp,
                                          G4double minKinEnergy,
                                          G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  const G4double tmax = MaxSecondaryEnergy(dp->GetDefinition(),kineticEnergy);
  const G4double maxKinEnergy = std::min(maxEnergy,tmax);
  if(minKinEnergy >= maxKinEnergy) { return; }

  //G4cout << "G4BetheBlochModel::SampleSecondaries Emin= " << minKinEnergy
  //         << " Emax= " << maxKinEnergy << G4endl;

  const G4double totEnergy = kineticEnergy + mass;
  const G4double etot2     = totEnergy*totEnergy;
  const G4double beta2     = kineticEnergy*(kineticEnergy + 2.0*mass)/etot2;

  G4double deltaKinEnergy, f; 
  G4double f1 = 0.0;
  G4double fmax = 1.0;
  if( 0.0 < spin ) { fmax += 0.5*maxKinEnergy*maxKinEnergy/etot2; }

  CLHEP::HepRandomEngine* rndmEngineMod = G4Random::getTheEngine();
  G4double rndm[2];

  // sampling without nuclear size effect
  do {
    rndmEngineMod->flatArray(2, rndm);
    deltaKinEnergy = minKinEnergy*maxKinEnergy
                    /(minKinEnergy*(1.0 - rndm[0]) + maxKinEnergy*rndm[0]);

    f = 1.0 - beta2*deltaKinEnergy/tmax;
    if( 0.0 < spin ) {
      f1 = 0.5*deltaKinEnergy*deltaKinEnergy/etot2;
      f += f1;
    }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while( fmax*rndm[1] > f);

  // projectile formfactor - suppresion of high energy
  // delta-electron production at high energy
  
  G4double x = formfact*deltaKinEnergy;
  if(x > 1.e-6) {

    G4double x1 = 1.0 + x;
    G4double grej  = 1.0/(x1*x1);
    if( 0.0 < spin ) {
      G4double x2 = 0.5*electron_mass_c2*deltaKinEnergy/(mass*mass);
      grej *= (1.0 + magMoment2*(x2 - f1/f)/(1.0 + x2));
    }
    if(grej > 1.1) {
      G4cout << "### G4BetheBlochModel WARNING: grej= " << grej
             << "  " << dp->GetDefinition()->GetParticleName()
             << " Ekin(MeV)= " <<  kineticEnergy
             << " delEkin(MeV)= " << deltaKinEnergy
             << G4endl;
    }
    if(rndmEngineMod->flat() > grej) { return; }
  }

  G4ThreeVector deltaDirection;

  if(UseAngularGeneratorFlag()) {
    const G4Material* mat = couple->GetMaterial();
    deltaDirection = 
      GetAngularDistribution()->SampleDirection(dp, deltaKinEnergy,
						SelectRandomAtomNumber(mat),
						mat);
  } else {
 
    G4double deltaMomentum =
      sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
    G4double cost = deltaKinEnergy * (totEnergy + electron_mass_c2) /
      (deltaMomentum * dp->GetTotalMomentum());
    cost = std::min(cost, 1.0);
    const G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
    const G4double phi = twopi*rndmEngineMod->flat();

    deltaDirection.set(sint*std::cos(phi),sint*std::sin(phi), cost) ;
    deltaDirection.rotateUz(dp->GetMomentumDirection());
  }  
  /*
    G4cout << "### G4BetheBlochModel " 
           << dp->GetDefinition()->GetParticleName()
           << " Ekin(MeV)= " <<  kineticEnergy
           << " delEkin(MeV)= " << deltaKinEnergy
           << " tmin(MeV)= " << minKinEnergy
           << " tmax(MeV)= " << maxKinEnergy
           << " dir= " << dp->GetMomentumDirection()
           << " dirDelta= " << deltaDirection
           << G4endl;
  */
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BetheBlochModel::MaxSecondaryEnergy(const G4ParticleDefinition* pd,
                                               G4double kinEnergy) 
{
  // here particle type is checked for the case, 
  // when this model is shared between particles
  if(pd != particle) { SetupParameters(pd); }
  G4double tau  = kinEnergy/mass;
  return 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.) /
    (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

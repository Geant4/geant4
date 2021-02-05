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
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4EmCorrections.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ICRU90StoppingData.hh"
#include "G4Log.hh"
#include "G4DeltaAngle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4BetheBlochModel::G4BetheBlochModel(const G4ParticleDefinition*, 
                                     const G4String& nam)
  : G4VEmModel(nam),
    particle(nullptr),
    fICRU90(nullptr),
    currentMaterial(nullptr),
    baseMaterial(nullptr),
    tlimit(DBL_MAX),
    twoln10(2.0*G4Log(10.0)),
    fAlphaTlimit(CLHEP::GeV),
    fProtonTlimit(10*CLHEP::GeV),
    iICRU90(-1),
    isIon(false)
{
  fParticleChange = nullptr;
  theElectron = G4Electron::Electron();
  SetParticle(theElectron);
  corr = G4LossTableManager::Instance()->EmCorrections();  
  nist = G4NistManager::Instance();
  SetLowEnergyLimit(2.0*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BetheBlochModel::~G4BetheBlochModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BetheBlochModel::Initialise(const G4ParticleDefinition* p,
                                   const G4DataVector&)
{
  SetGenericIon(p);
  SetParticle(p);

  //G4cout << "G4BetheBlochModel::Initialise for " << p->GetParticleName()
  //         << "  isIon= " << isIon 
  //         << G4endl;

  // always false before the run
  SetDeexcitationFlag(false);

  if(IsMaster() && G4EmParameters::Instance()->UseICRU90Data()) {
    if(!fICRU90) { fICRU90 = nist->GetICRU90StoppingData(); } 
    else if(particle->GetPDGMass() < GeV) { fICRU90->Initialise(); }
  }

  if(nullptr == fParticleChange) {
    fParticleChange = GetParticleChangeForLoss();
    if(UseAngularGeneratorFlag() && !GetAngularDistribution()) {
      SetAngularDistribution(new G4DeltaAngle());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::GetChargeSquareRatio(const G4ParticleDefinition* p,
                                                 const G4Material* mat,
                                                 G4double kineticEnergy)
{
  // this method is called only for ions
  G4double q2 = corr->EffectiveChargeSquareRatio(p,mat,kineticEnergy);
  corrFactor = q2*corr->EffectiveChargeCorrection(p,mat,kineticEnergy);
  return corrFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::GetParticleCharge(const G4ParticleDefinition* p,
                                              const G4Material* mat,
                                              G4double kineticEnergy)
{
  //G4cout<<"G4BetheBlochModel::GetParticleCharge e= "<<kineticEnergy <<
  //  " q= " <<  corr->GetParticleCharge(p,mat,kineticEnergy) <<G4endl;
  // this method is called only for ions, so no check if it is an ion
  return corr->GetParticleCharge(p,mat,kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BetheBlochModel::SetupParameters()
{
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  G4double q = particle->GetPDGCharge()*inveplus;
  chargeSquare = q*q;
  corrFactor = chargeSquare;
  ratio = electron_mass_c2/mass;
  static const G4double aMag = 1./(0.5*eplus*hbar_Planck*c_squared);
  G4double magmom = particle->GetPDGMagneticMoment()*mass*aMag;
  magMoment2 = magmom*magmom - 1.0;
  formfact = 0.0;
  tlimit = DBL_MAX;
  if(particle->GetLeptonNumber() == 0) {
    G4double x = 0.8426*CLHEP::GeV;
    if(spin == 0.0 && mass < GeV) { x = 0.736*CLHEP::GeV; }
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
  G4double maxEnergy = std::min(tmax,maxKinEnergy);
  if(cutEnergy < maxEnergy) {

    G4double totEnergy = kineticEnergy + mass;
    G4double energy2   = totEnergy*totEnergy;
    G4double beta2     = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;

    cross = (maxEnergy - cutEnergy)/(cutEnergy*maxEnergy) 
      - beta2*G4Log(maxEnergy/cutEnergy)/tmax;

    // +term for spin=1/2 particle
    if( 0.0 < spin ) { cross += 0.5*(maxEnergy - cutEnergy)/energy2; }

    cross *= twopi_mc2_rcl2*chargeSquare/beta2;
  }
  
   // G4cout << "BB: e= " << kineticEnergy << " tmin= " << cutEnergy 
   //        << " tmax= " << tmax << " cross= " << cross << G4endl;
  
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  return
    Z*ComputeCrossSectionPerElectron(p,kineticEnergy,cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::CrossSectionPerVolume(
                                           const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  return material->GetElectronDensity()
    *ComputeCrossSectionPerElectron(p,kineticEnergy,cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BetheBlochModel::ComputeDEDXPerVolume(const G4Material* material,
                                                 const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cut)
{
  G4double tmax      = MaxSecondaryEnergy(p, kineticEnergy);
  // projectile formfactor limit energy loss 
  G4double cutEnergy = std::min(std::min(cut,tmax), tlimit);

  G4double tau   = kineticEnergy/mass;
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau+2.0);
  G4double beta2 = bg2/(gam*gam);
  G4double xc    = cutEnergy/tmax;

  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc2 = eexc*eexc;

  G4double eDensity = material->GetElectronDensity();

  // added ICRU90 stopping data for limited list of materials
  if(fICRU90) {
    if(material != currentMaterial) {
      currentMaterial = material;
      baseMaterial = material->GetBaseMaterial() 
        ? material->GetBaseMaterial() : material;
      iICRU90 = fICRU90->GetIndex(baseMaterial);
    }
    if(iICRU90 >= 0) {
      G4double e = kineticEnergy*proton_mass_c2/mass;
      G4double dedx = 0.0;
      if(chargeSquare > 1.1 && e < fAlphaTlimit) {
        dedx = fICRU90->GetElectronicDEDXforAlpha(iICRU90, e)
          *material->GetDensity()*0.25;
      } else if(chargeSquare < 1.1 && e < fProtonTlimit) {
        dedx = fICRU90->GetElectronicDEDXforProton(iICRU90, e)
          *material->GetDensity();
      }
      if(cutEnergy < tmax) {
        dedx += (G4Log(xc) + (1.0 - xc)*beta2)*twopi_mc2_rcl2
          *eDensity/beta2;
        return std::max(chargeSquare*dedx, 0.0);
      }
    } 
  }
  // general Bethe-Bloch formula
  G4double dedx = G4Log(2.0*electron_mass_c2*bg2*cutEnergy/eexc2)
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
  dedx *= twopi_mc2_rcl2*chargeSquare*eDensity/beta2;

  //High order correction different for hadrons and ions
  if(isIon) {
    dedx += corr->IonBarkasCorrection(p,material,kineticEnergy);
  } else {      
    dedx += corr->HighOrderCorrections(p,material,kineticEnergy,cutEnergy);
  }

  dedx = std::max(dedx, 0.0); 

  //G4cout << "E(MeV)= " << kineticEnergy/MeV << " dedx= " << dedx 
  //         << "  " << material->GetName() << G4endl;

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BetheBlochModel::CorrectionsAlongStep(const G4MaterialCutsCouple* couple,
                                             const G4DynamicParticle* dp,
                                             G4double& eloss,
                                             G4double&,
                                             G4double length)
{
  if(isIon) {
    const G4Material* mat = couple->GetMaterial();
    const G4ParticleDefinition* p = dp->GetDefinition();
    G4double preKinEnergy = dp->GetKineticEnergy();
    G4double e = preKinEnergy - eloss*0.5;
    if(e < preKinEnergy*0.75) { e = preKinEnergy*0.75; }

    G4double q2 = corr->EffectiveChargeSquareRatio(p,mat,e);
    GetModelOfFluctuations()->SetParticleAndCharge(p, q2);
    G4double qfactor = q2*corr->EffectiveChargeCorrection(p,mat,e)/corrFactor;

    // no high order correction for ICRU90 data
    baseMaterial = mat->GetBaseMaterial() ? mat->GetBaseMaterial() : mat;
    G4double highOrder = 0.0;
    if(!fICRU90 || fICRU90->GetIndex(baseMaterial) < 0) {
      highOrder = length*corr->IonHighOrderCorrections(p,couple,e);
    }
    G4double elossnew  = eloss*qfactor + highOrder;
    eloss = std::max(std::min(elossnew,preKinEnergy),eloss*0.5);
    //G4cout << "G4BetheBlochModel::CorrectionsAlongStep: e= " << preKinEnergy
    //           << " qfactor= " << qfactor 
    //           << " highOrder= " << highOrder << " (" 
    // << highOrder/eloss << ")" << G4endl;    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BetheBlochModel::SampleSecondaries(vector<G4DynamicParticle*>* vdp,
                                          const G4MaterialCutsCouple* couple,
                                          const G4DynamicParticle* dp,
                                          G4double minKinEnergy,
                                          G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double tmax = MaxSecondaryEnergy(dp->GetDefinition(),kineticEnergy);

  G4double maxKinEnergy = std::min(maxEnergy,tmax);
  if(minKinEnergy >= maxKinEnergy) { return; }

  //G4cout << "G4BetheBlochModel::SampleSecondaries Emin= " << minKinEnergy
  //         << " Emax= " << maxKinEnergy << G4endl;

  G4double totEnergy     = kineticEnergy + mass;
  G4double etot2         = totEnergy*totEnergy;
  G4double beta2         = kineticEnergy*(kineticEnergy + 2.0*mass)/etot2;

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

    const G4Material* mat =  couple->GetMaterial();
    G4int Z = SelectRandomAtomNumber(mat);

    deltaDirection = 
      GetAngularDistribution()->SampleDirection(dp, deltaKinEnergy, Z, mat);

  } else {
 
    G4double deltaMomentum =
      sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
    G4double cost = deltaKinEnergy * (totEnergy + electron_mass_c2) /
      (deltaMomentum * dp->GetTotalMomentum());
    cost = std::min(cost, 1.0);
    G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));

    G4double phi = twopi*rndmEngineMod->flat();

    deltaDirection.set(sint*cos(phi),sint*sin(phi), cost) ;
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
  G4DynamicParticle* delta = 
    new G4DynamicParticle(theElectron,deltaDirection,deltaKinEnergy);

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
  SetParticle(pd);
  G4double tau  = kinEnergy/mass;
  return 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.) /
    (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

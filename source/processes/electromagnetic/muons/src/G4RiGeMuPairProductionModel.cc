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
// File name:     G4RiGeMuPairProductionModel
//
// Authors:       Girardo Depaola & Ricardo Pacheco
//
// Creation date: 29.10.2024
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RiGeMuPairProductionModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementDataRegistry.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4RiGeAngularGenerator.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4AutoLock.hh"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int G4RiGeMuPairProductionModel::ZDATPAIR[] = {1, 4, 13, 29, 92};

const G4double G4RiGeMuPairProductionModel::xgi[] = {
    0.0198550717512320, 0.1016667612931865, 0.2372337950418355, 0.4082826787521750,
    0.5917173212478250, 0.7627662049581645, 0.8983332387068135, 0.9801449282487680
  };

const G4double G4RiGeMuPairProductionModel::wgi[] = {
    0.0506142681451880, 0.1111905172266870, 0.1568533229389435, 0.1813418916891810,
    0.1813418916891810, 0.1568533229389435, 0.1111905172266870, 0.0506142681451880
  };

namespace
{
  G4Mutex theRiGeMuPairMutex = G4MUTEX_INITIALIZER;

  const G4double ak1 = 6.9;
  const G4double ak2 = 1.0;

  // Channel weights
  const G4double W[3] = {0.25, 0.5, 0.75};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4RiGeMuPairProductionModel::G4RiGeMuPairProductionModel(const G4ParticleDefinition* p)
  : G4VEmModel("muPairProdRiGe"),
    factorForCross(CLHEP::fine_structure_const*CLHEP::fine_structure_const*
		   CLHEP::classic_electr_radius*CLHEP::classic_electr_radius*
		   4./(3.*CLHEP::pi)),
    sqrte(std::sqrt(G4Exp(1.))),
    minPairEnergy(4.*CLHEP::electron_mass_c2),
    lowestKinEnergy(0.85*CLHEP::GeV)
{
  nist = G4NistManager::Instance();

  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();

  if (nullptr != p) { 
    SetParticle(p); 
    lowestKinEnergy = std::max(lowestKinEnergy, p->GetPDGMass()*8.0);  
  }
  emin = lowestKinEnergy;
  emax = emin*10000.;
  fAngularGenerator = new G4RiGeAngularGenerator();
  SetAngularDistribution(fAngularGenerator);
  for (G4int i=0; i<9; ++i) { randNumbs[i] = 0.0; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4RiGeMuPairProductionModel::MinPrimaryEnergy(const G4Material*,
                                              const G4ParticleDefinition*,
                                              G4double cut)
{
  return std::max(lowestKinEnergy, cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RiGeMuPairProductionModel::Initialise(const G4ParticleDefinition* p,
                                             const G4DataVector& cuts)
{ 
  SetParticle(p); 

  if (nullptr == fParticleChange) { 
    fParticleChange = GetParticleChangeForLoss();

    // define scale of internal table for each thread only once
    if (0 == nbine) {
      emin = std::max(lowestKinEnergy, LowEnergyLimit());
      emax = std::max(HighEnergyLimit(), emin*2);
      nbine = std::size_t(nYBinPerDecade*std::log10(emax/emin));
      if(nbine < 3) { nbine = 3; }

      ymin = G4Log(minPairEnergy/emin);
      dy = -ymin/G4double(nbiny);
    }
    if (p == particle) {
      G4int pdg = std::abs(p->GetPDGEncoding());
      if (pdg == 2212) {
        dataName = "pEEPairProd";
      } else if (pdg == 321) {
        dataName = "kaonEEPairProd";
      } else if (pdg == 211) {
        dataName = "pionEEPairProd";
      } else if (pdg == 11) {
        dataName = "eEEPairProd";
      } else if (pdg == 13) {
        if (GetName() == "muToMuonPairProd") {
          dataName = "muMuMuPairProd";
	} else {
	  dataName = "muEEPairProd";
	}
      } 
    }
  }

  // for low-energy application this process should not work
  if(lowestKinEnergy >= HighEnergyLimit()) { return; }

  if (p == particle) {
    auto data = G4ElementDataRegistry::Instance();
    fElementData = data->GetElementDataByName(dataName);
    if (nullptr == fElementData) { 
      G4AutoLock l(&theRiGeMuPairMutex);
      fElementData = data->GetElementDataByName(dataName);
      if (nullptr == fElementData) { 
        fElementData = new G4ElementData(NZDATPAIR);
        fElementData->SetName(dataName);
      }
      G4bool useDataFile = G4EmParameters::Instance()->RetrieveMuDataFromFile();
      if (useDataFile)  { useDataFile = RetrieveTables(); }
      if (!useDataFile) { MakeSamplingTables(); }
      if (fTableToFile) { StoreTables(); }
      l.unlock();
    }
    if (IsMaster()) {
      InitialiseElementSelectors(p, cuts);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4RiGeMuPairProductionModel::InitialiseLocal(const G4ParticleDefinition* p,
                                              G4VEmModel* masterModel)
{
  if(p == particle && lowestKinEnergy < HighEnergyLimit()) {
    SetElementSelectors(masterModel->GetElementSelectors());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4RiGeMuPairProductionModel::ComputeDEDXPerVolume(const G4Material* material,
                                                  const G4ParticleDefinition*,
                                                  G4double kineticEnergy,
                                                  G4double cutEnergy)
{
  G4double dedx = 0.0;
  if (cutEnergy <= minPairEnergy || kineticEnergy <= lowestKinEnergy)
    { return dedx; }

  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector =
                                   material->GetAtomicNumDensityVector();

  //  loop for elements in the material
  for (std::size_t i=0; i<material->GetNumberOfElements(); ++i) {
     G4double Z = (*theElementVector)[i]->GetZ();
     G4double tmax = MaxSecondaryEnergyForElement(kineticEnergy, Z);
     G4double loss = ComputMuPairLoss(Z, kineticEnergy, cutEnergy, tmax);
     dedx += loss*theAtomicNumDensityVector[i];
  }
  dedx = std::max(dedx, 0.0);
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4RiGeMuPairProductionModel::ComputMuPairLoss(G4double Z, G4double tkin,
                                                       G4double cutEnergy, 
                                                       G4double tmax)
{
  G4double loss = 0.0;

  G4double cut = std::min(cutEnergy, tmax);
  if(cut <= minPairEnergy) { return loss; }

  // calculate the rectricted loss
  // numerical integration in log(PairEnergy)
  G4double aaa = G4Log(minPairEnergy);
  G4double bbb = G4Log(cut);

  G4int kkk = std::min(std::max(G4lrint((bbb-aaa)/ak1 + ak2), 8), 1);
  G4double hhh = (bbb-aaa)/kkk;
  G4double x = aaa;

  for (G4int l=0 ; l<kkk; ++l) {
    for (G4int ll=0; ll<NINTPAIR; ++ll) {
      G4double ep = G4Exp(x+xgi[ll]*hhh);
      loss += wgi[ll]*ep*ep*ComputeDMicroscopicCrossSection(tkin, Z, ep);
    }
    x += hhh;
  }
  loss *= hhh;
  loss = std::max(loss, 0.0);
  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4RiGeMuPairProductionModel::ComputeMicroscopicCrossSection(G4double tkin,
							    G4double Z,
							    G4double cutEnergy)
{
  G4double cross = 0.;
  G4double tmax = MaxSecondaryEnergyForElement(tkin, Z);
  G4double cut  = std::max(cutEnergy, minPairEnergy);
  if (tmax <= cut) { return cross; }

  G4double aaa = G4Log(cut);
  G4double bbb = G4Log(tmax);
  G4int kkk = std::min(std::max(G4lrint((bbb-aaa)/ak1 + ak2), 8), 1);

  G4double hhh = (bbb-aaa)/(kkk);
  G4double x = aaa;

  for (G4int l=0; l<kkk; ++l) {
    for (G4int i=0; i<NINTPAIR; ++i) {
      G4double ep = G4Exp(x + xgi[i]*hhh);
      cross += ep*wgi[i]*ComputeDMicroscopicCrossSection(tkin, Z, ep);
    }
    x += hhh;
  }

  cross *= hhh;
  cross = std::max(cross, 0.0);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4RiGeMuPairProductionModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double pairEnergy)
// Calculates the  differential (D) microscopic cross section
// using the cross section formula of R.P. Kokoulin (18/01/98)
// Code modified by R.P. Kokoulin, V.N. Ivanchenko (27/01/04)
{
  static const G4double bbbtf= 183. ;
  static const G4double bbbh = 202.4 ;
  static const G4double g1tf = 1.95e-5 ;
  static const G4double g2tf = 5.3e-5 ;
  static const G4double g1h  = 4.4e-5 ;
  static const G4double g2h  = 4.8e-5 ;

  if (pairEnergy <= minPairEnergy)
    return 0.0;

  G4double totalEnergy  = tkin + particleMass;
  G4double residEnergy  = totalEnergy - pairEnergy;

  if (residEnergy <= 0.75*sqrte*z13*particleMass)
    return 0.0;

  G4double a0 = 1.0 / (totalEnergy * residEnergy);
  G4double alf = 4.0 * electron_mass_c2 / pairEnergy;
  G4double rt = std::sqrt(1.0 - alf);
  G4double delta = 6.0 * particleMass * particleMass * a0;
  G4double tmnexp = alf/(1.0 + rt) + delta*rt;

  if(tmnexp >= 1.0) { return 0.0; }

  G4double tmn = G4Log(tmnexp);

  G4double massratio = particleMass/CLHEP::electron_mass_c2;
  G4double massratio2 = massratio*massratio;
  G4double inv_massratio2 = 1.0 / massratio2;

  // zeta calculation
  G4double bbb,g1,g2;
  if( Z < 1.5 ) { bbb = bbbh ; g1 = g1h ; g2 = g2h ; }
  else          { bbb = bbbtf; g1 = g1tf; g2 = g2tf; }

  G4double zeta = 0.0;
  G4double z1exp = totalEnergy / (particleMass + g1*z23*totalEnergy);

  // 35.221047195922 is the root of zeta1(x) = 0.073 * log(x) - 0.26, so the
  // condition below is the same as zeta1 > 0.0, but without calling log(x)
  if (z1exp > 35.221047195922)
  {
    G4double z2exp = totalEnergy / (particleMass + g2*z13*totalEnergy);
    zeta = (0.073 * G4Log(z1exp) - 0.26) / (0.058 * G4Log(z2exp) - 0.14);
  }

  G4double z2 = Z*(Z+zeta);
  G4double screen0 = 2.*electron_mass_c2*sqrte*bbb/(z13*pairEnergy);
  G4double beta = 0.5*pairEnergy*pairEnergy*a0;
  G4double xi0 = 0.5*massratio2*beta;

  // Gaussian integration in ln(1-ro) ( with 8 points)
  G4double rho[NINTPAIR];
  G4double rho2[NINTPAIR];
  G4double xi[NINTPAIR];
  G4double xi1[NINTPAIR];
  G4double xii[NINTPAIR];

  for (G4int i = 0; i < NINTPAIR; ++i)
  {
    rho[i] = G4Exp(tmn*xgi[i]) - 1.0; // rho = -asymmetry
    rho2[i] = rho[i] * rho[i];
    xi[i] = xi0*(1.0-rho2[i]);
    xi1[i] = 1.0 + xi[i];
    xii[i] = 1.0 / xi[i];
  }

  G4double ye1[NINTPAIR];
  G4double ym1[NINTPAIR];

  G4double b40 = 4.0 * beta;
  G4double b62 = 6.0 * beta + 2.0;

  for (G4int i = 0; i < NINTPAIR; ++i)
  {
    G4double yeu = (b40 + 5.0) + (b40 - 1.0) * rho2[i];
    G4double yed = b62*G4Log(3.0 + xii[i]) + (2.0 * beta - 1.0)*rho2[i] - b40;

    G4double ymu = b62 * (1.0 + rho2[i]) + 6.0;
    G4double ymd = (b40 + 3.0)*(1.0 + rho2[i])*G4Log(3.0 + xi[i])
      + 2.0 - 3.0 * rho2[i];

    ye1[i] = 1.0 + yeu / yed;
    ym1[i] = 1.0 + ymu / ymd;
  }

  G4double be[NINTPAIR];
  G4double bm[NINTPAIR];

  for(G4int i = 0; i < NINTPAIR; ++i) {
    if(xi[i] <= 1000.0) {
      be[i] = ((2.0 + rho2[i])*(1.0 + beta) +
         xi[i]*(3.0 + rho2[i]))*G4Log(1.0 + xii[i]) +
  (1.0 - rho2[i] - beta)/xi1[i] - (3.0 + rho2[i]);
    } else {
      be[i] = 0.5*(3.0 - rho2[i] + 2.0*beta*(1.0 + rho2[i]))*xii[i];
    }

    if(xi[i] >= 0.001) {
      G4double a10 = (1.0 + 2.0 * beta) * (1.0 - rho2[i]);
      bm[i] = ((1.0 + rho2[i])*(1.0 + 1.5 * beta) - a10*xii[i])*G4Log(xi1[i]) +
                xi[i] * (1.0 - rho2[i] - beta)/xi1[i] + a10;
    } else {
      bm[i] = 0.5*(5.0 - rho2[i] + beta * (3.0 + rho2[i]))*xi[i];
    }
  }

  G4double sum = 0.0;

  for (G4int i = 0; i < NINTPAIR; ++i) {
    G4double screen = screen0*xi1[i]/(1.0 - rho2[i]);
    G4double ale = G4Log(bbb/z13*std::sqrt(xi1[i]*ye1[i])/(1. + screen*ye1[i]));
    G4double cre = 0.5*G4Log(1. + 2.25*z23*xi1[i]*ye1[i]*inv_massratio2);

    G4double fe = (ale-cre)*be[i];
    fe = std::max(fe, 0.0);

    G4double alm_crm = G4Log(bbb*massratio/(1.5*z23*(1. + screen*ym1[i])));
    G4double fm = std::max(alm_crm*bm[i], 0.0)*inv_massratio2;

    sum += wgi[i]*(1.0 + rho[i])*(fe + fm);
  }

  return -tmn*sum*factorForCross*z2*residEnergy/(totalEnergy*pairEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4RiGeMuPairProductionModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                        G4double kineticEnergy,
                                                        G4double Z, G4double,
                                                        G4double cutEnergy,
                                                        G4double maxEnergy)
{
  G4double cross = 0.0;
  if (kineticEnergy <= lowestKinEnergy) { return cross; }

  G4double maxPairEnergy = MaxSecondaryEnergyForElement(kineticEnergy, Z);
  G4double tmax = std::min(maxEnergy, maxPairEnergy);
  G4double cut  = std::max(cutEnergy, minPairEnergy);
  if (cut >= tmax) { return cross; }

  cross = ComputeMicroscopicCrossSection(kineticEnergy, Z, cut);
  if(tmax < kineticEnergy) {
    cross -= ComputeMicroscopicCrossSection(kineticEnergy, Z, tmax);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RiGeMuPairProductionModel::MakeSamplingTables()
{
  G4double factore = G4Exp(G4Log(emax/emin)/G4double(nbine));

  for (G4int iz=0; iz<NZDATPAIR; ++iz) {

    G4double Z = ZDATPAIR[iz];
    G4Physics2DVector* pv = new G4Physics2DVector(nbiny+1,nbine+1);
    G4double kinEnergy = emin;

    for (std::size_t it=0; it<=nbine; ++it) {

      pv->PutY(it, G4Log(kinEnergy/CLHEP::MeV));
      G4double maxPairEnergy = MaxSecondaryEnergyForElement(kinEnergy, Z);
      /*
      G4cout << "it= " << it << " E= " << kinEnergy 
             << "  " << particle->GetParticleName()   
             << " maxE= " << maxPairEnergy << "  minE= " << minPairEnergy 
             << " ymin= " << ymin << G4endl;
      */
      G4double coef = G4Log(minPairEnergy/kinEnergy)/ymin;
      G4double ymax = G4Log(maxPairEnergy/kinEnergy)/coef;
      G4double fac  = (ymax - ymin)/dy;
      std::size_t imax   = (std::size_t)fac;
      fac -= (G4double)imax;
   
      G4double xSec = 0.0;
      G4double x = ymin;
      /*
      G4cout << "Z= " << currentZ << " z13= " << z13 
             << " mE= " << maxPairEnergy << "  ymin= " << ymin 
             << " dy= " << dy << "  c= " << coef << G4endl;
      */
      // start from zero
      pv->PutValue(0, it, 0.0);
      if(0 == it) { pv->PutX(nbiny, 0.0); }

      for (std::size_t i=0; i<nbiny; ++i) {

        if(0 == it) { pv->PutX(i, x); }

        if(i < imax) {
          G4double ep = kinEnergy*G4Exp(coef*(x + dy*0.5));

          // not multiplied by interval, because table 
          // will be used only for sampling
          //G4cout << "i= " << i << " x= " << x << "E= " << kinEnergy  
          //         << " Egamma= " << ep << G4endl;
          xSec += ep*ComputeDMicroscopicCrossSection(kinEnergy, Z, ep);

          // last bin before the kinematic limit
        } else if(i == imax) {
          G4double ep = kinEnergy*G4Exp(coef*(x + fac*dy*0.5));
          xSec += ep*fac*ComputeDMicroscopicCrossSection(kinEnergy, Z, ep);
        }
        pv->PutValue(i + 1, it, xSec);
        x += dy;
      } 
      kinEnergy *= factore;

      // to avoid precision lost
      if(it+1 == nbine) { kinEnergy = emax; }
    }
    fElementData->InitialiseForElement(iz, pv);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RiGeMuPairProductionModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp, 
                                                const G4MaterialCutsCouple* couple,
                                                const G4DynamicParticle* aDynamicParticle,
                                                G4double tmin,
                                                G4double tmax)
{ 
  G4double eMass = CLHEP::electron_mass_c2;
  G4double eMass2 = eMass*eMass;
  
  // Energy and momentum of the pramary particle
  G4double kinEnergy = aDynamicParticle->GetKineticEnergy();
  G4double particleMomentum = aDynamicParticle->GetTotalMomentum();
  G4ThreeVector particleMomentumVector = aDynamicParticle->GetMomentum();
  G4ThreeVector partDirection = aDynamicParticle->GetMomentumDirection();

  G4double minQ2 = 4.*eMass2;
  G4double maxQ2 = (kinEnergy - particleMass)*(kinEnergy - particleMass);
  G4double intervalQ2 = G4Log(maxQ2/minQ2);
  
  // Square invariant of mass of the pair
  G4double Q2 = minQ2*G4Exp(intervalQ2*randNumbs[4]);

  G4double mingEnergy = std::sqrt(Q2);
  G4double maxgEnergy = kinEnergy - particleMass;
  G4double intervalgEnergy = maxgEnergy - mingEnergy;

  // Energy of virtual gamma
  G4double gEnergy = mingEnergy + intervalgEnergy*randNumbs[5];

  // Momentum module of the virtual gamma
  G4double gMomentum = std::sqrt(gEnergy*gEnergy - Q2);

  // Energy and momentum module of the outgoing parent particle
  G4double particleFinalEnergy = kinEnergy - gEnergy;
  G4double particleFinalMomentum = std::sqrt(particleFinalEnergy*particleFinalEnergy -
					     particleMass*particleMass);

  G4double mint3 = 0.;
  G4double maxt3 = CLHEP::pi;
  G4double Cmin = std::cos(maxt3);
  G4double Cmax = std::cos(mint3);

  //G4cout << "------- G4RiGeMuPairProductionModel::SampleSecondaries E(MeV)= " 
  //         << kinEnergy << "  " 
  //         << aDynamicParticle->GetDefinition()->GetParticleName() << G4endl;

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(couple,particle,kinEnergy);

  // define interval of energy transfer
  G4double maxPairEnergy = MaxSecondaryEnergyForElement(kinEnergy, 
                                                        anElement->GetZ());
  G4double maxEnergy = std::min(tmax, maxPairEnergy);
  G4double minEnergy = std::max(tmin, minPairEnergy);

  if (minEnergy >= maxEnergy) { return; }

  //G4cout << "emin= " << minEnergy << " emax= " << maxEnergy 
  // << " minPair= " << minPairEnergy << " maxpair= " << maxPairEnergy 
  //    << " ymin= " << ymin << " dy= " << dy << G4endl;

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  
  G4double coeff = G4Log(minPairEnergy/kinEnergy)/ymin;

  // compute limits 
  G4double yymin = G4Log(minEnergy/kinEnergy)/coeff;
  G4double yymax = G4Log(maxEnergy/kinEnergy)/coeff;
 
  //G4cout << "yymin= " << yymin << "  yymax= " << yymax << G4endl;

  // units should not be used, bacause table was built without
  G4double logTkin = G4Log(kinEnergy/CLHEP::MeV);

  // sample e-e+ energy, pair energy first

  // select sample table via Z
  G4int iz1(0), iz2(0);
  for (G4int iz=0; iz<NZDATPAIR; ++iz) { 
    if(currentZ == ZDATPAIR[iz]) {
      iz1 = iz2 = iz; 
      break;
    } else if(currentZ < ZDATPAIR[iz]) {
      iz2 = iz;
      if(iz > 0) { iz1 = iz-1; }
      else { iz1 = iz2; }
      break;
    } 
  }
  if (0 == iz1) { iz1 = iz2 = NZDATPAIR-1; }

  G4double pairEnergy = 0.0;
  G4int count = 0;
  //G4cout << "start loop Z1= " << iz1 << " Z2= " << iz2 << G4endl;
  do {
    ++count;
    // sampling using only one random number
    G4double rand = rndmEngine->flat();
  
    G4double x = FindScaledEnergy(iz1, rand, logTkin, yymin, yymax);
    if(iz1 != iz2) {
      G4double x2 = FindScaledEnergy(iz2, rand, logTkin, yymin, yymax);
      G4double lz1= nist->GetLOGZ(ZDATPAIR[iz1]);
      G4double lz2= nist->GetLOGZ(ZDATPAIR[iz2]);
      //G4cout << count << ".  x= " << x << "  x2= " << x2 
      //             << " Z1= " << iz1 << " Z2= " << iz2 << G4endl;
      x += (x2 - x)*(lnZ - lz1)/(lz2 - lz1);
    }
    //G4cout << "x= " << x << "  coeff= " << coeff << G4endl;
    pairEnergy = kinEnergy*G4Exp(x*coeff);
    
    // Loop checking, 30-Oct-2024, Vladimir Ivanchenko
  } while((pairEnergy < minEnergy || pairEnergy > maxEnergy) && 50 > count);

  //G4cout << "## pairEnergy(GeV)= " << pairEnergy/GeV 
  //         << " Etot(GeV)= " << totalEnergy/GeV << G4endl; 
  rndmEngine->flatArray(9, randNumbs);
  G4double phi3 = CLHEP::twopi*randNumbs[0];
  fAngularGenerator->PhiRotation(partDirection, phi3);

  G4LorentzVector muF;
  G4ThreeVector eDirection, pDirection;
  G4double eEnergy, pEnergy;
  
  if (randNumbs[7] < W[0]) {
    G4double A1 = -(Q2 - 2.*kinEnergy*gEnergy);
    G4double B1 = -(2.*gMomentum*particleMomentum);
    G4double tginterval = G4Log((A1 + B1)/(A1 - B1))/B1;
    
    G4double costg = (-A1 + (A1 - B1)*G4Exp(B1*tginterval*randNumbs[1]))/B1;
    G4double sintg = std::sqrt((1.0 - costg)*(1.0 + costg));
    G4double phig  = CLHEP::twopi*randNumbs[2];
    G4double sinpg = std::sin(phig);
    G4double cospg = std::cos(phig);  

    G4ThreeVector dirGamma;
    dirGamma.set(sintg*cospg, sintg*sinpg, costg);
    G4LorentzVector gFourMomentum(gEnergy, dirGamma*gMomentum);

    G4double Ap = particleMomentum*particleMomentum +
      particleFinalMomentum*particleFinalMomentum + gMomentum*gMomentum;
    G4double A = Ap - 2.*particleMomentum*gMomentum*costg;
    G4double B = 2.*particleMomentum*gMomentum*sintg*cospg;
    G4double C = 2.*particleFinalMomentum*gMomentum*costg -
      2.*particleMomentum*particleFinalMomentum;
    G4double absB = std::abs(B);
    G4double t1interval = (1./(A + C + absB*mint3) - 1./(A + C + absB*maxt3))/absB;
    G4double t1 = (-(A + C) + 1./(1./(A + C + absB*mint3) - absB*t1interval*randNumbs[0]))/absB;
    G4double sint1 = std::sin(t1);
    G4double cost1 = std::cos(t1);

    // Ingoing parent particle change
    G4double Phi = CLHEP::twopi*randNumbs[3];
    partDirection.set(sint1, 0., cost1);
    fAngularGenerator->PhiRotation(partDirection, Phi);
    kinEnergy = particleFinalEnergy;

    G4double cost5 = -1. + 2.*randNumbs[6];
    G4double phi5 = CLHEP::twopi*randNumbs[8];

    G4LorentzVector eFourMomentumMQ = fAngularGenerator->eDP2(Q2, eMass2, eMass2, cost5, phi5);
    G4LorentzVector pFourMomentumMQ = fAngularGenerator->pDP2(eMass2, eFourMomentumMQ);

    G4LorentzVector eFourMomentum = eFourMomentumMQ.boost(gFourMomentum.boostVector());
    G4LorentzVector pFourMomentum = pFourMomentumMQ.boost(gFourMomentum.boostVector());

    eEnergy = eFourMomentum.t();
    pEnergy = pFourMomentum.t();

  } else if (randNumbs[7] >= W[0] && randNumbs[7] < W[1]) {
    G4double A3 = Q2 + 2.*gEnergy*particleFinalEnergy;
    G4double B3 = -2.*gMomentum*particleFinalMomentum;
    
    G4double tQ3interval = G4Log((A3 + B3)/(A3 - B3))/B3;
    G4double tQMG = (-A3 + (A3 - B3)*G4Exp(B3*tQ3interval*randNumbs[0]))/B3;
    G4double phiQP = CLHEP::twopi*randNumbs[2];
    
    G4double sintQ3 = std::sqrt(1. - tQMG*tQMG);
    G4double cospQP = std::cos(phiQP);
    G4double sinpQP = std::sin(phiQP);
    
    G4double Ap = particleMomentum*particleMomentum +
      particleFinalMomentum*particleFinalMomentum + gMomentum*gMomentum;
    G4double A = Ap + 2.*particleFinalMomentum*gMomentum*tQMG;
    G4double B = -2.*particleMomentum*gMomentum*sintQ3*cospQP;
    G4double C = -2.*particleMomentum*gMomentum*tQMG - 2.*particleMomentum*particleFinalMomentum; 

    G4double absB = std::abs(B);
    G4double t3interval = (1./(A + C + absB*mint3) - 1./(A + C + absB*maxt3))/absB;
    G4double t3 = (-(A + C) + 1./(1./(A + C + absB*mint3) - absB*t3interval*randNumbs[0]))/absB;
    G4double sint3 = std::sin(t3);
    G4double cost3 = std::cos(t3);

    G4double cost = -sint3*sintQ3*cospQP + cost3*tQMG;
    G4double sint = std::sqrt((1. + cost)*(1. - cost));
    G4double cosp = (sintQ3*cospQP*cost3 + sint3*tQMG)/sint;
    G4double sinp = sintQ3*sinpQP/sint;
    
    G4ThreeVector dirGamma;
    dirGamma.set(sint*cosp, sint*sinp, cost);
    G4LorentzVector gFourMomentum(gEnergy, dirGamma*gMomentum);

    // Ingoing parent particle change
    G4double Phi = CLHEP::twopi*randNumbs[3];
    partDirection.set(sint3, 0., cost3);
    fAngularGenerator->PhiRotation(partDirection, Phi);
    kinEnergy = particleFinalEnergy;

    G4double cost5 = -1. + 2.*randNumbs[6];
    G4double phi5 = CLHEP::twopi*randNumbs[8];

    G4LorentzVector eFourMomentumMQ = fAngularGenerator->eDP2(Q2, eMass2, eMass2, cost5, phi5);
    G4LorentzVector pFourMomentumMQ = fAngularGenerator->pDP2(eMass2, eFourMomentumMQ);

    G4LorentzVector eFourMomentum = eFourMomentumMQ.boost(gFourMomentum.boostVector());
    G4LorentzVector pFourMomentum = pFourMomentumMQ.boost(gFourMomentum.boostVector());

    eEnergy = eFourMomentum.t();
    pEnergy = pFourMomentum.t();

  } else if (randNumbs[7] >= W[1] && randNumbs[7] < W[2]) {
    G4double phi5 = CLHEP::twopi*randNumbs[1];
    G4double phi6 = CLHEP::twopi*randNumbs[2];
    G4double muEnergyInterval = kinEnergy - 2.*eMass - particleMass;
    particleFinalEnergy = particleMass + muEnergyInterval*randNumbs[3];
    particleFinalMomentum = std::sqrt(particleFinalEnergy*particleFinalEnergy -
				      particleMass*particleMass);
    
    G4double mineEnergy = eMass;
    G4double maxeEnergy = kinEnergy - particleFinalEnergy - eMass;
    G4double eEnergyinterval = maxeEnergy - mineEnergy;
    eEnergy = mineEnergy + eEnergyinterval*randNumbs[4];
    
    G4double cosp3 = 1.;
    G4double sinp3 = 0.;
    G4double cosp5 = std::cos(phi5);
    G4double sinp5 = std::sin(phi5);
    G4double cosp6 = std::cos(phi6);
    G4double sinp6 = std::sin(phi6);

    G4double eMomentum = std::sqrt(eEnergy*eEnergy - eMass*eMass);
    pEnergy = kinEnergy - particleFinalEnergy - eEnergy;
    G4double pMomentum = std::sqrt(pEnergy*pEnergy - eMass*eMass);
    
    G4double A3 = -2.*particleMass*particleMass + 2.*kinEnergy*particleFinalEnergy;
    G4double B3 = -2.*particleMomentum*particleFinalMomentum;
    G4double cost3interval = G4Log((A3 + B3*Cmax)/(A3 + B3*Cmin))/B3;
    G4double expanCost3r6 = G4Exp(B3*cost3interval*randNumbs[5]);
    G4double cost3 = A3*(expanCost3r6 - 1.)/B3 + Cmin*expanCost3r6;
    G4double sint3 = std::sqrt((1. - cost3)*(1. + cost3));

    partDirection.set(sint3, 0., cost3);
    
    G4ThreeVector muFinalMomentumVector;
    muFinalMomentumVector.set(particleFinalMomentum*sint3, 0., particleFinalMomentum*cost3);

    G4LorentzVector muFourMomentum(particleMomentum, particleMomentumVector);
    G4LorentzVector muFinalFourMomentum(particleFinalEnergy, muFinalMomentumVector);
    G4LorentzVector auxVec1 = muFourMomentum - muFinalFourMomentum;
    G4double A5 = auxVec1.mag2() - 2.*eEnergy*(kinEnergy - particleFinalEnergy) +
      2.*particleMomentumVector[2]*eMomentum - 2.*particleFinalMomentum*eMomentum*cost3;
    G4double B5 = -2.*particleFinalMomentum*eMomentum*(sint3*cosp3*cosp5 + sint3*sinp3*sinp5);
    G4double absA5 = std::abs(A5);
    G4double absB5 = std::abs(B5);
    G4double mint5 = 0.;
    G4double maxt5 = CLHEP::pi;
    G4double t5interval = G4Log((absA5 + absB5*maxt5)/(absA5 + absB5*mint5))/absB5;
    G4double argexp = absB5*t5interval*randNumbs[6] + G4Log(absA5 + absB5*mint5);
    G4double t5 = -absA5/absB5 + G4Exp(argexp)/absB5;
    G4double sint5 = std::sin(t5);
    G4double cost5 = std::cos(t5);

    eDirection.set(sint5*cosp5, sint5*sinp5, cost5);
    G4ThreeVector eMomentumVector = eMomentum*eDirection;

    G4ThreeVector auxVec2 = particleMomentumVector - muFinalMomentumVector - eMomentumVector;
    G4double p1mp3mp52 = auxVec2.dot(auxVec2);
    G4double Bp = particleFinalMomentum*(sint3*cosp3*cosp6 + sint3*sinp3*sinp6) +
      eMomentum*(sint5*cosp5*cosp6 + sint5*sinp5*sinp6);
    G4double Cp = -particleMomentum + particleFinalMomentum*cost3 + eMomentum*cost5;
    G4double A6 = p1mp3mp52 + pMomentum*pMomentum;
    G4double B6 = 2.*pMomentum*Bp;
    G4double C6 = 2.*pMomentum*Cp;
    G4double mint6 = 0.;
    G4double maxt6 = CLHEP::pi;
    G4double absA6C6 = std::abs(A6 + C6);
    G4double absB6 = std::abs(B6);
    G4double t6interval = (1./(absA6C6 + absB6*mint6) - 1./(absA6C6 + absB6*maxt6))/absB6;
    G4double t6 = (-absA6C6 + 1./(1./(absA6C6 + absB6*mint6) - absB6*t6interval*randNumbs[8]))/absB6;
    G4double sint6 = std::sin(t6);
    G4double cost6 = std::cos(t6);
    
    pDirection.set(sint6*cosp6, sint6*sinp6, cost6);

  } else {
    G4double phi6 = CLHEP::twopi*randNumbs[1];
    G4double phi5 = CLHEP::twopi*randNumbs[2];
    G4double muFinalEnergyinterval = kinEnergy - 2.*eMass - particleMass;
    particleFinalEnergy = particleMass + muFinalEnergyinterval*randNumbs[3];
    particleFinalMomentum = std::sqrt(particleFinalEnergy*particleFinalEnergy -
				      particleMass*particleMass);
    
    G4double maxpEnergy = kinEnergy - particleFinalEnergy - eMass;
    G4double pEnergyinterval = maxpEnergy - eMass;
    pEnergy = eMass + pEnergyinterval*randNumbs[4];
    
    G4double cosp3 = 1.;
    G4double sinp3 = 0.;
    G4double cosp5 = std::cos(phi5);
    G4double sinp5 = std::sin(phi5);
    G4double cosp6 = std::cos(phi6);
    G4double sinp6 = std::sin(phi6);

    G4double pMomentum = std::sqrt(pEnergy*pEnergy - eMass*eMass);
    eEnergy = kinEnergy - particleFinalEnergy - pEnergy;
    G4double eMomentum = std::sqrt(eEnergy*eEnergy - eMass*eMass);
    
    G4double A3 = -2.*particleMass*particleMass + 2.*kinEnergy*particleFinalEnergy;
    G4double B3 = -2.*particleMomentum*particleFinalMomentum;
    G4double cost3interval = G4Log((A3 + B3*Cmax)/(A3 + B3*Cmin))/B3;
    G4double expanCost3r6 = G4Exp(B3*cost3interval*randNumbs[5]);
    G4double cost3 = A3*(expanCost3r6 - 1.)/B3 + Cmin*expanCost3r6;
    G4double sint3 = std::sqrt((1. - cost3)*(1. + cost3));

    partDirection.set(sint3*cosp3, sint3*sinp3, cost3);

    G4ThreeVector muFinalMomentumVector;
    muFinalMomentumVector.set(particleFinalMomentum*sint3*cosp3,
			      particleFinalMomentum*sint3*sinp3,
			      particleFinalMomentum*cost3);

    G4LorentzVector muFourMomentum(particleMomentum, particleMomentumVector);
    G4LorentzVector muFinalFourMomentum(particleFinalEnergy, muFinalMomentumVector);
    G4LorentzVector auxVec1 = muFourMomentum - muFinalFourMomentum;
    G4double A6 = auxVec1.mag2() -
      2.*pEnergy*(kinEnergy - particleFinalEnergy) + 2.*particleMomentumVector[2]*pMomentum -
      2.*particleFinalMomentum*pMomentum*cost3;
    G4double B6 = -2.*particleFinalMomentum*pMomentum*(sint3*cosp3*cosp6 + sint3*sinp3*sinp6);
    G4double absA6 = std::abs(A6);
    G4double absB6 = std::abs(B6);
    G4double mint6 = 0.;
    G4double maxt6 = CLHEP::pi;
    G4double t6interval = G4Log((absA6 + absB6*maxt6)/(absA6 + absB6*mint6))/absB6;
    G4double argexp = absB6*t6interval*randNumbs[6] + G4Log(absA6 + absB6*mint6);
    G4double t6 = -absA6/absB6 + G4Exp(argexp)/absB6;
    G4double sint6 = std::sin(t6);
    G4double cost6 = std::cos(t6);

    pDirection.set(sint6*cosp6, sint6*sinp6, cost6);
    G4ThreeVector pMomentumVector = pMomentum*pDirection;

    G4ThreeVector auxVec2 = particleMomentumVector - muFinalMomentumVector - pMomentumVector;
    G4double p1mp3mp62 = auxVec2.dot(auxVec2);
    G4double Bp = particleFinalMomentum*(sint3*cosp3*cosp5 + sint3*sinp3*sinp5) +
      pMomentum*(sint6*cosp6*cosp5 + sint6*sinp6*sinp5);
    G4double Cp = -particleMomentum + particleFinalMomentum*cost3 + pMomentum*cost6;
    G4double A5 = p1mp3mp62 + eMomentum*eMomentum;
    G4double B5 = 2.*eMomentum*Bp;
    G4double C5 = 2.*eMomentum*Cp;
    G4double mint5 = 0.;
    G4double maxt5 = CLHEP::pi;
    G4double absA5C5 = std::abs(A5 + C5);
    G4double absB5 = std::abs(B5);
    G4double t5interval = (1./(absA5C5 + absB5*mint5) - 1./(absA5C5 + absB5*maxt5))/absB5;
    G4double t5 = (-absA5C5 + 1./(1./(absA5C5 + absB5*mint5) - absB5*t5interval*randNumbs[8]))/absB5;
    G4double sint5 = std::sin(t5);
    G4double cost5 = std::cos(t5);

    eDirection.set(sint5*cosp5, sint5*sinp5, cost5);
  }

  fAngularGenerator->Sample5DPairDirections(aDynamicParticle, eDirection, pDirection,
					    gEnergy, Q2, gMomentum,
					    particleFinalMomentum,
					    particleFinalEnergy,
					    randNumbs, W);
  
  // create G4DynamicParticle object for e+e-
  auto aParticle1 = new G4DynamicParticle(theElectron, eDirection, eEnergy);
  auto aParticle2 = new G4DynamicParticle(thePositron, pDirection, pEnergy);

  // Fill output vector
  vdp->push_back(aParticle1);
  vdp->push_back(aParticle2); 

  // if energy transfer is higher than threshold (very high by default)
  // then stop tracking the primary particle and create a new secondary
  if (pairEnergy > SecondaryThreshold()) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    auto newdp = new G4DynamicParticle(particle, muF);
    vdp->push_back(newdp);
  } else { // continue tracking the primary e-/e+ otherwise
    fParticleChange->SetProposedMomentumDirection(muF.vect().unit());
    G4double ekin = std::max(muF.e() - particleMass, 0.0);
    fParticleChange->SetProposedKineticEnergy(ekin);
  }
  //G4cout << "-- G4RiGeMuPairProductionModel::SampleSecondaries done" << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4RiGeMuPairProductionModel::FindScaledEnergy(G4int iz, G4double rand,
					      G4double logTkin,
					      G4double yymin, G4double yymax)
{
  G4double res = yymin;
  G4Physics2DVector* pv = fElementData->GetElement2DData(iz);
  if (nullptr != pv) { 
    G4double pmin = pv->Value(yymin, logTkin);
    G4double pmax = pv->Value(yymax, logTkin);
    G4double p0   = pv->Value(0.0, logTkin);
    if(p0 <= 0.0) { DataCorrupted(ZDATPAIR[iz], logTkin); }
    else { res = pv->FindLinearX((pmin + rand*(pmax - pmin))/p0, logTkin); }
  } else {
    DataCorrupted(ZDATPAIR[iz], logTkin); 
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RiGeMuPairProductionModel::DataCorrupted(G4int Z, G4double logTkin) const
{
  G4ExceptionDescription ed;
  ed << "G4ElementData is not properly initialized Z= " << Z
     << " Ekin(MeV)= " << G4Exp(logTkin)
     << " IsMasterThread= " << IsMaster() 
     << " Model " << GetName();
  G4Exception("G4RiGeMuPairProductionModel::()", "em0033", FatalException, ed, "");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RiGeMuPairProductionModel::StoreTables() const
{
  for (G4int iz=0; iz<NZDATPAIR; ++iz) {
    G4int Z = ZDATPAIR[iz];
    G4Physics2DVector* pv = fElementData->GetElement2DData(Z);
    if(nullptr == pv) { 
      DataCorrupted(Z, 1.0);
      return;
    }
    std::ostringstream ss;
    ss << "mupair/" << particle->GetParticleName() << Z << ".dat";
    std::ofstream outfile(ss.str());
    pv->Store(outfile);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4RiGeMuPairProductionModel::RetrieveTables()
{
  for (G4int iz=0; iz<NZDATPAIR; ++iz) {
    G4double Z = ZDATPAIR[iz];
    G4Physics2DVector* pv = new G4Physics2DVector(nbiny+1,nbine+1);
    std::ostringstream ss;
    ss << G4EmParameters::Instance()->GetDirLEDATA() << "/mupair/"
       << particle->GetParticleName() << Z << ".dat";
    std::ifstream infile(ss.str(), std::ios::in);
    if(!pv->Retrieve(infile)) { 
      delete pv;
      return false; 
    }
    fElementData->InitialiseForElement(iz, pv);
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

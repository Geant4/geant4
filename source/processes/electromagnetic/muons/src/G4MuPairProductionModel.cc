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
// File name:     G4MuPairProductionModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 24.06.2002
//
// Modifications:
//
// 04-12-02 Change G4DynamicParticle constructor in PostStep (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 24-01-03 Fix for compounds (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add model (V.Ivanchenko)
// 06-06-03 Fix in cross section calculation for high energy (V.Ivanchenko)
// 20-10-03 2*xi in ComputeDDMicroscopicCrossSection   (R.Kokoulin)
//          8 integration points in ComputeDMicroscopicCrossSection
// 12-01-04 Take min cut of e- and e+ not its sum (V.Ivanchenko)
// 10-02-04 Update parameterisation using R.Kokoulin model (V.Ivanchenko)
// 28-04-04 For complex materials repeat calculation of max energy for each
//          material (V.Ivanchenko)
// 01-11-04 Fix bug inside ComputeDMicroscopicCrossSection (R.Kokoulin)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 03-08-05 Add SetParticle method (V.Ivantchenko)
// 23-10-05 Add protection in sampling of e+e- pair energy needed for 
//          low cuts (V.Ivantchenko)
// 13-02-06 Add ComputeCrossSectionPerAtom (mma)
// 24-04-07 Add protection in SelectRandomAtom method (V.Ivantchenko)
// 12-05-06 Updated sampling (use cut) in SelectRandomAtom (A.Bogdanov) 
// 11-10-07 Add ignoreCut flag (V.Ivanchenko) 

//
// Class Description:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuPairProductionModel.hh"
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
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ModifiedMephi.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// static members
//
static const G4double ak1 = 6.9;
static const G4double ak2 = 1.0;
static const G4int    nzdat = 5;
static const G4int    zdat[5] = {1, 4, 13, 29, 92};

static const G4double xgi[] =
{ 0.0198550717512320, 0.1016667612931865, 0.2372337950418355, 0.4082826787521750,
  0.5917173212478250, 0.7627662049581645, 0.8983332387068135, 0.9801449282487680 };

static const G4double wgi[] =
{ 0.0506142681451880, 0.1111905172266870, 0.1568533229389435, 0.1813418916891810,
  0.1813418916891810, 0.1568533229389435, 0.1111905172266870, 0.0506142681451880 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuPairProductionModel::G4MuPairProductionModel(const G4ParticleDefinition* p,
                                                 const G4String& nam)
  : G4VEmModel(nam),
    factorForCross(CLHEP::fine_structure_const*CLHEP::fine_structure_const*
		   CLHEP::classic_electr_radius*CLHEP::classic_electr_radius*
		   4./(3.*CLHEP::pi)),
    sqrte(sqrt(G4Exp(1.))),
    minPairEnergy(4.*CLHEP::electron_mass_c2),
    lowestKinEnergy(0.85*CLHEP::GeV)
{
  nist = G4NistManager::Instance();

  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();

  if(nullptr != p) { 
    SetParticle(p); 
    lowestKinEnergy = std::max(lowestKinEnergy, p->GetPDGMass()*8.0);  
  }
  emin = lowestKinEnergy;
  emax = emin*10000.;
  SetAngularDistribution(new G4ModifiedMephi());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuPairProductionModel::MinPrimaryEnergy(const G4Material*,
                                                   const G4ParticleDefinition*,
                                                   G4double cut)
{
  return std::max(lowestKinEnergy, cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuPairProductionModel::Initialise(const G4ParticleDefinition* p,
                                         const G4DataVector& cuts)
{ 
  SetParticle(p); 

  if(nullptr == fParticleChange) { 
    fParticleChange = GetParticleChangeForLoss();
  }

  // for low-energy application this process should not work
  if(lowestKinEnergy >= HighEnergyLimit()) { return; }

  // define scale of internal table for each thread only once
  if(0 == nbine) {
    emin = std::max(lowestKinEnergy, LowEnergyLimit());
    emax = std::max(HighEnergyLimit(), emin*2);
    nbine = size_t(nYBinPerDecade*std::log10(emax/emin));
    if(nbine < 3) { nbine = 3; }

    ymin = G4Log(minPairEnergy/emin);
    dy   = -ymin/G4double(nbiny);
  }

  if(IsMaster() && p == particle) { 
    if(nullptr == fElementData) { 
      fElementData = new G4ElementData();
      G4bool dataFile = G4EmParameters::Instance()->RetrieveMuDataFromFile();
      if(dataFile)  { dataFile = RetrieveTables(); }
      if(!dataFile) { MakeSamplingTables(); }
      if(fTableToFile) { StoreTables(); }
    }    
    InitialiseElementSelectors(p, cuts); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuPairProductionModel::InitialiseLocal(const G4ParticleDefinition* p,
                                              G4VEmModel* masterModel)
{
  if(p == particle && lowestKinEnergy < HighEnergyLimit()) {
    SetElementSelectors(masterModel->GetElementSelectors());
    fElementData = masterModel->GetElementData();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuPairProductionModel::ComputeDEDXPerVolume(
                                              const G4Material* material,
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
  for (size_t i=0; i<material->GetNumberOfElements(); ++i) {
     G4double Z = (*theElementVector)[i]->GetZ();
     G4double tmax = MaxSecondaryEnergyForElement(kineticEnergy, Z);
     G4double loss = ComputMuPairLoss(Z, kineticEnergy, cutEnergy, tmax);
     dedx += loss*theAtomicNumDensityVector[i];
  }
  dedx = std::max(dedx, 0.0);
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuPairProductionModel::ComputMuPairLoss(G4double Z, 
                                                   G4double tkin,
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

  G4int kkk = G4lrint((bbb-aaa)/ak1+ak2);
  if(kkk > 8) { kkk = 8; }
  else if (kkk < 1) { kkk = 1; }
  G4double hhh = (bbb-aaa)/kkk;
  G4double x = aaa;

  for (G4int l=0 ; l<kkk; ++l) {
    for (G4int ll=0; ll<8; ++ll) {
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

G4double G4MuPairProductionModel::ComputeMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double cutEnergy)
{
  G4double cross = 0.;
  G4double tmax = MaxSecondaryEnergyForElement(tkin, Z);
  G4double cut  = std::max(cutEnergy, minPairEnergy);
  if (tmax <= cut) { return cross; }

  G4double aaa = G4Log(cut);
  G4double bbb = G4Log(tmax);
  G4int kkk = G4lrint((bbb-aaa)/ak1 + ak2);
  if(kkk > 8) { kkk = 8; }
  else if (kkk < 1) { kkk = 1; }

  G4double hhh = (bbb-aaa)/(kkk);
  G4double x = aaa;

  for(G4int l=0; l<kkk; ++l) {
    for(G4int i=0; i<8; ++i) {
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

G4double G4MuPairProductionModel::ComputeDMicroscopicCrossSection(
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
  G4double rt = sqrt(1.0 - alf);
  G4double delta = 6.0 * particleMass * particleMass * a0;
  G4double tmnexp = alf/(1.0 + rt) + delta*rt;

  if(tmnexp >= 1.0) { return 0.0; }

  G4double tmn = G4Log(tmnexp);

  G4double massratio      = particleMass/electron_mass_c2;
  G4double massratio2     = massratio*massratio;
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
  G4double rho[8];
  G4double rho2[8];
  G4double xi[8];
  G4double xi1[8];
  G4double xii[8];

  for (G4int i = 0; i < 8; ++i)
  {
    rho[i] = G4Exp(tmn*xgi[i]) - 1.0; // rho = -asymmetry
    rho2[i] = rho[i] * rho[i];
    xi[i] = xi0*(1.0-rho2[i]);
    xi1[i] = 1.0 + xi[i];
    xii[i] = 1.0 / xi[i];
  }

  G4double ye1[8];
  G4double ym1[8];

  G4double b40 = 4.0 * beta;
  G4double b62 = 6.0 * beta + 2.0;

  for (G4int i = 0; i < 8; ++i)
  {
    G4double yeu = (b40 + 5.0) + (b40 - 1.0) * rho2[i];
    G4double yed = b62*G4Log(3.0 + xii[i]) + (2.0 * beta - 1.0)*rho2[i] - b40;

    G4double ymu = b62 * (1.0 + rho2[i]) + 6.0;
    G4double ymd = (b40 + 3.0)*(1.0 + rho2[i])*G4Log(3.0 + xi[i])
      + 2.0 - 3.0 * rho2[i];

    ye1[i] = 1.0 + yeu / yed;
    ym1[i] = 1.0 + ymu / ymd;
  }

  G4double be[8];
  G4double bm[8];

  for(G4int i = 0; i < 8; ++i) {
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

  for (G4int i = 0; i < 8; ++i) {
    G4double screen = screen0*xi1[i]/(1.0 - rho2[i]);
    G4double ale = G4Log(bbb/z13*sqrt(xi1[i]*ye1[i])/(1. + screen*ye1[i]));
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

G4double G4MuPairProductionModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition*,
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

void G4MuPairProductionModel::MakeSamplingTables()
{
  G4double factore = G4Exp(G4Log(emax/emin)/G4double(nbine));

  for (G4int iz=0; iz<nzdat; ++iz) {

    G4double Z = zdat[iz];
    G4Physics2DVector* pv = new G4Physics2DVector(nbiny+1,nbine+1);
    G4double kinEnergy = emin;

    for (size_t it=0; it<=nbine; ++it) {

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
      size_t imax   = (size_t)fac;
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

      for (size_t i=0; i<nbiny; ++i) {

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
    fElementData->InitialiseForElement(zdat[iz], pv);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuPairProductionModel::SampleSecondaries(
                              std::vector<G4DynamicParticle*>* vdp, 
                              const G4MaterialCutsCouple* couple,
                              const G4DynamicParticle* aDynamicParticle,
                              G4double tmin,
                              G4double tmax)
{
  G4double kinEnergy = aDynamicParticle->GetKineticEnergy();
  //G4cout << "------- G4MuPairProductionModel::SampleSecondaries E(MeV)= " 
  //         << kinEnergy << "  " 
  //         << aDynamicParticle->GetDefinition()->GetParticleName() << G4endl;
  G4double totalEnergy   = kinEnergy + particleMass;
  G4double totalMomentum = 
    sqrt(kinEnergy*(kinEnergy + 2.0*particleMass));

  G4ThreeVector partDirection = aDynamicParticle->GetMomentumDirection();

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(couple,particle,kinEnergy);

  // define interval of energy transfer
  G4double maxPairEnergy = MaxSecondaryEnergyForElement(kinEnergy, 
                                                        anElement->GetZ());
  G4double maxEnergy = std::min(tmax, maxPairEnergy);
  G4double minEnergy = std::max(tmin, minPairEnergy);

  if(minEnergy >= maxEnergy) { return; }
  //G4cout << "emin= " << minEnergy << " emax= " << maxEnergy 
  // << " minPair= " << minPairEnergy << " maxpair= " << maxPairEnergy 
  //    << " ymin= " << ymin << " dy= " << dy << G4endl;

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
  for(G4int iz=0; iz<nzdat; ++iz) { 
    if(currentZ == zdat[iz]) {
      iz1 = iz2 = currentZ; 
      break;
    } else if(currentZ < zdat[iz]) {
      iz2 = zdat[iz];
      if(iz > 0) { iz1 = zdat[iz-1]; }
      else { iz1 = iz2; }
      break;
    } 
  }
  if(0 == iz1) { iz1 = iz2 = zdat[nzdat-1]; }

  G4double pairEnergy = 0.0;
  G4int count = 0;
  //G4cout << "start loop Z1= " << iz1 << " Z2= " << iz2 << G4endl;
  do {
    ++count;
    // sampling using only one random number
    G4double rand = G4UniformRand();
  
    G4double x = FindScaledEnergy(iz1, rand, logTkin, yymin, yymax);
    if(iz1 != iz2) {
      G4double x2 = FindScaledEnergy(iz2, rand, logTkin, yymin, yymax);
      G4double lz1= nist->GetLOGZ(iz1);
      G4double lz2= nist->GetLOGZ(iz2);
      //G4cout << count << ".  x= " << x << "  x2= " << x2 
      //             << " Z1= " << iz1 << " Z2= " << iz2 << G4endl;
      x += (x2 - x)*(lnZ - lz1)/(lz2 - lz1);
    }
    //G4cout << "x= " << x << "  coeff= " << coeff << G4endl;
    pairEnergy = kinEnergy*G4Exp(x*coeff);
    
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while((pairEnergy < minEnergy || pairEnergy > maxEnergy) && 10 > count);

  //G4cout << "## pairEnergy(GeV)= " << pairEnergy/GeV 
  //         << " Etot(GeV)= " << totalEnergy/GeV << G4endl; 

  // sample r=(E+-E-)/pairEnergy  ( uniformly .....)
  G4double rmax =
    (1.-6.*particleMass*particleMass/(totalEnergy*(totalEnergy-pairEnergy)))
                                       *sqrt(1.-minPairEnergy/pairEnergy);
  G4double r = rmax * (-1.+2.*G4UniformRand()) ;

  // compute energies from pairEnergy,r
  G4double eEnergy = (1.-r)*pairEnergy*0.5;
  G4double pEnergy = pairEnergy - eEnergy;

  // Sample angles 
  G4ThreeVector eDirection, pDirection;
  //
  GetAngularDistribution()->SamplePairDirections(aDynamicParticle, 
                                                 eEnergy, pEnergy,
                                                 eDirection, pDirection);
  // create G4DynamicParticle object for e+e-
  eEnergy = std::max(eEnergy - CLHEP::electron_mass_c2, 0.0);
  pEnergy = std::max(pEnergy - CLHEP::electron_mass_c2, 0.0);
  G4DynamicParticle* aParticle1 =
    new G4DynamicParticle(theElectron,eDirection,eEnergy);
  G4DynamicParticle* aParticle2 = 
    new G4DynamicParticle(thePositron,pDirection,pEnergy);
  // Fill output vector
  vdp->push_back(aParticle1);
  vdp->push_back(aParticle2);

  // primary change
  kinEnergy -= pairEnergy;
  partDirection *= totalMomentum;
  partDirection -= (aParticle1->GetMomentum() + aParticle2->GetMomentum());
  partDirection = partDirection.unit();

  // if energy transfer is higher than threshold (very high by default)
  // then stop tracking the primary particle and create a new secondary
  if (pairEnergy > SecondaryThreshold()) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    G4DynamicParticle* newdp = 
      new G4DynamicParticle(particle, partDirection, kinEnergy);
    vdp->push_back(newdp);
  } else { // continue tracking the primary e-/e+ otherwise
    fParticleChange->SetProposedMomentumDirection(partDirection);
    fParticleChange->SetProposedKineticEnergy(kinEnergy);
  }
  //G4cout << "-- G4MuPairProductionModel::SampleSecondaries done" << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4MuPairProductionModel::FindScaledEnergy(G4int Z, G4double rand,
					  G4double logTkin,
					  G4double yymin, G4double yymax)
{
  G4double res = yymin;
  G4Physics2DVector* pv = fElementData->GetElement2DData(Z);
  if(nullptr != pv) { 
    G4double pmin = pv->Value(yymin, logTkin);
    G4double pmax = pv->Value(yymax, logTkin);
    G4double p0   = pv->Value(0.0, logTkin);
    if(p0 <= 0.0) { DataCorrupted(Z, logTkin); }
    else { res = pv->FindLinearX((pmin + rand*(pmax - pmin))/p0, logTkin); }
  } else {
    DataCorrupted(Z, logTkin); 
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuPairProductionModel::DataCorrupted(G4int Z, G4double logTkin) const
{
  G4ExceptionDescription ed;
  ed << "G4ElementData is not properly initialized Z= " << Z
     << " Ekin(MeV)= " << G4Exp(logTkin)
     << " IsMasterThread= " << IsMaster() 
     << " Model " << GetName();
  G4Exception("G4MuPairProductionModel::()", "em0033", FatalException, ed, "");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuPairProductionModel::StoreTables() const
{
  for (G4int iz=0; iz<nzdat; ++iz) {
    G4int Z = zdat[iz];
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

G4bool G4MuPairProductionModel::RetrieveTables()
{
  const char* path = G4FindDataDir("G4LEDATA");
  G4String dir("");
  if (path) { 
    std::ostringstream ost;
    ost << path << "/mupair/";
    dir = ost.str(); 
  } else {
    dir = "./mupair/";
  }

  for (G4int iz=0; iz<nzdat; ++iz) {
    G4double Z = zdat[iz];
    G4Physics2DVector* pv = new G4Physics2DVector(nbiny+1,nbine+1);
    std::ostringstream ss;
    ss << dir << particle->GetParticleName() << Z << ".dat";
    std::ifstream infile(ss.str(), std::ios::in);
    if(!pv->Retrieve(infile)) { 
      delete pv;
      return false; 
    }
    fElementData->InitialiseForElement(Z, pv);
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

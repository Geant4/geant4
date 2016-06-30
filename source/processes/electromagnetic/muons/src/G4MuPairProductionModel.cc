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
// $Id: G4MuPairProductionModel.cc 97392 2016-06-02 10:10:32Z gcosmo $
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
#include "G4ParticleChangeForGamma.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// static members
//
static const G4double ak1 = 6.9;
static const G4double ak2 = 1.0;
static const G4int    zdat[5] = {1, 4, 13, 29, 92};

const G4double G4MuPairProductionModel::xgi[] = 
  { 0.0199, 0.1017, 0.2372, 0.4083, 0.5917, 0.7628, 0.8983, 0.9801 };
const G4double G4MuPairProductionModel::wgi[8] = 
  { 0.0506, 0.1112, 0.1569, 0.1813, 0.1813, 0.1569, 0.1112, 0.0506 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuPairProductionModel::G4MuPairProductionModel(const G4ParticleDefinition* p,
                                                 const G4String& nam)
  : G4VEmModel(nam),
    particle(nullptr),
    factorForCross(4.*fine_structure_const*fine_structure_const
                   *classic_electr_radius*classic_electr_radius/(3.*pi)),
    sqrte(sqrt(G4Exp(1.))),
    currentZ(0),
    fParticleChange(nullptr),
    minPairEnergy(4.*electron_mass_c2),
    lowestKinEnergy(1.0*GeV),
    nzdat(5),
    nYBinPerDecade(4),
    nbiny(1000),
    nbine(0),
    ymin(-5.),
    dy(0.005)
{
  nist = G4NistManager::Instance();

  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();

  particleMass = lnZ = z13 = z23 = 0.;

  // setup lowest limit dependent on particle mass
  if(p) { 
    SetParticle(p); 
    lowestKinEnergy = std::max(lowestKinEnergy,p->GetPDGMass()*8.0); 
  }
  emin = lowestKinEnergy;
  emax = 10.*TeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuPairProductionModel::~G4MuPairProductionModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuPairProductionModel::MinPrimaryEnergy(const G4Material*,
                                                   const G4ParticleDefinition*,
                                                   G4double cut)
{
  return std::max(lowestKinEnergy,cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuPairProductionModel::Initialise(const G4ParticleDefinition* p,
                                         const G4DataVector& cuts)
{ 
  SetParticle(p); 
  if(!fParticleChange) { fParticleChange = GetParticleChangeForLoss(); }

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
    
    if(!fElementData) { 
      fElementData = new G4ElementData();
      MakeSamplingTables(); 
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

  G4double cut = std::min(cutEnergy,tmax);
  if(cut <= minPairEnergy) { return loss; }

  // calculate the rectricted loss
  // numerical integration in log(PairEnergy)
  G4double aaa = G4Log(minPairEnergy);
  G4double bbb = G4Log(cut);

  G4int    kkk = (G4int)((bbb-aaa)/ak1+ak2);
  if (kkk > 8)      { kkk = 8; }
  else if (kkk < 1) { kkk = 1; }

  G4double hhh = (bbb-aaa)/(G4double)kkk;
  G4double x = aaa;

  for (G4int l=0 ; l<kkk; l++)
  {

    for (G4int ll=0; ll<8; ll++)
    {
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

  //  G4double ak1=6.9 ;
  // G4double ak2=1.0 ;
  G4double aaa = G4Log(cut);
  G4double bbb = G4Log(tmax);
  G4int kkk = (G4int)((bbb-aaa)/ak1 + ak2);
  if(kkk > 8) { kkk = 8; }
  else if (kkk < 1) { kkk = 1; }

  G4double hhh = (bbb-aaa)/G4double(kkk);
  G4double x = aaa;

  for(G4int l=0; l<kkk; ++l)
  {
    for(G4int i=0; i<8; ++i)
    {
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

  G4double totalEnergy  = tkin + particleMass;
  G4double residEnergy  = totalEnergy - pairEnergy;
  G4double massratio    = particleMass/electron_mass_c2;
  G4double massratio2   = massratio*massratio;
  G4double cross = 0.;

  G4double c3 = 0.75*sqrte*particleMass;
  if (residEnergy <= c3*z13) { return cross; }

  G4double c7 = 4.*CLHEP::electron_mass_c2;
  G4double c8 = 6.*particleMass*particleMass;
  G4double alf = c7/pairEnergy;
  G4double a3 = 1. - alf;
  if (a3 <= 0.) { return cross; }

  // zeta calculation
  G4double bbb,g1,g2;
  if( Z < 1.5 ) { bbb = bbbh ; g1 = g1h ; g2 = g2h ; }
  else          { bbb = bbbtf; g1 = g1tf; g2 = g2tf; }

  G4double zeta = 0;
  G4double zeta1 = 
    0.073*G4Log(totalEnergy/(particleMass+g1*z23*totalEnergy))-0.26;
  if ( zeta1 > 0.)
  {
    G4double zeta2 = 
      0.058*G4Log(totalEnergy/(particleMass+g2*z13*totalEnergy))-0.14;
    zeta  = zeta1/zeta2 ;
  }

  G4double z2 = Z*(Z+zeta);
  G4double screen0 = 2.*electron_mass_c2*sqrte*bbb/(z13*pairEnergy);
  G4double a0 = totalEnergy*residEnergy;
  G4double a1 = pairEnergy*pairEnergy/a0;
  G4double bet = 0.5*a1;
  G4double xi0 = 0.25*massratio2*a1;
  G4double del = c8/a0;

  G4double rta3 = sqrt(a3);
  G4double tmnexp = alf/(1. + rta3) + del*rta3;
  if(tmnexp >= 1.0) { return cross; }

  G4double tmn = G4Log(tmnexp);
  G4double sum = 0.;

  // Gaussian integration in ln(1-ro) ( with 8 points)
  for (G4int i=0; i<8; ++i)
  {
    G4double a4 = G4Exp(tmn*xgi[i]);     // a4 = (1.-asymmetry)
    G4double a5 = a4*(2.-a4) ;
    G4double a6 = 1.-a5 ;
    G4double a7 = 1.+a6 ;
    G4double a9 = 3.+a6 ;
    G4double xi = xi0*a5 ;
    G4double xii = 1./xi ;
    G4double xi1 = 1.+xi ;
    G4double screen = screen0*xi1/a5 ;
    G4double yeu = 5.-a6+4.*bet*a7 ;
    G4double yed = 2.*(1.+3.*bet)*G4Log(3.+xii)-a6-a1*(2.-a6) ;
    G4double ye1 = 1.+yeu/yed ;
    G4double ale = G4Log(bbb/z13*sqrt(xi1*ye1)/(1.+screen*ye1)) ;
    G4double cre = 0.5*G4Log(1.+2.25*z23*xi1*ye1/massratio2) ;
    G4double be;

    if (xi <= 1.e3) { 
      be = ((2.+a6)*(1.+bet)+xi*a9)*G4Log(1.+xii)+(a5-bet)/xi1-a9;
    } else {           
      be = (3.-a6+a1*a7)/(2.*xi);
    }
    G4double fe = (ale-cre)*be;
    if ( fe < 0.) fe = 0. ;

    G4double ymu = 4.+a6 +3.*bet*a7 ;
    G4double ymd = a7*(1.5+a1)*G4Log(3.+xi)+1.-1.5*a6 ;
    G4double ym1 = 1.+ymu/ymd ;
    G4double alm_crm = G4Log(bbb*massratio/(1.5*z23*(1.+screen*ym1)));
    G4double a10,bm;
    if ( xi >= 1.e-3)
    {
      a10 = (1.+a1)*a5 ;
      bm  = (a7*(1.+1.5*bet)-a10*xii)*G4Log(xi1)+xi*(a5-bet)/xi1+a10;
    } else {
      bm = (5.-a6+bet*a9)*(xi/2.);
    }

    G4double fm = alm_crm*bm;
    if ( fm < 0.) { fm = 0.; }

    sum += wgi[i]*a4*(fe+fm/massratio2);
  }

  cross = -tmn*sum*factorForCross*z2*residEnergy/(totalEnergy*pairEnergy);
  cross = std::max(cross, 0.0); 
  return cross;
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

      pv->PutY(it, G4Log(kinEnergy/MeV));
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
  G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();
  //G4cout << "------- G4MuPairProductionModel::SampleSecondaries E(MeV)= " 
  //         << kineticEnergy << "  " 
  //         << aDynamicParticle->GetDefinition()->GetParticleName() << G4endl;
  G4double totalEnergy   = kineticEnergy + particleMass;
  G4double totalMomentum = 
    sqrt(kineticEnergy*(kineticEnergy + 2.0*particleMass));

  G4ThreeVector partDirection = aDynamicParticle->GetMomentumDirection();

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(couple,particle,kineticEnergy);

  // define interval of energy transfer
  G4double maxPairEnergy = MaxSecondaryEnergyForElement(kineticEnergy, 
                                                        anElement->GetZ());
  G4double maxEnergy     = std::min(tmax, maxPairEnergy);
  G4double minEnergy     = std::max(tmin, minPairEnergy);

  if(minEnergy >= maxEnergy) { return; }
  //G4cout << "emin= " << minEnergy << " emax= " << maxEnergy 
  // << " minPair= " << minPairEnergy << " maxpair= " << maxPairEnergy 
  //    << " ymin= " << ymin << " dy= " << dy << G4endl;

  G4double coeff = G4Log(minPairEnergy/kineticEnergy)/ymin;

  // compute limits 
  G4double yymin = G4Log(minEnergy/kineticEnergy)/coeff;
  G4double yymax = G4Log(maxEnergy/kineticEnergy)/coeff;
 
  //G4cout << "yymin= " << yymin << "  yymax= " << yymax << G4endl;

  // units should not be used, bacause table was built without
  G4double logTkin = G4Log(kineticEnergy/MeV);

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

  G4double PairEnergy = 0.0;
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
    PairEnergy = kineticEnergy*G4Exp(x*coeff);
    
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while((PairEnergy < minEnergy || PairEnergy > maxEnergy) && 10 > count);

  //G4cout << "## PairEnergy(GeV)= " << PairEnergy/GeV 
  //         << " Etot(GeV)= " << totalEnergy/GeV << G4endl; 

  // sample r=(E+-E-)/PairEnergy  ( uniformly .....)
  G4double rmax =
    (1.-6.*particleMass*particleMass/(totalEnergy*(totalEnergy-PairEnergy)))
                                       *sqrt(1.-minPairEnergy/PairEnergy);
  G4double r = rmax * (-1.+2.*G4UniformRand()) ;

  // compute energies from PairEnergy,r
  G4double ElectronEnergy = (1.-r)*PairEnergy*0.5;
  G4double PositronEnergy = PairEnergy - ElectronEnergy;

  // The angle of the emitted virtual photon is sampled
  // according to the muon bremsstrahlung model
 
  G4double gam  = totalEnergy/particleMass;
  G4double gmax = gam*std::min(1.0, totalEnergy/PairEnergy - 1.0);
  G4double gmax2= gmax*gmax;
  G4double x = G4UniformRand()*gmax2/(1.0 + gmax2);

  G4double theta = sqrt(x/(1.0 - x))/gam;
  G4double sint  = sin(theta);
  G4double phi   = twopi * G4UniformRand() ;
  G4double dirx  = sint*cos(phi), diry = sint*sin(phi), dirz = cos(theta) ;

  G4ThreeVector gDirection(dirx, diry, dirz);
  gDirection.rotateUz(partDirection);

  // the angles of e- and e+ assumed to be the same as virtual gamma

  // create G4DynamicParticle object for the particle1
  G4DynamicParticle* aParticle1 = 
    new G4DynamicParticle(theElectron, gDirection, 
                          ElectronEnergy - electron_mass_c2);

  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2 = 
    new G4DynamicParticle(thePositron, gDirection,
                          PositronEnergy - electron_mass_c2);

  // primary change
  kineticEnergy -= (ElectronEnergy + PositronEnergy);
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);

  partDirection *= totalMomentum;
  partDirection -= (aParticle1->GetMomentum() + aParticle2->GetMomentum());
  partDirection = partDirection.unit();
  fParticleChange->SetProposedMomentumDirection(partDirection);

  // add secondary
  vdp->push_back(aParticle1);
  vdp->push_back(aParticle2);
  //G4cout << "-- G4MuPairProductionModel::SampleSecondaries done" << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuPairProductionModel::DataCorrupted(G4int Z, G4double logTkin)
{
  G4ExceptionDescription ed;
  ed << "G4ElementData is not properly initialized Z= " << Z
     << " Ekin(MeV)= " << G4Exp(logTkin)
     << " IsMasterThread= " << IsMaster() 
     << " Model " << GetName();
  G4Exception("G4MuPairProductionModel::()","em0033",FatalException,
              ed,"");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



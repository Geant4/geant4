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
// $Id: G4BraggIonModel.cc 83008 2014-07-24 14:49:52Z gcosmo $
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
#include "G4DeltaAngle.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4ASTARStopping* G4BraggIonModel::fASTAR = 0;

G4BraggIonModel::G4BraggIonModel(const G4ParticleDefinition* p,
                                 const G4String& nam)
  : G4VEmModel(nam),
    corr(0),
    particle(0),
    fParticleChange(0),
    currentMaterial(0),
    iMolecula(-1),
    iASTAR(-1),
    isIon(false),
    isInitialised(false)
{
  SetHighEnergyLimit(2.0*MeV);

  HeMass           = 3.727417*GeV;
  rateMassHe2p     = HeMass/proton_mass_c2;
  lowestKinEnergy  = 1.0*keV/rateMassHe2p;
  massFactor       = 1000.*amu_c2/HeMass;
  theZieglerFactor = eV*cm2*1.0e-15;
  theElectron      = G4Electron::Electron();
  corrFactor       = 1.0;
  if(p) { SetParticle(p); }
  else  { SetParticle(theElectron); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BraggIonModel::~G4BraggIonModel()
{
  if(IsMaster()) { delete fASTAR; fASTAR = 0; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggIonModel::Initialise(const G4ParticleDefinition* p,
                                 const G4DataVector&)
{
  if(p != particle) { SetParticle(p); }

  corrFactor = chargeSquare;

  // always false before the run
  SetDeexcitationFlag(false);

  if(!isInitialised) {
    isInitialised = true;

    if(UseAngularGeneratorFlag() && !GetAngularDistribution()) {
      SetAngularDistribution(new G4DeltaAngle());
    }
    G4String pname = particle->GetParticleName();
    if(particle->GetParticleType() == "nucleus" &&
       pname != "deuteron" && pname != "triton" &&
       pname != "alpha+"   && pname != "helium" &&
       pname != "hydrogen") { isIon = true; }

    corr = G4LossTableManager::Instance()->EmCorrections();

    fParticleChange = GetParticleChangeForLoss();
    if(!fASTAR) { fASTAR = new G4ASTARStopping(); }
  }
  if(IsMaster() && particle->GetPDGMass() < GeV) { fASTAR->Initialise(); }
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
  //G4cout<<"G4BraggIonModel::GetChargeSquareRatio e= "<<kineticEnergy<<G4endl;
  // this method is called only for ions
  G4double q2 = corr->EffectiveChargeSquareRatio(p,mat,kineticEnergy);
  corrFactor  = q2*corr->EffectiveChargeCorrection(p,mat,kineticEnergy); 
  return corrFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::GetParticleCharge(const G4ParticleDefinition* p,
					    const G4Material* mat,
					    G4double kineticEnergy)
{
  //G4cout<<"G4BraggIonModel::GetParticleCharge e= "<<kineticEnergy <<
  //  " q= " <<  corr->GetParticleCharge(p,mat,kineticEnergy) <<G4endl;
  // this method is called only for ions
  return corr->GetParticleCharge(p,mat,kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::ComputeCrossSectionPerElectron(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxKinEnergy)
{
  G4double cross     = 0.0;
  G4double tmax      = MaxSecondaryEnergy(p, kineticEnergy);
  G4double maxEnergy = std::min(tmax,maxKinEnergy);
  if(cutEnergy < tmax) {

    G4double energy  = kineticEnergy + mass;
    G4double energy2 = energy*energy;
    G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;
    cross = (maxEnergy - cutEnergy)/(cutEnergy*maxEnergy) 
      - beta2*G4Log(maxEnergy/cutEnergy)/tmax;

    if( 0.5 == spin ) { cross += 0.5*(maxEnergy - cutEnergy)/energy2; }

    cross *= twopi_mc2_rcl2*chargeSquare/beta2;
  }
 //   G4cout << "BR: e= " << kineticEnergy << " tmin= " << cutEnergy 
 //          << " tmax= " << tmax << " cross= " << cross << G4endl;
 
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
						 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double cross = Z*ComputeCrossSectionPerElectron
                                         (p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::CrossSectionPerVolume(
					   const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double eDensity = material->GetElectronDensity();
  G4double cross = eDensity*ComputeCrossSectionPerElectron
                                         (p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::ComputeDEDXPerVolume(const G4Material* material,
					       const G4ParticleDefinition* p,
					       G4double kineticEnergy,
					       G4double cutEnergy)
{
  G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tmin  = min(cutEnergy, tmax);
  G4double tkin  = kineticEnergy/massRate;
  G4double dedx  = 0.0;

  if(tkin < lowestKinEnergy) {
    dedx = DEDX(material, lowestKinEnergy)*sqrt(tkin/lowestKinEnergy);
  } else {
    dedx = DEDX(material, tkin); 
  }

  if (cutEnergy < tmax) {

    G4double tau   = kineticEnergy/mass;
    G4double gam   = tau + 1.0;
    G4double bg2   = tau * (tau+2.0);
    G4double beta2 = bg2/(gam*gam);
    G4double x     = tmin/tmax;

    dedx += (G4Log(x) + (1.0 - x)*beta2) * twopi_mc2_rcl2
          * (material->GetElectronDensity())/beta2;
  }

  // now compute the total ionization loss

  if (dedx < 0.0) dedx = 0.0 ;

  dedx *= chargeSquare;

  //G4cout << " tkin(MeV) = " << tkin/MeV << " dedx(MeVxcm^2/g) = " 
  //       << dedx*gram/(MeV*cm2*material->GetDensity()) 
  //       << " q2 = " << chargeSquare <<  G4endl;

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggIonModel::CorrectionsAlongStep(const G4MaterialCutsCouple* couple,
					   const G4DynamicParticle* dp,
					   G4double& eloss,
					   G4double&,
					   G4double /*length*/)
{
  // this method is called only for ions
  const G4ParticleDefinition* p = dp->GetDefinition();
  const G4Material* mat = couple->GetMaterial();
  G4double preKinEnergy = dp->GetKineticEnergy();
  G4double e = preKinEnergy - eloss*0.5;
  if(e < 0.0) { e = preKinEnergy*0.5; }

  G4double q2 = corr->EffectiveChargeSquareRatio(p,mat,e);
  GetModelOfFluctuations()->SetParticleAndCharge(p, q2);
  G4double qfactor = q2*corr->EffectiveChargeCorrection(p,mat,e)/corrFactor; 
  eloss *= qfactor; 

  //G4cout << "G4BraggIonModel::CorrectionsAlongStep e= " <<  e 
  //	 << " qfactor= " << qfactor << " " << p->GetParticleName() <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BraggIonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
					const G4MaterialCutsCouple* couple,
					const G4DynamicParticle* dp,
					G4double xmin,
					G4double maxEnergy)
{
  G4double tmax = MaxSecondaryKinEnergy(dp);
  G4double xmax = std::min(tmax, maxEnergy);
  if(xmin >= xmax) { return; }

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double energy  = kineticEnergy + mass;
  G4double energy2 = energy*energy;
  G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;
  G4double grej    = 1.0;
  G4double deltaKinEnergy, f;

  // sampling follows ...
  do {
    G4double q = G4UniformRand();
    deltaKinEnergy = xmin*xmax/(xmin*(1.0 - q) + xmax*q);

    f = 1.0 - beta2*deltaKinEnergy/tmax;

    if(f > grej) {
        G4cout << "G4BraggIonModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for e= " << deltaKinEnergy
               << G4endl;
    }

  } while( grej*G4UniformRand() >= f );

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

    G4double phi = twopi * G4UniformRand() ;

    deltaDirection.set(sint*cos(phi),sint*sin(phi), cost) ;
    deltaDirection.rotateUz(dp->GetMomentumDirection());
  }  

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

G4bool G4BraggIonModel::HasMaterial(const G4Material*)
{
  return false;
  /*
  G4String chFormula = material->GetChemicalFormula();
  if("" == chFormula) { return false; }

  // ICRU Report N49, 1993. Ziegler model for He.
  
  static const size_t numberOfMolecula = 11;
  static const G4String molName[numberOfMolecula] = {
    "CaF_2",  "Cellulose_Nitrate",  "LiF", "Policarbonate",  
    "(C_2H_4)_N-Polyethylene",  "(C_2H_4)_N-Polymethly_Methacralate",
    "Polysterene", "SiO_2", "NaI", "H_2O",
    "Graphite" };

  // Search for the material in the table
  for (size_t i=0; i<numberOfMolecula; ++i) {
    if (chFormula == molName[i]) {
      iASTAR = fASTAR->GetIndex(matName[i]);  
      break;
    }
  }
  return (iASTAR >= 0);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::StoppingPower(const G4Material* material,
					G4double kineticEnergy) 
{
  G4double ionloss = 0.0 ;

  if (iMolecula >= 0) {
  
    // The data and the fit from: 
    // ICRU Report N49, 1993. Ziegler's model for alpha
    // He energy in internal units of parametrisation formula (MeV)

    G4double T = kineticEnergy*rateMassHe2p/MeV ;

    static const G4double a[11][5] = {
       {9.43672, 0.54398, 84.341, 1.3705, 57.422},
       {67.1503, 0.41409, 404.512, 148.97, 20.99},
       {5.11203, 0.453,  36.718,  50.6,  28.058}, 
       {61.793, 0.48445, 361.537, 57.889, 50.674},
       {7.83464, 0.49804, 160.452, 3.192, 0.71922},
       {19.729, 0.52153, 162.341, 58.35, 25.668}, 
       {26.4648, 0.50112, 188.913, 30.079, 16.509},
       {7.8655, 0.5205, 63.96, 51.32, 67.775},
       {8.8965, 0.5148, 339.36, 1.7205, 0.70423},
       {2.959, 0.53255, 34.247, 60.655, 15.153}, 
       {3.80133, 0.41590, 12.9966, 117.83, 242.28} };   

    static const G4double atomicWeight[11] = {
       101.96128, 44.0098, 16.0426, 28.0536, 42.0804,
       104.1512, 44.665, 60.0843, 18.0152, 18.0152, 12.0};       

    G4int i = iMolecula;

    // Free electron gas model
    if ( T < 0.001 ) {
      G4double slow  = a[i][0] ;
      G4double shigh = G4Log( 1.0 + a[i][3]*1000.0 + a[i][4]*0.001 )
	 * a[i][2]*1000.0 ;
      ionloss  = slow*shigh / (slow + shigh) ;
      ionloss *= sqrt(T*1000.0) ;

      // Main parametrisation
    } else {
      G4double slow  = a[i][0] * G4Exp(G4Log(T*1000.0)*a[i][1]) ;
      G4double shigh = G4Log( 1.0 + a[i][3]/T + a[i][4]*T ) * a[i][2]/T ;
      ionloss = slow*shigh / (slow + shigh) ;
       /*
	 G4cout << "## " << i << ". T= " << T << " slow= " << slow
	 << " a0= " << a[i][0] << " a1= " << a[i][1] 
	 << " shigh= " << shigh 
	 << " dedx= " << ionloss << " q^2= " <<  HeEffChargeSquare(z, T*MeV)
	 << G4endl;
       */
    }
    if ( ionloss < 0.0) ionloss = 0.0 ;

    // He effective charge
    G4double aa = atomicWeight[iMolecula];
    ionloss /= (HeEffChargeSquare(0.5*aa, T)*aa);

  // pure material (normally not the case for this function)
  } else if(1 == (material->GetNumberOfElements())) {
    G4double z = material->GetZ() ;
    ionloss = ElectronicStoppingPower( z, kineticEnergy ) ;  
  }
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::ElectronicStoppingPower(G4double z,
                                                  G4double kineticEnergy) const
{
  G4double ionloss ;
  G4int i = G4int(z)-1 ;  // index of atom
  if(i < 0)  i = 0 ;
  if(i > 91) i = 91 ;

  // The data and the fit from:
  // ICRU Report 49, 1993. Ziegler's type of parametrisations.
  // Proton kinetic energy for parametrisation (keV/amu)

   // He energy in internal units of parametrisation formula (MeV)
  G4double T = kineticEnergy*rateMassHe2p/MeV ;

  static const G4double a[92][5] = {
    {0.35485, 0.6456, 6.01525,  20.8933, 4.3515
   },{ 0.58,    0.59,   6.3,     130.0,   44.07
   },{ 1.42,    0.49,   12.25,    32.0,    9.161
   },{ 2.206,   0.51,   15.32,    0.25,    8.995 //Be Ziegler77
       // },{ 2.1895,  0.47183,7.2362,   134.30,  197.96 //Be from ICRU
   },{ 3.691,   0.4128, 18.48,    50.72,   9.0
   },{ 3.83523, 0.42993,12.6125,  227.41,  188.97
   },{ 1.9259,  0.5550, 27.15125, 26.0665, 6.2768
   },{ 2.81015, 0.4759, 50.0253,  10.556,  1.0382
   },{ 1.533,   0.531,  40.44,    18.41,   2.718
   },{ 2.303,   0.4861, 37.01,    37.96,   5.092
       // Z= 11-20
   },{ 9.894,   0.3081, 23.65,    0.384,   92.93
   },{ 4.3,     0.47,   34.3,     3.3,     12.74
   },{ 2.5,     0.625,  45.7,     0.1,     4.359
   },{ 2.1,     0.65,   49.34,    1.788,   4.133
   },{ 1.729,   0.6562, 53.41,    2.405,   3.845
   },{ 1.402,   0.6791, 58.98,    3.528,   3.211
   },{ 1.117,   0.7044, 69.69,    3.705,    2.156
   },{ 2.291,   0.6284, 73.88,    4.478,    2.066
   },{ 8.554,   0.3817, 83.61,    11.84,    1.875
   },{ 6.297,   0.4622, 65.39,    10.14,    5.036
       // Z= 21-30     
   },{ 5.307,   0.4918, 61.74,    12.4,    6.665
   },{ 4.71,    0.5087, 65.28,    8.806,    5.948
   },{ 6.151,   0.4524, 83.0,    18.31,    2.71
   },{ 6.57,    0.4322, 84.76,    15.53,    2.779
   },{ 5.738,   0.4492, 84.6,    14.18,    3.101
   },{ 5.013,   0.4707, 85.8,    16.55,    3.211
   },{ 4.32,    0.4947, 76.14,    10.85,    5.441
   },{ 4.652,   0.4571, 80.73,    22.0,    4.952
   },{ 3.114,   0.5236, 76.67,    7.62,    6.385
   },{ 3.114,   0.5236, 76.67,    7.62,    7.502
       // Z= 31-40
   },{ 3.114,   0.5236, 76.67,    7.62,    8.514
   },{ 5.746,   0.4662, 79.24,    1.185,    7.993
   },{ 2.792,   0.6346, 106.1,    0.2986,   2.331
   },{ 4.667,   0.5095, 124.3,    2.102,    1.667
   },{ 2.44,    0.6346, 105.0,    0.83,    2.851
   },{ 1.413,   0.7377, 147.9,    1.466,    1.016
   },{ 11.72,   0.3826, 102.8,    9.231,    4.371
   },{ 7.126,   0.4804, 119.3,    5.784,    2.454
   },{ 11.61,   0.3955, 146.7,    7.031,    1.423
   },{ 10.99,   0.41,   163.9,   7.1,      1.052
       // Z= 41-50
   },{ 9.241,   0.4275, 163.1,    7.954,    1.102
   },{ 9.276,   0.418,  157.1,   8.038,    1.29
   },{ 3.999,   0.6152, 97.6,    1.297,    5.792
   },{ 4.306,   0.5658, 97.99,    5.514,    5.754
   },{ 3.615,   0.6197, 86.26,    0.333,    8.689
   },{ 5.8,     0.49,   147.2,   6.903,    1.289
   },{ 5.6,     0.49,   130.0,   10.0,     2.844
   },{ 3.55,    0.6068, 124.7,    1.112,    3.119
   },{ 3.6,     0.62,   105.8,   0.1692,   6.026
   },{ 5.4,     0.53,   103.1,   3.931,    7.767
       // Z= 51-60
   },{ 3.97,    0.6459, 131.8,    0.2233,   2.723
   },{ 3.65,    0.64,   126.8,   0.6834,   3.411
   },{ 3.118,   0.6519, 164.9,    1.208,    1.51
   },{ 3.949,   0.6209, 200.5,    1.878,    0.9126
   },{ 14.4,    0.3923, 152.5,    8.354,    2.597
   },{ 10.99,   0.4599, 138.4,    4.811,    3.726
   },{ 16.6,    0.3773, 224.1,    6.28,    0.9121
   },{ 10.54,   0.4533, 159.3,   4.832,    2.529
   },{ 10.33,   0.4502, 162.0,   5.132,    2.444
   },{ 10.15,   0.4471, 165.6,   5.378,    2.328
       // Z= 61-70
   },{ 9.976,   0.4439, 168.0,   5.721,    2.258
   },{ 9.804,   0.4408, 176.2,   5.675,    1.997
   },{ 14.22,   0.363,  228.4,   7.024,    1.016
   },{ 9.952,   0.4318, 233.5,   5.065,    0.9244
   },{ 9.272,   0.4345, 210.0,   4.911,    1.258
   },{ 10.13,   0.4146, 225.7,   5.525,    1.055
   },{ 8.949,   0.4304, 213.3,   5.071,    1.221
   },{ 11.94,   0.3783, 247.2,   6.655,    0.849
   },{ 8.472,   0.4405, 195.5,   4.051,    1.604
   },{ 8.301,   0.4399, 203.7,   3.667,    1.459
       // Z= 71-80
   },{ 6.567,   0.4858, 193.0,   2.65,     1.66
   },{ 5.951,   0.5016, 196.1,   2.662,    1.589
   },{ 7.495,   0.4523, 251.4,   3.433,    0.8619
   },{ 6.335,   0.4825, 255.1,   2.834,    0.8228
   },{ 4.314,   0.5558, 214.8,   2.354,    1.263
   },{ 4.02,    0.5681, 219.9,   2.402,    1.191
   },{ 3.836,   0.5765, 210.2,   2.742,    1.305
   },{ 4.68,    0.5247, 244.7,   2.749,    0.8962
   },{ 2.892,   0.6204, 208.6,   2.415,    1.416 //Au Z77
       // },{ 3.223,   0.5883, 232.7,   2.954,    1.05  //Au ICRU
   },{ 2.892,   0.6204, 208.6,   2.415,    1.416
       // Z= 81-90
   },{ 4.728,   0.5522, 217.0,   3.091,    1.386
   },{ 6.18,    0.52,   170.0,   4.0,      3.224
   },{ 9.0,     0.47,   198.0,   3.8,      2.032
   },{ 2.324,   0.6997, 216.0,   1.599,    1.399
   },{ 1.961,   0.7286, 223.0,   1.621,    1.296
   },{ 1.75,    0.7427, 350.1,   0.9789,   0.5507
   },{ 10.31,   0.4613, 261.2,   4.738,    0.9899
   },{ 7.962,   0.519,  235.7,   4.347,    1.313
   },{ 6.227,   0.5645, 231.9,   3.961,    1.379
   },{ 5.246,   0.5947, 228.6,   4.027,    1.432
       // Z= 91-92
   },{ 5.408,   0.5811, 235.7,   3.961,    1.358
   },{ 5.218,   0.5828, 245.0,   3.838,    1.25}
  };

  // Free electron gas model
  if ( T < 0.001 ) {
    G4double slow  = a[i][0] ;
    G4double shigh = G4Log( 1.0 + a[i][3]*1000.0 + a[i][4]*0.001 )
                   * a[i][2]*1000.0 ;
    ionloss  = slow*shigh / (slow + shigh) ;
    ionloss *= sqrt(T*1000.0) ;

  // Main parametrisation
  } else {
    G4double slow  = a[i][0] * G4Exp(G4Log(T*1000.0)*a[i][1]) ;
    G4double shigh = G4Log( 1.0 + a[i][3]/T + a[i][4]*T ) * a[i][2]/T ;
    ionloss = slow*shigh / (slow + shigh) ;
    /*
    G4cout << "## " << i << ". T= " << T << " slow= " << slow
           << " a0= " << a[i][0] << " a1= " << a[i][1] 
           << " shigh= " << shigh 
           << " dedx= " << ionloss << " q^2= " <<  HeEffChargeSquare(z, T*MeV) 
	   << G4endl;
    */
  }
  if ( ionloss < 0.0) { ionloss = 0.0; }

  // He effective charge
  ionloss /= HeEffChargeSquare(z, T);

  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::DEDX(const G4Material* material,
                                     G4double kineticEnergy)
{
  G4double eloss = 0.0;
  // check DB
  if(material != currentMaterial) {
    currentMaterial = material;
    iASTAR    = -1;
    iMolecula = -1;
    if( !HasMaterial(material) ) { iASTAR = fASTAR->GetIndex(material); }
  }

  const G4int numberOfElements = material->GetNumberOfElements();
  const G4double* theAtomicNumDensityVector =
                                 material->GetAtomicNumDensityVector();

  if( iASTAR >= 0 ) {
    G4double T = kineticEnergy*rateMassHe2p;
    return fASTAR->GetElectronicDEDX(iASTAR, T)*material->GetDensity()/
      HeEffChargeSquare(fASTAR->GetEffectiveZ(iASTAR), T/MeV);

  } else if(iMolecula >= 0) {

    eloss = StoppingPower(material, kineticEnergy)*
      material->GetDensity()/amu;

  // pure material
  } else if(1 == numberOfElements) {

    G4double z = material->GetZ();
    eloss = ElectronicStoppingPower(z, kineticEnergy)
                               * (material->GetTotNbOfAtomsPerVolume());

  // Brugg's rule calculation
  } else {
    const G4ElementVector* theElementVector =
                           material->GetElementVector() ;

    //  loop for the elements in the material
    for (G4int i=0; i<numberOfElements; i++)
    {
      const G4Element* element = (*theElementVector)[i] ;
      eloss   += ElectronicStoppingPower(element->GetZ(), kineticEnergy)
                                   * theAtomicNumDensityVector[i];
    }
  }
  return eloss*theZieglerFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BraggIonModel::HeEffChargeSquare(G4double z, 
                                            G4double kinEnergyHeInMeV) const
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
  for (G4int i=1; i<6; i++) {
    y *= e ;
    x += y * c[i] ;
  }

  G4double w = 7.6 -  e ;
  w = 1.0 + (0.007 + 0.00005*z) * G4Exp( -w*w ) ;
  w = 4.0 * (1.0 - G4Exp(-x)) * w * w ;

  return w;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


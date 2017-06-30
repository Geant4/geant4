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
// $Id: G4PairProductionRelModel.cc 104555 2017-06-06 07:31:32Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PairProductionRelModel
//
// Author:        Andreas Schaelicke
//
// Creation date: 02.04.2009
//
// Modifications:
//
// 20.03.17    change LPMconstant such that it gives suppression variable 's'
//             that consistent to Migdal's one; fix a small bug in 'logTS1'
//             computation; suppression is consistent now with the one in the
//             brem. model (F.Hariri)
//
// Class Description:
//
// Main References:
//  J.W.Motz et.al., Rev. Mod. Phys. 41 (1969) 581.
//  S.Klein,  Rev. Mod. Phys. 71 (1999) 1501.
//  T.Stanev et.al., Phys. Rev. D25 (1982) 1291.
//  M.L.Ter-Mikaelian, High-energy Electromagnetic Processes in Condensed Media,
//                     Wiley, 1972.
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4PairProductionRelModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

const G4double G4PairProductionRelModel::facFel = G4Log(184.15);
const G4double G4PairProductionRelModel::facFinel = G4Log(1194.); // 1440.

const G4double G4PairProductionRelModel::preS1 = 1./(184.15*184.15);
const G4double G4PairProductionRelModel::logTwo = G4Log(2.);

const G4double G4PairProductionRelModel::xgi[]={ 0.0199, 0.1017, 0.2372, 0.4083,
						 0.5917, 0.7628, 0.8983, 0.9801 };
const G4double G4PairProductionRelModel::wgi[]={ 0.0506, 0.1112, 0.1569, 0.1813,
						 0.1813, 0.1569, 0.1112, 0.0506 };
const G4double G4PairProductionRelModel::Fel_light[]  = {0., 5.31  , 4.79  , 4.74 ,  4.71};
const G4double G4PairProductionRelModel::Finel_light[] = {0., 6.144 , 5.621 , 5.805 , 5.924};

const G4double G4PairProductionRelModel::xsfactor =
  4*CLHEP::fine_structure_const*CLHEP::classic_electr_radius*CLHEP::classic_electr_radius;
const G4double G4PairProductionRelModel::Egsmall = 2.*CLHEP::MeV;
const G4double G4PairProductionRelModel::Eghigh = 100.*CLHEP::GeV;

G4PairProductionRelModel::G4PairProductionRelModel(const G4ParticleDefinition*,
						   const G4String& nam)
  : G4VEmModel(nam),
    fLPMconstant(CLHEP::fine_structure_const*CLHEP::electron_mass_c2*CLHEP::electron_mass_c2/
		 (4.*CLHEP::pi*CLHEP::hbarc)),
    fLPMflag(true),
    lpmEnergy(0.),
    use_completescreening(false)
{
  fParticleChange = nullptr;
  theGamma    = G4Gamma::Gamma();
  thePositron = G4Positron::Positron();
  theElectron = G4Electron::Electron();

  nist = G4NistManager::Instance();

  currentZ = z13 = z23 = lnZ = Fel = Finel = fCoulomb = phiLPM = gLPM = xiLPM = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PairProductionRelModel::~G4PairProductionRelModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PairProductionRelModel::Initialise(const G4ParticleDefinition* p,
					  const G4DataVector& cuts)
{
  if(!fParticleChange) { fParticleChange = GetParticleChangeForGamma(); }
  if(IsMaster() && LowEnergyLimit() < HighEnergyLimit()) {
    InitialiseElementSelectors(p, cuts);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PairProductionRelModel::InitialiseLocal(const G4ParticleDefinition*,
					       G4VEmModel* masterModel)
{
  if(LowEnergyLimit() < HighEnergyLimit()) {
    SetElementSelectors(masterModel->GetElementSelectors());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PairProductionRelModel::ComputeXSectionPerAtom(G4double totalEnergy, G4double Z)
{
  G4double cross = 0.0;

  // number of intervals and integration step
  G4double vcut = electron_mass_c2/totalEnergy ;

  // limits by the screening variable
  G4double dmax = DeltaMax();
  G4double dmin = std::min(DeltaMin(totalEnergy),dmax);
  G4double vcut1 = 0.5 - 0.5*sqrt(1. - dmin/dmax);
  vcut = max(vcut, vcut1);

  G4double vmax = 0.5;
  G4int n = 1;  // needs optimisation

  G4double delta = (vmax - vcut)*totalEnergy/G4double(n);

  G4double e0 = vcut*totalEnergy;

  // simple integration
  for(G4int l=0; l<n; ++l) {
    e0 += delta;
    for(G4int i=0; i<8; ++i) {

      G4double eg = (e0 + xgi[i]*delta);
      G4double xs = (fLPMflag && totalEnergy > Eghigh)
	? ComputeRelDXSectionPerAtom(eg,totalEnergy,Z)
	: ComputeDXSectionPerAtom(eg,totalEnergy,Z);
      cross += wgi[i]*xs;
    }
  }

  cross *= delta*2.;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4PairProductionRelModel::ComputeDXSectionPerAtom(G4double eplusEnergy,
						  G4double totalEnergy,
						  G4double /*Z*/)
{
  // most simple case - complete screening:

  //  dsig/dE+ = 4 * alpha * Z**2 * r0**2 / k
  //     * [ (y**2 + (1-y**2) + 2/3*y*(1-y) ) * ( log (183 * Z**-1/3)  + 1/9 * y*(1-y) ]
  // y = E+/k
  G4double yp=eplusEnergy/totalEnergy;
  G4double ym=1.-yp;

  G4double cross = 0.;
  if (use_completescreening)
    cross = (yp*yp + ym*ym + 2./3.*ym*yp)*(Fel - fCoulomb) + yp*ym/9.;
  else {
    G4double delta = 0.25*DeltaMin(totalEnergy)/(yp*ym);
    cross = (yp*yp + ym*ym)*(0.25*Phi1(delta) - lnZ/3. - fCoulomb)
      + 2./3.*ym*yp*(0.25*Phi2(delta) - lnZ/3. - fCoulomb);
  }
  return cross/totalEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4PairProductionRelModel::ComputeRelDXSectionPerAtom(G4double eplusEnergy,
						     G4double totalEnergy,
						     G4double /*Z*/)
{
  // most simple case - complete screening:

  //  dsig/dE+ = 4 * alpha * Z**2 * r0**2 / k
  //     * [ (y**2 + (1-y**2) + 2/3*y*(1-y) ) * ( log (183 * Z**-1/3)  + 1/9 * y*(1-y) ]
  // y = E+/k
  G4double yp=eplusEnergy/totalEnergy;
  G4double ym=1.-yp;

  CalcLPMFunctions(totalEnergy,eplusEnergy); // gamma

  G4double cross = 0.;
  if (use_completescreening)
    cross = xiLPM*(2./3.*phiLPM*(yp*yp + ym*ym) + gLPM)*(Fel - fCoulomb);
  else {
    G4double delta = 0.25*DeltaMin(totalEnergy)/(yp*ym);
    cross = (1./3.*gLPM + 2./3.*phiLPM)*(yp*yp + ym*ym)
                             *(0.25*Phi1(delta) - lnZ/3. - fCoulomb)
           + 2./3.*gLPM*ym*yp*(0.25*Phi2(delta) - lnZ/3. - fCoulomb);
    cross *= xiLPM;
  }
  return cross/totalEnergy;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4PairProductionRelModel::CalcLPMFunctions(G4double k, G4double eplusEnergy)
{
  // *** calculate lpm variable s & sprime ***
  // Klein eqs. (78) & (79)
  G4double sprime = sqrt(0.125*k*lpmEnergy/(eplusEnergy*(k-eplusEnergy)));

  G4double s1 = preS1*z23;
  G4double logS1 = 2./3.*lnZ-2.*facFel;
  G4double logTS1 = 0.5*logTwo+logS1;

  xiLPM = 2.;

  if (sprime>1)
    xiLPM = 1.;
  else if (sprime>sqrt(2.)*s1) {
    G4double h  = G4Log(sprime)/logTS1;
    xiLPM = 1+h-0.08*(1-h)*(1-sqr(1-h))/logTS1;
  }

  G4double s0 = sprime/sqrt(xiLPM);
  //   G4cout<<"k="<<k<<" y="<<eplusEnergy/k<<G4endl;
  //   G4cout<<"s0="<<s0<<G4endl;

  // *** calculate supression functions phi and G ***
  // Klein eqs. (77)
  G4double s2=s0*s0;
  G4double s3=s0*s2;
  G4double s4=s2*s2;

  if (s0<0.1) {
    // high suppression limit
    phiLPM = 6.*s0 - 18.84955592153876*s2 + 39.47841760435743*s3
      - 57.69873135166053*s4;
    gLPM = 37.69911184307752*s2 - 236.8705056261446*s3 + 807.7822389*s4;
  }
  else if (s0<1.9516) {
    // intermediate suppression
    // using eq.77 approxim. valid s0<2.
    phiLPM = 1.-G4Exp(-6.*s0*(1.+(3.-pi)*s0)
		+s3/(0.623+0.795*s0+0.658*s2));
    if (s0<0.415827397755) {
      // using eq.77 approxim. valid 0.07<s<2
      G4double psiLPM = 1-G4Exp(-4*s0-8*s2/(1+3.936*s0+4.97*s2-0.05*s3+7.50*s4));
      gLPM = 3*psiLPM-2*phiLPM;
    }
    else {
      // using alternative parametrisiation
      G4double pre = -0.16072300849123999 + s0*3.7550300067531581 + s2*-1.7981383069010097
	+ s3*0.67282686077812381 + s4*-0.1207722909879257;
      gLPM = std::tanh(pre);
    }
  }
  else {
    // low suppression limit valid s>2.
    phiLPM = 1. - 0.0119048/s4;
    gLPM = 1. - 0.0230655/s4;
  }

  // *** make sure suppression is smaller than 1 ***
  // *** caused by Migdal approximation in xi    ***
  if (xiLPM*phiLPM>1. || s0>0.57)  { xiLPM=1./phiLPM; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4PairProductionRelModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
      G4double gammaEnergy, G4double Z, G4double, G4double, G4double)
{
  G4double crossSection = 0.0 ;
  if ( gammaEnergy <= 2.0*electron_mass_c2 ) { return crossSection; }

  SetCurrentElement(Z);
  // choose calculator according to parameters and switches
  // in the moment only one calculator:
  crossSection=ComputeXSectionPerAtom(gammaEnergy,Z);

  G4double xi = Finel/(Fel - fCoulomb); // inelastic contribution
  crossSection *= xsfactor*Z*(Z+xi);

  return crossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4PairProductionRelModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					    const G4MaterialCutsCouple* couple,
					    const G4DynamicParticle* aDynamicGamma,
					    G4double,
					    G4double)
// The secondaries e+e- energies are sampled using the Bethe - Heitler
// cross sections with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at
//          low energy are ignored.
// Note 2 : The differential cross section implicitly takes account of
//          pair creation in both nuclear and atomic electron fields.
//          However triplet prodution is not generated.
{
  const G4Material* aMaterial = couple->GetMaterial();

  G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
  G4ParticleMomentum GammaDirection = aDynamicGamma->GetMomentumDirection();

  G4double epsil ;
  G4double epsil0 = electron_mass_c2/GammaEnergy ;
  if(epsil0 > 1.0) { return; }

  SetupForMaterial(theGamma, aMaterial, GammaEnergy);

  // select randomly one element constituing the material
  const G4Element* anElement =
    SelectRandomAtom(aMaterial, theGamma, GammaEnergy);

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  // do it fast if GammaEnergy < 2. MeV
  if (GammaEnergy < Egsmall) {
    epsil = epsil0 + (0.5-epsil0)*rndmEngine->flat();

  } else {
    // now comes the case with GammaEnergy >= 2. MeV
    // Extract Coulomb factor for this Element
    G4double FZ = 8.*(anElement->GetIonisation()->GetlogZ3());
    static const G4double midEnergy = 50.*CLHEP::MeV;
    if (GammaEnergy > midEnergy) { FZ += 8.*(anElement->GetfCoulomb()); }

    // limits of the screening variable
    G4double screenfac = 136.*epsil0/(anElement->GetIonisation()->GetZ3());
    //F.Hariri : correct sign of last term
    G4double screenmax = G4Exp ((42.24 - FZ)/8.368) + 0.952 ;
    G4double screenmin = std::min(4.*screenfac, screenmax);

    // limits of the energy sampling
    G4double epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
    G4double epsilmin = std::max(epsil0, epsil1); 
    G4double epsilrange = 0.5 - epsilmin;

    //
    // sample the energy rate of the created electron (or positron)
    //
    //G4double epsil, screenvar, greject ;
    G4double  screenvar, greject ;

    G4double F10 = ScreenFunction1(screenmin) - FZ;
    G4double F20 = ScreenFunction2(screenmin) - FZ;
    G4double NormF1 = std::max(F10*epsilrange*epsilrange,0.);
    G4double NormF2 = std::max(1.5*F20,0.);

    do {
      if ( NormF1/(NormF1+NormF2) > rndmEngine->flat() ) {
	epsil = 0.5 - epsilrange*nist->GetZ13(rndmEngine->flat());
	screenvar = screenfac/(epsil*(1-epsil));
	if (fLPMflag && GammaEnergy > Eghigh) {
	  CalcLPMFunctions(GammaEnergy,GammaEnergy*epsil);
	  greject = xiLPM*((gLPM+2.*phiLPM)*Phi1(screenvar) -
			   gLPM*Phi2(screenvar) - phiLPM*FZ)/F10;
	}
	else {
	  greject = (ScreenFunction1(screenvar) - FZ)/F10;
	}

      } else {
	epsil = epsilmin + epsilrange*rndmEngine->flat();
	screenvar = screenfac/(epsil*(1-epsil));
	if (fLPMflag && GammaEnergy > Eghigh) {
	  CalcLPMFunctions(GammaEnergy,GammaEnergy*epsil);
	  greject = xiLPM*((0.5*gLPM+phiLPM)*Phi1(screenvar) +
			   0.5*gLPM*Phi2(screenvar) - 0.5*(gLPM+phiLPM)*FZ)/F20;
	}
	else {
	  greject = (ScreenFunction2(screenvar) - FZ)/F20;
	}
      }

      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while( greject < rndmEngine->flat());

  }   //  end of epsil sampling

  //
  // fixe charges randomly
  //

  G4double ElectTotEnergy, PositTotEnergy;
  if (rndmEngine->flat() > 0.5) {
    ElectTotEnergy = (1.-epsil)*GammaEnergy;
    PositTotEnergy = epsil*GammaEnergy;

  } else {
    PositTotEnergy = (1.-epsil)*GammaEnergy;
    ElectTotEnergy = epsil*GammaEnergy;
  }

  //
  // scattered electron (positron) angles. ( Z - axis along the parent photon)
  //
  //  universal distribution suggested by L. Urban
  // (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  static const G4double a1 = 1.6;
  static const G4double a2 = a1/3.;
  G4double uu = -G4Log(rndmEngine->flat()*rndmEngine->flat());
  G4double u = (0.25 > rndmEngine->flat()) ? uu*a1 : uu*a2;

  G4double thetaEle = u*electron_mass_c2/ElectTotEnergy;
  G4double sinte = std::sin(thetaEle);
  G4double coste = std::cos(thetaEle);

  G4double thetaPos = u*electron_mass_c2/PositTotEnergy;
  G4double sintp = std::sin(thetaPos);
  G4double costp = std::cos(thetaPos);

  G4double phi  = twopi * rndmEngine->flat();
  G4double sinp = std::sin(phi);
  G4double cosp = std::cos(phi);

  //
  // kinematic of the created pair
  //
  // the electron and positron are assumed to have a symetric
  // angular distribution with respect to the Z axis along the parent photon.

  G4double ElectKineEnergy = max(0.,ElectTotEnergy - electron_mass_c2);

  G4ThreeVector ElectDirection (sinte*cosp, sinte*sinp, coste);
  ElectDirection.rotateUz(GammaDirection);

  // create G4DynamicParticle object for the particle1
  G4DynamicParticle* aParticle1= new G4DynamicParticle(
		     theElectron,ElectDirection,ElectKineEnergy);

  // the e+ is always created (even with Ekine=0) for further annihilation.

  G4double PositKineEnergy = std::max(0.,PositTotEnergy - electron_mass_c2);

  G4ThreeVector PositDirection (-sintp*cosp, -sintp*sinp, costp);
  PositDirection.rotateUz(GammaDirection);

  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2= new G4DynamicParticle(
                      thePositron,PositDirection,PositKineEnergy);

  // Fill output vector
  fvect->push_back(aParticle1);
  fvect->push_back(aParticle2);

  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4PairProductionRelModel::SetupForMaterial(const G4ParticleDefinition*,
						const G4Material* mat, G4double)
{
  lpmEnergy = mat->GetRadlen()*fLPMconstant;
  //  G4cout<<" lpmEnergy="<<lpmEnergy<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

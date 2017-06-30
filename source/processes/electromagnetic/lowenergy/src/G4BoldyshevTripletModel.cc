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
// Author: Sebastien Incerti
//         22 January 2012
//         on base of G4BoldyshevTripletModel (original version)
//         and G4LivermoreRayleighModel (MT version)

#include "G4BoldyshevTripletModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4BoldyshevTripletModel::maxZ = 99;
G4LPhysicsFreeVector* G4BoldyshevTripletModel::data[] = {0};

G4BoldyshevTripletModel::G4BoldyshevTripletModel(const G4ParticleDefinition*, const G4String& nam)
  :G4VEmModel(nam),smallEnergy(4.*MeV)
{
  fParticleChange = nullptr;
  
  lowEnergyLimit = 4.0*electron_mass_c2;
  momentumThreshold_c = energyThreshold = xb = xn = lowEnergyLimit; 
  
  verboseLevel= 0;
  // Verbosity scale for debugging purposes:
  // 0 = nothing 
  // 1 = calculation of cross sections, file openings...
  // 2 = entering in methods

  if(verboseLevel > 0) 
  {
    G4cout << "G4BoldyshevTripletModel is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BoldyshevTripletModel::~G4BoldyshevTripletModel()
{
  if(IsMaster()) {
    for(G4int i=0; i<maxZ; ++i) {
      if(data[i]) { 
	delete data[i];
	data[i] = nullptr;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BoldyshevTripletModel::Initialise(const G4ParticleDefinition*,
				         const G4DataVector&)
{
  if (verboseLevel > 1) 
  {
    G4cout << "Calling Initialise() of G4BoldyshevTripletModel." 
	   << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / MeV << " MeV - "
	   << HighEnergyLimit() / GeV << " GeV isMaster: " << IsMaster()
	   << G4endl;
  }
  // compute values only once
  energyThreshold = 1.1*electron_mass_c2; 
  momentumThreshold_c = std::sqrt(energyThreshold * energyThreshold 
                                - electron_mass_c2*electron_mass_c2); 
  G4double momentumThreshold_N = momentumThreshold_c/electron_mass_c2; 
  G4double t = 0.5*G4Log(momentumThreshold_N + 
			 std::sqrt(momentumThreshold_N*momentumThreshold_N + 1.0));
  //G4cout << 0.5*asinh(momentumThreshold_N) << "  " << t << G4endl;
  G4double sinht = std::sinh(t);                         
  G4double cosht = std::cosh(t);  
  G4double logsinht = G4Log(2.*sinht);                     
  G4double J1 = 0.5*(t*cosht/sinht - logsinht);
  G4double J2 = (-2./3.)*logsinht + t*cosht/sinht 
    + (sinht - t*cosht*cosht*cosht)/(3.*sinht*sinht*sinht);

  xb = 2.*(J1-J2)/J1; 
  xn = 1. - xb/6.;

  if(IsMaster()) 
  {
    // Access to elements  
    char* path = getenv("G4LEDATA");

    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
  
    G4int numOfCouples = theCoupleTable->GetTableSize();
  
    for(G4int i=0; i<numOfCouples; ++i) 
    {
      const G4Material* material = 
        theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
      const G4ElementVector* theElementVector = material->GetElementVector();
      G4int nelm = material->GetNumberOfElements();
    
      for (G4int j=0; j<nelm; ++j) 
      {
        G4int Z = std::min((*theElementVector)[j]->GetZasInt(), maxZ);
        if(!data[Z]) { ReadData(Z, path); }
      }
    }
  }
  if(!fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4BoldyshevTripletModel::MinPrimaryEnergy(const G4Material*,
					  const G4ParticleDefinition*,
					  G4double)
{
  return lowEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BoldyshevTripletModel::ReadData(size_t Z, const char* path)
{
  if (verboseLevel > 1) 
  {
    G4cout << "Calling ReadData() of G4BoldyshevTripletModel" 
	   << G4endl;
  }

  if(data[Z]) { return; }
  
  const char* datadir = path;

  if(!datadir) 
  {
    datadir = getenv("G4LEDATA");
    if(!datadir) 
    {
      G4Exception("G4BoldyshevTripletModel::ReadData()",
		  "em0006",FatalException,
		  "Environment variable G4LEDATA not defined");
      return;
    }
  }
  
  data[Z] = new G4LPhysicsFreeVector();
  std::ostringstream ost;
  ost << datadir << "/livermore/tripdata/pp-trip-cs-" << Z <<".dat";
  std::ifstream fin(ost.str().c_str());
  
  if( !fin.is_open()) 
  {
    G4ExceptionDescription ed;
    ed << "G4BoldyshevTripletModel data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4BoldyshevTripletModel::ReadData()",
		"em0003",FatalException,
		ed,"G4LEDATA version should be G4EMLOW6.27 or later.");
    return;
  } 
  
  else 
  {
    
    if(verboseLevel > 3) { G4cout << "File " << ost.str() 
	     << " is opened by G4BoldyshevTripletModel" << G4endl;}
    
    data[Z]->Retrieve(fin, true);
  } 

  // Activation of spline interpolation
  data[Z]->SetSpline(true);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BoldyshevTripletModel::ComputeCrossSectionPerAtom(
         const G4ParticleDefinition* part,
	 G4double GammaEnergy, G4double Z, G4double, G4double, G4double)
{
  if (verboseLevel > 1) 
  {
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4BoldyshevTripletModel" 
	   << G4endl;
  }

  if (GammaEnergy < lowEnergyLimit) { return 0.0; } 

  G4double xs = 0.0;  
  G4int intZ = std::max(1, std::min(G4lrint(Z), maxZ));
  G4LPhysicsFreeVector* pv = data[intZ];

  // if element was not initialised
  // do initialisation safely for MT mode
  if(!pv) 
  {
    InitialiseForElement(part, intZ);
    pv = data[intZ];
    if(!pv) { return xs; }
  }
  // x-section is taken from the table
  xs = pv->Value(GammaEnergy); 

  if(verboseLevel > 1)
  {
    G4cout  <<  "*** Triplet conversion xs for Z=" << Z << " at energy E(MeV)=" 
	    << GammaEnergy/MeV <<  "  cs=" << xs/millibarn << " mb" << G4endl;
  }
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BoldyshevTripletModel::SampleSecondaries(
                                 std::vector<G4DynamicParticle*>* fvect,
				 const G4MaterialCutsCouple* /*couple*/,
				 const G4DynamicParticle* aDynamicGamma,
				 G4double, G4double)
{

  // The energies of the secondary particles are sampled using
  // a modified Wheeler-Lamb model (see PhysRevD 7 (1973), 26)
  if (verboseLevel > 1) {
    G4cout << "Calling SampleSecondaries() of G4BoldyshevTripletModel" 
	   << G4endl;
  }

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();
  G4ParticleMomentum photonDirection = aDynamicGamma->GetMomentumDirection();

  G4double epsilon;

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  // recoil electron thould be 3d particle
  G4DynamicParticle* particle3 = nullptr;
  static const G4double costlim = std::cos(4.47*CLHEP::pi/180.);
  
  G4double loga, f1_re, greject, cost;
  G4double cosThetaMax = (energyThreshold - electron_mass_c2 
     + electron_mass_c2*(energyThreshold + electron_mass_c2)/photonEnergy )
	/momentumThreshold_c;
  if (cosThetaMax > 1.) {
    //G4cout << "G4BoldyshevTripletModel::SampleSecondaries: ERROR cosThetaMax= " 
    //	   << cosThetaMax << G4endl;
    cosThetaMax = 1.0;
  }

  G4double logcostm = G4Log(cosThetaMax);
  G4int nn = 0;
  do {
    cost = G4Exp(logcostm*rndmEngine->flat());
    G4double are = 1./(14.*cost*cost);
    G4double bre = (1.-5.*cost*cost)/(2.*cost);
    loga = G4Log((1.+ cost)/(1.- cost));
    f1_re = 1. - bre*loga;
    greject = (cost < costlim) ? are*f1_re : 1.0;
    // G4cout << nn << ". step of the 1st loop greject= " << greject << G4endl;
    ++nn;
  } while(greject < rndmEngine->flat());
      
  // Calculo de phi - elecron de recoil
  G4double sint2 = (1. - cost)*(1. + cost);
  G4double fp = 1. - sint2*loga/(2.*cost) ;
  G4double rt, phi_re;
  nn = 0;
  do {
    phi_re = twopi*rndmEngine->flat();
    rt = (1. - std::cos(2.*phi_re)*fp/f1_re)/twopi;
    //G4cout << nn << ". step of the 2nd loop greject= " << rt << G4endl;
    ++nn;
  } while(rt < rndmEngine->flat());

  // Calculo de la energia - elecron de recoil - relacion momento maximo <-> angulo
  G4double S  = electron_mass_c2*(2.* photonEnergy + electron_mass_c2);
  G4double P2 = S - electron_mass_c2*electron_mass_c2;

  G4double D2 = 4.*S * electron_mass_c2*electron_mass_c2 + P2*P2*sint2;
  G4double ener_re = electron_mass_c2 * (S + electron_mass_c2*electron_mass_c2)/sqrt(D2);
      
  if(ener_re >= energyThreshold) 
    {
      G4double electronRKineEnergy = ener_re - electron_mass_c2;
      G4double sint = std::sqrt(sint2);
      G4ThreeVector electronRDirection (sint*std::cos(phi_re), sint*std::sin(phi_re), cost);
      electronRDirection.rotateUz(photonDirection);
      particle3 = new G4DynamicParticle (G4Electron::Electron(),
					 electronRDirection,
					 electronRKineEnergy);
    }
  else
    {
      // deposito la energia  ener_re - electron_mass_c2
      // G4cout << "electron de retroceso " << ener_re << G4endl;
      fParticleChange->ProposeLocalEnergyDeposit(std::max(0.0, ener_re - electron_mass_c2));
      ener_re = 0.0;
    }
  
  // Depaola (2004) suggested distribution for e+e- energy
  // VI: very suspect that 1 random number is not enough
  //     and sampling below is not correct - should be fixed
  G4double re = rndmEngine->flat();
  
  G4double a  = std::sqrt(16./xb - 3. - 36.*re*xn + 36.*re*re*xn*xn + 6.*xb*re*xn);
  G4double c1 = G4Exp(G4Log((-6. + 12.*re*xn + xb + 2*a)*xb*xb)/3.);
  epsilon = c1/(2.*xb) + (xb - 4.)/(2.*c1) + 0.5;
  
  G4double photonEnergy1 = photonEnergy - ener_re ; 
  // resto al foton la energia del electron de retro.
  G4double positronTotEnergy = std::max(epsilon*photonEnergy1, electron_mass_c2);
  G4double electronTotEnergy = std::max(photonEnergy1 - positronTotEnergy, electron_mass_c2);
  
  static const G4double a1 = 1.6;
  static const G4double a2 = 0.5333333333;
  G4double uu = -G4Log(rndmEngine->flat()*rndmEngine->flat());
  G4double u = (0.25 > rndmEngine->flat()) ? uu*a1 : uu*a2;

  G4double thetaEle = u*electron_mass_c2/electronTotEnergy;
  G4double sinte = std::sin(thetaEle);
  G4double coste = std::cos(thetaEle);

  G4double thetaPos = u*electron_mass_c2/positronTotEnergy;
  G4double sintp = std::sin(thetaPos);
  G4double costp = std::cos(thetaPos);

  G4double phi  = twopi * rndmEngine->flat();
  G4double sinp = std::sin(phi);
  G4double cosp = std::cos(phi);

  // Kinematics of the created pair:
  // the electron and positron are assumed to have a symetric angular
  // distribution with respect to the Z axis along the parent photon

  G4double electronKineEnergy = electronTotEnergy - electron_mass_c2;

  G4ThreeVector electronDirection (sinte*cosp, sinte*sinp, coste);
  electronDirection.rotateUz(photonDirection);

  G4DynamicParticle* particle1 = new G4DynamicParticle (G4Electron::Electron(),
                                                        electronDirection,
							electronKineEnergy);

  G4double positronKineEnergy = positronTotEnergy - electron_mass_c2;

  G4ThreeVector positronDirection (-sintp*cosp, -sintp*sinp, costp);
  positronDirection.rotateUz(photonDirection);

  // Create G4DynamicParticle object for the particle2      
  G4DynamicParticle* particle2 = new G4DynamicParticle(G4Positron::Positron(),
                                                       positronDirection, positronKineEnergy);
  // Fill output vector       
  
  fvect->push_back(particle1);
  fvect->push_back(particle2);

  if(particle3) { fvect->push_back(particle3); }
  
  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AutoLock.hh"
namespace { G4Mutex BoldyshevTripletModelMutex = G4MUTEX_INITIALIZER; }

void G4BoldyshevTripletModel::InitialiseForElement(
     const G4ParticleDefinition*, G4int Z)
{
  G4AutoLock l(&BoldyshevTripletModelMutex);
  // G4cout << "G4BoldyshevTripletModel::InitialiseForElement Z= " 
  //	  << Z << G4endl;
  if(!data[Z]) { ReadData(Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

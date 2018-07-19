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
//         on base of G4LivermoreGammaConversionModel (original version)
//         and G4LivermoreRayleighModel (MT version)

#include "G4LivermoreGammaConversionModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4EmParameters.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreGammaConversionModel::lowEnergyLimit = 2.*CLHEP::electron_mass_c2;
G4double G4LivermoreGammaConversionModel::tripletLowEnergy = 0.0;
G4double G4LivermoreGammaConversionModel::tripletHighEnergy = 100.0*CLHEP::GeV;
G4int G4LivermoreGammaConversionModel::verboseLevel = 0;
G4int G4LivermoreGammaConversionModel::nbinsTriplet = 0;
G4int G4LivermoreGammaConversionModel::maxZ = 99;
G4LPhysicsFreeVector* G4LivermoreGammaConversionModel::data[] = {nullptr};
G4PhysicsLogVector* G4LivermoreGammaConversionModel::probTriplet[] = {nullptr};

G4LivermoreGammaConversionModel::G4LivermoreGammaConversionModel
(const G4ParticleDefinition*, const G4String& nam)
  : G4VEmModel(nam),fParticleChange(nullptr)
{
  // Verbosity scale for debugging purposes:
  // 0 = nothing 
  // 1 = calculation of cross sections, file openings...
  // 2 = entering in methods

  if(verboseLevel > 0) 
  {
    G4cout << "G4LivermoreGammaConversionModel is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreGammaConversionModel::~G4LivermoreGammaConversionModel()
{
  if(IsMaster()) {
    for(G4int i=0; i<maxZ; ++i) {
      if(data[i]) { 
	delete data[i];
	data[i] = nullptr;
      }
      if(probTriplet[i]) { 
	delete probTriplet[i];
	probTriplet[i] = nullptr;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreGammaConversionModel::Initialise(
                                const G4ParticleDefinition* particle,
				const G4DataVector& cuts)
{
  if (verboseLevel > 1) 
  {
    G4cout << "Calling Initialise() of G4LivermoreGammaConversionModel." 
	   << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / MeV << " MeV - "
	   << HighEnergyLimit() / GeV << " GeV isMater: " << IsMaster() 
	   << G4endl;
  }

  if(!fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
    if(GetTripletModel()) {
      GetTripletModel()->SetParticleChange(fParticleChange);
    }
  }
  if(GetTripletModel()) { GetTripletModel()->Initialise(particle, cuts); }

  if(IsMaster()) 
  {
    // Initialise element selector
    InitialiseElementSelectors(particle, cuts);

    // Access to elements
    char* path = getenv("G4LEDATA");

    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
  
    G4int numOfCouples = theCoupleTable->GetTableSize();
  
    for(G4int i=0; i<numOfCouples; ++i) 
    {
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      SetCurrentCouple(couple);
      const G4Material* mat = couple->GetMaterial();
      const G4ElementVector* theElementVector = mat->GetElementVector();
      G4int nelm = mat->GetNumberOfElements();
    
      for (G4int j=0; j<nelm; ++j) 
      {
        G4int Z = std::min((*theElementVector)[j]->GetZasInt(), maxZ);
        if(!data[Z]) { ReadData(Z, path); }
        if(GetTripletModel()) { InitialiseProbability(particle, Z); }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreGammaConversionModel::InitialiseLocal(
     const G4ParticleDefinition*, G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4LivermoreGammaConversionModel::MinPrimaryEnergy(const G4Material*,
						  const G4ParticleDefinition*,
						  G4double)
{
  return lowEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreGammaConversionModel::ReadData(size_t Z, const char* path)
{
  if (verboseLevel > 1) 
  {
    G4cout << "Calling ReadData() of G4LivermoreGammaConversionModel" 
	   << G4endl;
  }

  if(data[Z]) { return; }
  
  const char* datadir = path;

  if(!datadir) 
  {
    datadir = getenv("G4LEDATA");
    if(!datadir) 
    {
      G4Exception("G4LivermoreGammaConversionModel::ReadData()",
		  "em0006",FatalException,
		  "Environment variable G4LEDATA not defined");
      return;
    }
  }
  data[Z] = new G4LPhysicsFreeVector();
  std::ostringstream ost;
  ost << datadir << "/livermore/pair/pp-cs-" << Z <<".dat";
  std::ifstream fin(ost.str().c_str());
  
  if( !fin.is_open()) 
  {
    G4ExceptionDescription ed;
    ed << "G4LivermoreGammaConversionModel data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4LivermoreGammaConversionModel::ReadData()",
		"em0003",FatalException,
		ed,"G4LEDATA version should be G4EMLOW6.27 or later.");
    return;
  }   
  else 
  {
    
    if(verboseLevel > 1) { G4cout << "File " << ost.str() 
	     << " is opened by G4LivermoreGammaConversionModel" << G4endl;}
    
    data[Z]->Retrieve(fin, true);
  } 
  // Activation of spline interpolation
  data[Z] ->SetSpline(true);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreGammaConversionModel::ComputeCrossSectionPerAtom(
           const G4ParticleDefinition* particle,
	   G4double GammaEnergy, G4double Z, G4double, G4double, G4double)
{
  if (verboseLevel > 1) 
  {
    G4cout << "G4LivermoreGammaConversionModel::ComputeCrossSectionPerAtom() Z= " 
	   << Z << G4endl;
  }

  if (GammaEnergy < lowEnergyLimit) { return 0.0; } 

  G4double xs = 0.0;
  
  G4int intZ = std::max(1, std::min(G4lrint(Z), maxZ));

  G4LPhysicsFreeVector* pv = data[intZ];

  // if element was not initialised
  // do initialisation safely for MT mode
  if(!pv) 
  {
    InitialiseForElement(particle, intZ);
    pv = data[intZ];
    if(!pv) { return xs; }
  }
  // x-section is taken from the table
  xs = pv->Value(GammaEnergy); 

  if(verboseLevel > 0)
  {
    G4cout  <<  "*** Gamma conversion xs for Z=" << Z << " at energy E(MeV)=" 
	    << GammaEnergy/MeV <<  "  cs=" << xs/millibarn << " mb" << G4endl;
  }

  return xs;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreGammaConversionModel::SampleSecondaries(
                                 std::vector<G4DynamicParticle*>* fvect,
				 const G4MaterialCutsCouple* couple,
				 const G4DynamicParticle* aDynamicGamma,
				 G4double, G4double)
{

  // The energies of the e+ e- secondaries are sampled using the Bethe - Heitler
  // cross sections with Coulomb correction. A modified version of the random
  // number techniques of Butcher & Messel is used (Nuc Phys 20(1960),15).
  
  // Note 1 : Effects due to the breakdown of the Born approximation at low
  // energy are ignored.
  // Note 2 : The differential cross section implicitly takes account of
  // pair creation in both nuclear and atomic electron fields. However triplet
  // prodution is not generated.

  if (verboseLevel > 1) {
    G4cout << "Calling SampleSecondaries() of G4LivermoreGammaConversionModel" 
	   << G4endl;
  }

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();
  G4ParticleMomentum photonDirection = aDynamicGamma->GetMomentumDirection();

  G4double epsilon ;
  G4double epsilon0Local = electron_mass_c2 / photonEnergy ;

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  // Do it fast if photon energy < 2. MeV
  static const G4double smallEnergy = 2.*CLHEP::MeV;
  if (photonEnergy < smallEnergy )
  {
    epsilon = epsilon0Local + (0.5 - epsilon0Local) * rndmEngine->flat();
  }
  else
  {
    // Select randomly one element in the current material

    const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
    const G4Element* element = SelectRandomAtom(couple,particle,photonEnergy);
    G4int Z = element->GetZasInt();

    // triplet production
    if(GetTripletModel()) {
      if(!probTriplet[Z]) { InitialiseForElement(particle, Z); }
      /*
      G4cout << "Liv: E= " << photonEnergy 
	     << " prob= " << probTriplet[Z]->Value(photonEnergy) 
	     << G4endl;
      */
      if(probTriplet[Z] && 
	 rndmEngine->flat() < probTriplet[Z]->Value(photonEnergy)) { 
	GetTripletModel()->SampleSecondaries(fvect, couple, aDynamicGamma);
	return;
      }
    }

    G4IonisParamElm* ionisation = element->GetIonisation();

    // Extract Coulomb factor for this Element
    G4double fZ = 8. * (ionisation->GetlogZ3());
    static const G4double midEnergy = 50.*CLHEP::MeV;
    if (photonEnergy > midEnergy) { fZ += 8. * (element->GetfCoulomb()); }

    // Limits of the screening variable
    G4double screenFactor = 136. * epsilon0Local / (element->GetIonisation()->GetZ3());
    G4double screenMax = G4Exp((42.24 - fZ)/8.368) + 0.952;
    G4double screenMin = std::min(4.*screenFactor,screenMax);

    // Limits of the energy sampling
    G4double epsilon1 = 0.5 - 0.5 * std::sqrt(1. - screenMin / screenMax) ;
    G4double epsilonMin = std::max(epsilon0Local,epsilon1);
    G4double epsilonRange = 0.5 - epsilonMin ;

    // Sample the energy rate of the created electron (or positron)
    G4double screen;
    G4double gReject;

    G4double f10 = ScreenFunction1(screenMin) - fZ;
    G4double f20 = ScreenFunction2(screenMin) - fZ;
    G4double normF1 = std::max(f10 * epsilonRange * epsilonRange,0.);
    G4double normF2 = std::max(1.5 * f20,0.);

    do 
      {
	if (normF1 > (normF1 + normF2)*rndmEngine->flat() )
	  {
	    epsilon = 0.5 - epsilonRange *G4Exp(G4Log(rndmEngine->flat())/3.);
	    screen = screenFactor / (epsilon * (1. - epsilon));
	    gReject = (ScreenFunction1(screen) - fZ) / f10 ;
	  }
	else
	  {
	    epsilon = epsilonMin + epsilonRange * rndmEngine->flat();
	    screen = screenFactor / (epsilon * (1 - epsilon));
	    gReject = (ScreenFunction2(screen) - fZ) / f20 ;
	  }
      } while ( gReject < rndmEngine->flat() );
    
  }   //  End of epsilon sampling

  // Fix charges randomly

  G4double electronTotEnergy;
  G4double positronTotEnergy;

  if (rndmEngine->flat() > 0.5)
    {
      electronTotEnergy = (1. - epsilon) * photonEnergy;
      positronTotEnergy = epsilon * photonEnergy;
    }
  else
    {
      positronTotEnergy = (1. - epsilon) * photonEnergy;
      electronTotEnergy = epsilon * photonEnergy;
    }

  // Scattered electron (positron) angles. ( Z - axis along the parent photon)
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211),
  // derived from Tsai distribution (Rev. Mod. Phys. 49, 421 (1977)

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
  
  G4double electronKineEnergy = std::max(0.,electronTotEnergy - electron_mass_c2) ;
  
  G4ThreeVector electronDirection (sinte*cosp, sinte*sinp, coste);
  electronDirection.rotateUz(photonDirection);
      
  G4DynamicParticle* particle1 = new G4DynamicParticle (G4Electron::Electron(),
							electronDirection,
							electronKineEnergy);

  // The e+ is always created 
  G4double positronKineEnergy = std::max(0.,positronTotEnergy - electron_mass_c2) ;

  G4ThreeVector positronDirection (-sintp*cosp, -sintp*sinp, costp);
  positronDirection.rotateUz(photonDirection);   
  
  // Create G4DynamicParticle object for the particle2 
  G4DynamicParticle* particle2 = new G4DynamicParticle(G4Positron::Positron(),
						       positronDirection, 
						       positronKineEnergy);
  // Fill output vector
  fvect->push_back(particle1);
  fvect->push_back(particle2);

  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);   

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AutoLock.hh"
namespace { G4Mutex LivermoreGammaConversionModelMutex = G4MUTEX_INITIALIZER; }

void G4LivermoreGammaConversionModel::InitialiseForElement(
				      const G4ParticleDefinition* part, 
				      G4int Z)
{
  if(GetTripletModel()) { GetTripletModel()->InitialiseForElement(part, Z); }
  G4AutoLock l(&LivermoreGammaConversionModelMutex);
  //  G4cout << "G4LivermoreGammaConversionModel::InitialiseForElement Z= " 
  //   << Z << G4endl;
  if(!data[Z]) { ReadData(Z); }
  if(GetTripletModel() && !probTriplet[Z]) { InitialiseProbability(part, Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LivermoreGammaConversionModel::InitialiseProbability(
                     const G4ParticleDefinition* part, G4int Z)
{
  if(!probTriplet[Z]) {
    const G4Material* mat = (CurrentCouple()) ? CurrentCouple()->GetMaterial()
      : nullptr;
    if(0 == nbinsTriplet) {
      tripletLowEnergy = GetTripletModel()->MinPrimaryEnergy(mat, part, 0.0);
      tripletHighEnergy = 
	std::max(GetTripletModel()->HighEnergyLimit(), 10*tripletLowEnergy);
      G4int nbins = G4EmParameters::Instance()->NumberOfBinsPerDecade();
      nbinsTriplet = std::max(3, 
        (G4int)(nbins*G4Log(tripletHighEnergy/tripletLowEnergy)/(6*G4Log(10.))));
    }
    /*
    G4cout << "G4LivermoreGammaConversionModel::InitialiseProbability Z= " 
	   << Z << " Nbin= " << nbinsTriplet 
	   << "  Emin(MeV)= " << tripletLowEnergy 
	   << "  Emax(MeV)= " << tripletHighEnergy <<  G4endl;
    */
    probTriplet[Z] = 
      new G4PhysicsLogVector(tripletLowEnergy,tripletHighEnergy,nbinsTriplet);
    probTriplet[Z]->SetSpline(true);
    G4double zz = (G4double)Z;
    // loop over bins
    for(G4int j=0; j<=nbinsTriplet; ++j) {
      G4double e = (probTriplet[Z])->Energy(j);
      SetupForMaterial(part, mat, e);
      G4double cross = ComputeCrossSectionPerAtom(part, e, zz);
      G4double tcross = 
	GetTripletModel()->ComputeCrossSectionPerAtom(part, e, zz);
      tcross = (0.0 < cross) ? tcross/cross : 0.0;
      (probTriplet[Z])->PutValue(j, tcross);
      //G4cout << j << ".   E= " << e << " prob= " << tcross << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

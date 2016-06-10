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
// $Id: G4LivermorePolarizedRayleighModel.cc 93810 2015-11-02 11:27:56Z gcosmo $
//
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyPolarizedRayleigh developed by R. Capra
//
// History:
// --------
// 02 May 2009   S Incerti as V. Ivanchenko proposed in G4LivermoreRayleighModel.cc
//
// Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - remove GetMeanFreePath method and table
//                  - remove initialisation of element selector 
//                  - use G4ElementSelector

#include "G4LivermorePolarizedRayleighModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogLogInterpolation.hh"
#include "G4CompositeEMDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4LivermorePolarizedRayleighModel::maxZ = 100;
G4LPhysicsFreeVector* G4LivermorePolarizedRayleighModel::dataCS[] = {0};
G4VEMDataSet* G4LivermorePolarizedRayleighModel::formFactorData = 0;

G4LivermorePolarizedRayleighModel::G4LivermorePolarizedRayleighModel(const G4ParticleDefinition*,
									 const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),isInitialised(false)
{
  fParticleChange =0;
  lowEnergyLimit = 250 * eV; 
  //SetLowEnergyLimit(lowEnergyLimit);
  //SetHighEnergyLimit(highEnergyLimit);
  //
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(verboseLevel > 0) {
    G4cout << "Livermore Polarized Rayleigh is constructed " << G4endl
         << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
         << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedRayleighModel::~G4LivermorePolarizedRayleighModel()
{  
 if(IsMaster()) {
   for(G4int i=0; i<maxZ; ++i) {
     if(dataCS[i]) { 
       delete dataCS[i];
       dataCS[i] = 0;
     }
   }
   delete formFactorData;
   formFactorData = 0; 
   
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedRayleighModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts)
{
// Rayleigh process:                      The Quantum Theory of Radiation
//                                        W. Heitler,       Oxford at the Clarendon Press, Oxford (1954)                                                 
// Scattering function:                   A simple model of photon transport
//                                        D.E. Cullen,      Nucl. Instr. Meth. in Phys. Res. B 101 (1995) 499-510                                       
// Polarization of the outcoming photon:  Beam test of a prototype detector array for the PoGO astronomical hard X-ray/soft gamma-ray polarimeter
//                                        T. Mizuno et al., Nucl. Instr. Meth. in Phys. Res. A 540 (2005) 158-168                                        

  if (verboseLevel > 3)
    G4cout << "Calling G4LivermorePolarizedRayleighModel::Initialise()" << G4endl;


  if(IsMaster()) {
    
    // Form Factor 
    
    G4VDataSetAlgorithm* ffInterpolation = new G4LogLogInterpolation;
    G4String formFactorFile = "rayl/re-ff-";
    formFactorData = new G4CompositeEMDataSet(ffInterpolation,1.,1.);
    formFactorData->LoadData(formFactorFile);
    
    // Initialise element selector
    InitialiseElementSelectors(particle, cuts);
    
    // Access to elements
    char* path = getenv("G4LEDATA");
    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = theCoupleTable->GetTableSize();
   
     for(G4int i=0; i<numOfCouples; ++i) 
       {
	 const G4MaterialCutsCouple* couple = 
	   theCoupleTable->GetMaterialCutsCouple(i);
	 const G4Material* material = couple->GetMaterial();
	 const G4ElementVector* theElementVector = material->GetElementVector();
	 G4int nelm = material->GetNumberOfElements();
	 
	 for (G4int j=0; j<nelm; ++j) 
	   {
	     G4int Z = G4lrint((*theElementVector)[j]->GetZ());
	     if(Z < 1)          { Z = 1; }
	     else if(Z > maxZ)  { Z = maxZ; }
	     if( (!dataCS[Z]) ) { ReadData(Z, path); }
	   }
       }
  }
  
  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedRayleighModel::InitialiseLocal(const G4ParticleDefinition*,
							  G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedRayleighModel::ReadData(size_t Z, const char* path)
{
  if (verboseLevel > 1) 
    {
      G4cout << "Calling ReadData() of G4LivermoreRayleighModel" 
	     << G4endl;
    }
  
  if(dataCS[Z]) { return; }
  
  const char* datadir = path;
  
  if(!datadir) 
    {
      datadir = getenv("G4LEDATA");
      if(!datadir) 
	{
	  G4Exception("G4LivermoreRayleighModelModel::ReadData()","em0006",
		      FatalException,
		      "Environment variable G4LEDATA not defined");
	  return;
	}
    }
  
  //
  
  dataCS[Z] = new G4LPhysicsFreeVector();
  
  // Activation of spline interpolation
  //dataCS[Z] ->SetSpline(true);
  
  std::ostringstream ostCS;
  ostCS << datadir << "/livermore/rayl/re-cs-" << Z <<".dat";
  std::ifstream finCS(ostCS.str().c_str());
  
  if( !finCS .is_open() ) 
    {
     G4ExceptionDescription ed;
     ed << "G4LivermorePolarizedRayleighModel data file <" << ostCS.str().c_str()
        << "> is not opened!" << G4endl;
     G4Exception("G4LivermorePolarizedRayleighModel::ReadData()","em0003",FatalException,
		 ed,"G4LEDATA version should be G4EMLOW6.27 or later.");
     return;
   } 
   else 
   {
     if(verboseLevel > 3) { 
       G4cout << "File " << ostCS.str() 
        << " is opened by G4LivermoreRayleighModel" << G4endl;
     }
     dataCS[Z]->Retrieve(finCS, true);
   } 
 }
 
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
 G4double G4LivermorePolarizedRayleighModel::ComputeCrossSectionPerAtom(
                                        const G4ParticleDefinition*,
					G4double GammaEnergy,
					G4double Z, G4double,
					G4double, G4double)
 {
   if (verboseLevel > 1) 
     {
       G4cout << "G4LivermoreRayleighModel::ComputeCrossSectionPerAtom()" 
	      << G4endl;
     }
 
   if(GammaEnergy < lowEnergyLimit) { return 0.0; }
   
   G4double xs = 0.0;
   
   G4int intZ = G4lrint(Z);
   
   if(intZ < 1 || intZ > maxZ) { return xs; }
   
   G4LPhysicsFreeVector* pv = dataCS[intZ];
 
   // if element was not initialised
   // do initialisation safely for MT mode
   if(!pv) { 
     InitialiseForElement(0, intZ);
     pv = dataCS[intZ];
     if(!pv) { return xs; }
   }
 
   G4int n = pv->GetVectorLength() - 1;
   G4double e = GammaEnergy/MeV;
   if(e >= pv->Energy(n)) {
     xs = (*pv)[n]/(e*e);  
   } else if(e >= pv->Energy(0)) {
     xs = pv->Value(e)/(e*e);  
   }
 
   /*   if(verboseLevel > 0)
	{
	G4cout  <<  "****** DEBUG: tcs value for Z=" << Z << " at energy (MeV)=" 
     << e << G4endl;
       G4cout  <<  "  cs (Geant4 internal unit)=" << xs << G4endl;
     G4cout  <<  "    -> first E*E*cs value in CS data file (iu) =" << (*pv)[0] 
     << G4endl;
       G4cout  <<  "    -> last  E*E*cs value in CS data file (iu) =" << (*pv)[n] 
     << G4endl;
       G4cout  <<  "*********************************************************" 
     << G4endl;
       }
   */
   
   return xs;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedRayleighModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
					      const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4LivermorePolarizedRayleighModel" << G4endl;

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();
  
  if (photonEnergy0 <= lowEnergyLimit)
  {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
  }

  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  // G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy0);
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = (G4int)elm->GetZ();

  G4double outcomingPhotonCosTheta = GenerateCosTheta(photonEnergy0, Z);
  G4double outcomingPhotonPhi = GeneratePhi(outcomingPhotonCosTheta);
  G4double beta=GeneratePolarizationAngle();
 
  // incomingPhoton reference frame:
  // z = versor parallel to the incomingPhotonDirection
  // x = versor parallel to the incomingPhotonPolarization
  // y = defined as z^x
 
  // outgoingPhoton reference frame:
  // z' = versor parallel to the outgoingPhotonDirection
  // x' = defined as x-x*z'z' normalized
  // y' = defined as z'^x'
 
  G4ThreeVector z(aDynamicGamma->GetMomentumDirection().unit()); 
  G4ThreeVector x(GetPhotonPolarization(*aDynamicGamma));
  G4ThreeVector y(z.cross(x));
 
  // z' = std::cos(phi)*std::sin(theta) x + std::sin(phi)*std::sin(theta) y + std::cos(theta) z
  G4double xDir;
  G4double yDir;
  G4double zDir;
  zDir=outcomingPhotonCosTheta;
  xDir=std::sqrt(1-outcomingPhotonCosTheta*outcomingPhotonCosTheta);
  yDir=xDir;
  xDir*=std::cos(outcomingPhotonPhi);
  yDir*=std::sin(outcomingPhotonPhi);
 
  G4ThreeVector zPrime((xDir*x + yDir*y + zDir*z).unit());
  G4ThreeVector xPrime(x.perpPart(zPrime).unit());
  G4ThreeVector yPrime(zPrime.cross(xPrime));
 
  // outgoingPhotonPolarization is directed as x' std::cos(beta) + y' std::sin(beta)
  G4ThreeVector outcomingPhotonPolarization(xPrime*std::cos(beta) + yPrime*std::sin(beta));
 
  fParticleChange->ProposeMomentumDirection(zPrime);
  fParticleChange->ProposePolarization(outcomingPhotonPolarization);
  fParticleChange->SetProposedKineticEnergy(photonEnergy0); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedRayleighModel::GenerateCosTheta(G4double incomingPhotonEnergy, G4int zAtom) const
{
  //  d sigma                                                                    k0
  // --------- =  r0^2 * pi * F^2(x, Z) * ( 2 - sin^2 theta) * std::sin (theta), x = ---- std::sin(theta/2)
  //  d theta                                                                    hc
 
  //  d sigma                                             k0          1 - y
  // --------- = r0^2 * pi * F^2(x, Z) * ( 1 + y^2), x = ---- std::sqrt ( ------- ), y = std::cos(theta)
  //    d y                                               hc            2

  //              Z
  // F(x, Z) ~ --------
  //            a + bx
  //
  // The time to exit from the outer loop grows as ~ k0
  // On pcgeant2 the time is ~ 1 s for k0 ~ 1 MeV on the oxygen element. A 100 GeV
  // event will take ~ 10 hours.
  //
  // On the avarage the inner loop does 1.5 iterations before exiting
 
  const G4double xFactor = (incomingPhotonEnergy*cm)/(h_Planck*c_light);
  //const G4VEMDataSet * formFactorData = GetScatterFunctionData();

  G4double cosTheta;
  G4double fCosTheta;
  G4double x;
  G4double fValue;

  if (incomingPhotonEnergy > 5.*MeV)
  {
    cosTheta = 1.;
  }
  else
  {
    do
    {
      do
	{
	  cosTheta = 2.*G4UniformRand()-1.;
	  fCosTheta = (1.+cosTheta*cosTheta)/2.;
	}
      while (fCosTheta < G4UniformRand());
  
      x = xFactor*std::sqrt((1.-cosTheta)/2.);
  
      if (x > 1.e+005)
	fValue = formFactorData->FindValue(x, zAtom-1);
      else
	fValue = formFactorData->FindValue(0., zAtom-1);
   
      fValue/=zAtom;
      fValue*=fValue;
    }
    while(fValue < G4UniformRand());
  }

  return cosTheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedRayleighModel::GeneratePhi(G4double cosTheta) const
{
  //  d sigma
  // --------- = alpha * ( 1 - sin^2 (theta) * cos^2 (phi) )
  //   d phi
 
  // On the average the loop takes no more than 2 iterations before exiting 

  G4double phi;
  G4double cosPhi;
  G4double phiProbability;
  G4double sin2Theta;
 
  sin2Theta=1.-cosTheta*cosTheta;
 
  do
    {
      phi = twopi * G4UniformRand();
      cosPhi = std::cos(phi);
      phiProbability= 1. - sin2Theta*cosPhi*cosPhi;
    }
  while (phiProbability < G4UniformRand());
 
  return phi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedRayleighModel::GeneratePolarizationAngle(void) const
{
  // Rayleigh polarization is always on the x' direction

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedRayleighModel::GetPhotonPolarization(const G4DynamicParticle&  photon)
{

// SI - From G4VLowEnergyDiscretePhotonProcess.cc
 
  G4ThreeVector photonMomentumDirection;
  G4ThreeVector photonPolarization;

  photonPolarization = photon.GetPolarization(); 
  photonMomentumDirection = photon.GetMomentumDirection();

  if ((!photonPolarization.isOrthogonal(photonMomentumDirection, 1e-6)) || photonPolarization.mag()==0.)
    {
      // if |photonPolarization|==0. or |photonPolarization * photonDirection0| > 1e-6 * |photonPolarization ^ photonDirection0|
      // then polarization is choosen randomly.
  
      G4ThreeVector e1(photonMomentumDirection.orthogonal().unit());
      G4ThreeVector e2(photonMomentumDirection.cross(e1).unit());
  
      G4double angle(G4UniformRand() * twopi);
  
      e1*=std::cos(angle);
      e2*=std::sin(angle);
  
      photonPolarization=e1+e2;
    }
  else if (photonPolarization.howOrthogonal(photonMomentumDirection) != 0.)
    {
      // if |photonPolarization * photonDirection0| != 0.
      // then polarization is made orthonormal;
  
      photonPolarization=photonPolarization.perpPart(photonMomentumDirection);
    }
 
  return photonPolarization.unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
 #include "G4AutoLock.hh"
 namespace { G4Mutex LivermorePolarizedRayleighModelMutex = G4MUTEX_INITIALIZER; }
 
void  G4LivermorePolarizedRayleighModel::InitialiseForElement(const G4ParticleDefinition*, 
                  G4int Z)
 {
   G4AutoLock l(&LivermorePolarizedRayleighModelMutex);
   //  G4cout << "G4LivermoreRayleighModel::InitialiseForElement Z= " 
   //   << Z << G4endl;
   if(!dataCS[Z]) { ReadData(Z); }
   l.unlock();
 }

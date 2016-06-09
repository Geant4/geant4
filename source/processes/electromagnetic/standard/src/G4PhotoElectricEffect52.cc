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
// $Id: G4PhotoElectricEffect52.cc,v 1.2 2006/06/29 19:53:28 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 12-06-96, Added SelectRandomAtom() method, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, PartialSumSigma(i)
//           split of ComputeBindingEnergy, M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 04-06-98, in DoIt, secondary production condition:
//                        range > std::min(threshold,safety)
// 13-08-98, new methods SetBining() PrintInfo()
// 17-11-98, use table of Atomic shells in PostStepDoIt
// 06-01-99, use Sandia crossSection below 50 keV, V.Grichine mma
// 20-05-99, protection against very low energy photons ,L.Urban
// 08-06-99, removed this above protection from the DoIt. mma
// 21-06-00, in DoIt, killing photon: aParticleChange.SetEnergyChange(0.); mma
// 22-06-00, in DoIt, absorbe very low energy photon (back to 20-05-99); mma
// 22-02-01, back to 08-06-99 after correc in SandiaTable (materials-V03-00-05)
// 28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation
// 13-07-01, DoIt: suppression of production cut of the electron (mma)
// 06-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 17-09-01, migration of Materials to pure STL (mma)
// 20-09-01, DoIt: fminimalEnergy of generated electron = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)       
// 10-01-02, moved few function from icc to cc
// 17-04-02, Keep only Sandia crossSections. Remove BuildPhysicsTables.
//           Simplify public interface (mma)
// 29-04-02, Generate theta angle of the photoelectron from Sauter-Gavrila
//           distribution (mma) 
// 15-01-03, photoelectron theta ditribution : return costeta=1 if gamma>5
//           (helmut burkhardt)
// 04-05-05, Add 52 to class name (V.Ivanchenko)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhotoElectricEffect52.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4PhotoElectricEffect52::G4PhotoElectricEffect52(const G4String& processName,
    G4ProcessType type):G4VDiscreteProcess (processName, type),
    fminimalEnergy(1*eV)
    
{ PrintInfoDefinition();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
// destructor
 
G4PhotoElectricEffect52::~G4PhotoElectricEffect52()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4bool G4PhotoElectricEffect52::IsApplicable(const G4ParticleDefinition&
                                                        particle)
{
   return ( &particle == G4Gamma::Gamma() ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4PhotoElectricEffect52::ComputeCrossSectionPerAtom(G4double GammaEnergy,
                                                          G4double AtomicNumber)

// returns the photoElectric cross Section in GEANT4 internal units
{
 G4double* SandiaCof 
           = G4SandiaTable::GetSandiaCofPerAtom((int)AtomicNumber,GammaEnergy);
				    
 G4double energy2 = GammaEnergy*GammaEnergy, energy3 = GammaEnergy*energy2, 
          energy4 = energy2*energy2;

 return SandiaCof[0]/GammaEnergy + SandiaCof[1]/energy2 +
        SandiaCof[2]/energy3     + SandiaCof[3]/energy4; 
          
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4PhotoElectricEffect52::ComputeMeanFreePath(G4double GammaEnergy,
                                                    G4Material* aMaterial)

// returns the gamma mean free path in GEANT4 internal units
{
 G4double* SandiaCof = aMaterial->GetSandiaTable()
                                ->GetSandiaCofForMaterial(GammaEnergy);
				    
 G4double energy2 = GammaEnergy*GammaEnergy, energy3 = GammaEnergy*energy2, 
          energy4 = energy2*energy2;

 
 G4double SIGMA = SandiaCof[0]/GammaEnergy + SandiaCof[1]/energy2 +
                  SandiaCof[2]/energy3     + SandiaCof[3]/energy4; 
          
 return SIGMA > DBL_MIN ? 1./SIGMA : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4PhotoElectricEffect52::GetMeanFreePath(const G4Track& aTrack,
                                                       G4double,
                                                       G4ForceCondition*)

// returns the gamma mean free path in GEANT4 internal units
{
 G4double  GammaEnergy = aTrack.GetDynamicParticle()->GetKineticEnergy();
 G4double* SandiaCof   = aTrack.GetMaterial()->GetSandiaTable()
                                   ->GetSandiaCofForMaterial(GammaEnergy);
				    
 G4double energy2 = GammaEnergy*GammaEnergy, energy3 = GammaEnergy*energy2, 
          energy4 = energy2*energy2;

 
 G4double SIGMA = SandiaCof[0]/GammaEnergy + SandiaCof[1]/energy2 +
                  SandiaCof[2]/energy3     + SandiaCof[3]/energy4; 
          
 MeanFreePath = SIGMA > DBL_MIN ? 1./SIGMA : DBL_MAX;
 return MeanFreePath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VParticleChange* G4PhotoElectricEffect52::PostStepDoIt(const G4Track& aTrack,
                                                      const G4Step&  aStep)
//
// Generate an electron resulting of a photo electric effect.
// The incident photon disappear.
// GEANT4 internal units
//
 
{  aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();

   const G4DynamicParticle* aDynamicPhoton = aTrack.GetDynamicParticle();

   G4double PhotonEnergy = aDynamicPhoton->GetKineticEnergy();
   G4ParticleMomentum PhotonDirection = aDynamicPhoton->GetMomentumDirection();
   
   // select randomly one element constituing the material.
   G4Element* anElement = SelectRandomAtom(aDynamicPhoton, aMaterial);

   //
   // Photo electron
   //

   G4int NbOfShells = anElement->GetNbOfAtomicShells();
   G4int i=0;
   while ((i<NbOfShells)&&(PhotonEnergy<anElement->GetAtomicShell(i))) i++;

   if (i==NbOfShells) return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
      
   G4double ElecKineEnergy = PhotonEnergy - anElement->GetAtomicShell(i);

   if (ElecKineEnergy > fminimalEnergy)
     {
      // direction of the photo electron
      //
      G4double cosTeta = ElecThetaDistribution(ElecKineEnergy);
      G4double sinTeta = sqrt(1.-cosTeta*cosTeta);
      G4double Phi     = twopi * G4UniformRand();
      G4double dirx = sinTeta*cos(Phi),diry = sinTeta*sin(Phi),dirz = cosTeta;
      G4ThreeVector ElecDirection(dirx,diry,dirz);
      ElecDirection.rotateUz(PhotonDirection);
      // 
      G4DynamicParticle* aElectron = new G4DynamicParticle (
                        G4Electron::Electron(),ElecDirection, ElecKineEnergy);
      aParticleChange.SetNumberOfSecondaries(1);
      aParticleChange.AddSecondary(aElectron); 
     }
   else
     {
      ElecKineEnergy = 0.;
      aParticleChange.SetNumberOfSecondaries(0);
     }

   //
   // Kill the incident photon 
   //
   aParticleChange.ProposeLocalEnergyDeposit(PhotonEnergy-ElecKineEnergy);
   aParticleChange.ProposeEnergy(0.);  
   aParticleChange.ProposeTrackStatus(fStopAndKill); 

   //  Reset NbOfInteractionLengthLeft and return aParticleChange
   return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Element* G4PhotoElectricEffect52::SelectRandomAtom(
                                     const G4DynamicParticle* aDynamicPhoton,
                                           G4Material* aMaterial)
{
  // select randomly 1 element within the material

  const G4int NumberOfElements            = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  if (NumberOfElements == 1) return (*theElementVector)[0];

  G4double GammaEnergy = aDynamicPhoton->GetKineticEnergy();
  const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

  G4double PartialSumSigma = 0. ;
  G4double rval = G4UniformRand();
 
  for ( G4int elm=0 ; elm < NumberOfElements ; elm++ )
     {PartialSumSigma += NbOfAtomsPerVolume[elm] *
                         ComputeCrossSectionPerAtom(GammaEnergy,
                                         (*theElementVector)[elm]->GetZ());
      if (rval<=PartialSumSigma*MeanFreePath) return ((*theElementVector)[elm]);
     }
  return ((*theElementVector)[NumberOfElements-1]);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4PhotoElectricEffect52::ElecThetaDistribution(G4double kineEnergy)
{
 // Compute Theta distribution of the emitted electron, with respect to the
 // incident Gamma.
 // The Sauter-Gavrila distribution for the K-shell is used.
 //
 G4double costeta = 1.;
 G4double gamma   = 1. + kineEnergy/electron_mass_c2;
 if (gamma > 5.) return costeta;
 G4double beta  = sqrt(gamma*gamma-1.)/gamma;
 G4double b     = 0.5*gamma*(gamma-1.)*(gamma-2);
    
 G4double rndm,term,greject,grejsup;
 if (gamma < 2.) grejsup = gamma*gamma*(1.+b-beta*b);
 else            grejsup = gamma*gamma*(1.+b+beta*b);
  
 do { rndm = 1.-2*G4UniformRand();
      costeta = (rndm+beta)/(rndm*beta+1.);
      term = 1.-beta*costeta;
      greject = (1.-costeta*costeta)*(1.+b*term)/(term*term);
 } while(greject < G4UniformRand()*grejsup);
       
 return costeta;      
     
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PhotoElectricEffect52::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from Sandia parametrisation. ";
	             
  G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

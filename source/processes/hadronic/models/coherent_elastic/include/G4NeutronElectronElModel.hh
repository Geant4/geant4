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
// $Id: G4NeutronElectronElModel.hh 90228 2015-05-21 08:49:57Z gcosmo $
//
// Geant4 Header : G4NeutronElectronElModel
//
// 16.5.17 V.Grichine 
//  
// Modified:
//
// Class Description
// Default model for neutron-electron elastic  scattering; 
// Class Description - End

#ifndef G4NeutronElectronElModel_h
#define G4NeutronElectronElModel_h 1
 
#include "globals.hh"
#include "G4HadronElastic.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"

using namespace std;
using namespace CLHEP;

class G4ParticleDefinition;
class G4PhysicsLogVector;
class G4PhysicsTable;

class G4NeutronElectronElModel : public G4HadronElastic
{
public:

  G4NeutronElectronElModel(const G4String& name = "nu-e-elastic");

  virtual ~G4NeutronElectronElModel();

  virtual G4bool IsApplicable(const G4HadProjectile & aTrack, 
  			      G4Nucleus & targetNucleus);

  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
					  G4Nucleus & targetNucleus);

  void Initialise();

  G4double XscIntegrand(G4double x);
  
  G4double CalculateAm( G4double momentum);

  inline G4double GetAm(){return fAm;};


  // sample recoil electron energy in lab frame

  G4double SampleSin2HalfTheta(G4double Tkin);
  G4double GetTransfer( G4int iTkin, G4int iTransfer, G4double position );

  void SetCutEnergy(G4double ec){fCutEnergy=ec;};
  G4double GetCutEnergy(){return fCutEnergy;};


  
  virtual void ModelDescription(std::ostream&) const;

private:

  G4double fAm, fM, fM2, fMv2, fme, fme2, fee, fee2;
  G4double fMinEnergy, fMaxEnergy;
  G4int fEnergyBin, fAngleBin;
  G4ParticleDefinition* theElectron; 
  G4double fCutEnergy; // minimal recoil electron energy detected
  G4PhysicsLogVector* fEnergyVector;
  G4PhysicsTable* fAngleTable;
};

////////////////////////////////////////////////////////////////////
//
// return Wentzel atom screening correction for neutron-electron scattering

inline  G4double G4NeutronElectronElModel::CalculateAm( G4double Tkin)
{
  fee      = (Tkin+fM)*fme/fM;
    // G4cout<<"fee = "<<fee<<" MeV"<<G4endl;
  fee2     = fee*fee;
  G4double   momentum = sqrt( fee2 - fme2 );
  G4double k   = momentum/CLHEP::hbarc;
  G4double ch  = 1.13;
  G4double zn  = 1.77*k*CLHEP::Bohr_radius;
  G4double zn2 = zn*zn;
  fAm          = ch/zn2;

  return fAm;
}


#endif

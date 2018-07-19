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
// Neutron-electron elastic cross section base on the integration of
// the Rosenbluth differential xsc
//
// 16.05.17 V. Grichine
//
//

#ifndef G4NeutronElectronElXsc_h
#define G4NeutronElectronElXsc_h


#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"

using namespace std;
using namespace CLHEP;

// class G4ParticleDefinition;
class G4PhysicsLogVector;
class G4PhysicsTable;

class G4NeutronElectronElXsc : public G4VCrossSectionDataSet
{
public:
   
  G4NeutronElectronElXsc();
  ~G4NeutronElectronElXsc();

  void Initialise();

  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z, const G4Material*);


  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

  G4double GetRosenbluthXsc(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

  G4double XscIntegrand(G4double x);

   G4double GetElementNonRelXsc(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

  G4double CalculateAm( G4double momentum);

  inline G4double GetAm(){return fAm;};

  void SetCutEnergy(G4double ec){fCutEnergy=ec;};
  G4double GetCutEnergy(){return fCutEnergy;};

  void SetBiasingFactor(G4double bf){fBiasingFactor=bf;};

protected:
  G4double fM, fM2, fMv2, fme, fme2, fee, fee2;
  G4double fCofXsc;    // 
  G4double fAm;    //
  G4int fEnergyBin; 
  G4double fMinEnergy, fMaxEnergy, fCutEnergy; // minimal recoil electron energy detected
  G4double fBiasingFactor; // biasing xsc up

  G4PhysicsLogVector* fEnergyXscVector;
  static const G4double fXscArray[200];
};



////////////////////////////////////////////////////////////////////
//
// return Wentzel atom screening correction for neutron-electron scattering

inline  G4double G4NeutronElectronElXsc::CalculateAm( G4double momentum)
{
  G4double k   = momentum/CLHEP::hbarc;
  G4double ch  = 1.13;
  G4double zn  = 1.77*k*CLHEP::Bohr_radius;
  G4double zn2 = zn*zn;
  fAm          = ch/zn2;

  return fAm;
}

////////////////////////////////////////////////////
//
// Slow electron (Tkin << me_c2) in the neutron rest frame

inline G4double G4NeutronElectronElXsc::
GetElementNonRelXsc(const G4DynamicParticle* aPart, G4int ZZ,  
		       const G4Material*) 
{
  G4double result(0.), te(0.), momentum(0.);

  te = aPart->GetKineticEnergy()*fme/fM;
  momentum = sqrt( te*(te + 2.*fme) );
  fAm = CalculateAm(momentum);

  result = 1. + log(1. +1./fAm);
  result *= fCofXsc; //*energy;
  result *= ZZ;  // incoherent sum over  all element electrons

  return result;
}

#endif

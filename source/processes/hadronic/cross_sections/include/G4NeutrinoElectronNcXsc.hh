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
// Neutrino electron elastic (neutral current) cross sections
//
// 02.04.17 V. Grichine
//
//

#ifndef G4NeutrinoElectronNcXsc_h
#define G4NeutrinoElectronNcXsc_h


#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class G4NeutrinoElectronNcXsc : public G4VCrossSectionDataSet
{
public:
   
  G4NeutrinoElectronNcXsc();
  ~G4NeutrinoElectronNcXsc();

  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z, const G4Material*);


  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

  void SetCutEnergy(G4double ec){fCutEnergy=ec;};
  G4double GetCutEnergy(){return fCutEnergy;};

  void SetBiasingFactor(G4double bf){fBiasingFactor=bf;};

protected:

  G4double fCofXsc;    // 2*Gf*Gf*MeC2/pi
  G4double fSin2tW;    // sin^2theta_Weinberg
  G4double fCutEnergy; // minimal recoil electron energy detected
  G4double fBiasingFactor; // biasing xsc up
};

#endif

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
// Neutrino electron total (Cc+Nc) cross sections
//
// 27.11.17 V. Grichine
//
//

#ifndef G4NeutrinoElectronTotXsc_h
#define G4NeutrinoElectronTotXsc_h


#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class G4NeutrinoElectronCcXsc;
class G4NeutrinoElectronNcXsc;

class G4NeutrinoElectronTotXsc : public G4VCrossSectionDataSet
{
public:
   
  G4NeutrinoElectronTotXsc();
  ~G4NeutrinoElectronTotXsc();

  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z, const G4Material*);


  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

  void SetCutEnergy(G4double ec){fCutEnergy=ec;};
  G4double GetCutEnergy(){return fCutEnergy;};

  void SetBiasingFactor(G4double bf);
  void SetBiasingFactors(G4double bfCc, G4double bfNc); // for separate testing
  G4double GetBiasingFactor(){return fBiasingFactor;};
  G4double GetCcRatio(){return fCcRatio;};

protected:

  G4NeutrinoElectronCcXsc* fCcXsc;
  G4NeutrinoElectronNcXsc* fNcXsc;
  G4double fCutEnergy; // min detected recoil energy
  G4double fBiasingFactor; // biasing xsc up
  G4double fCcRatio; // biasing xsc up

};

#endif

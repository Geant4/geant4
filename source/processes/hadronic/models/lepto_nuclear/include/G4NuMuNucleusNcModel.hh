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
// $Id: G4NuMuNucleusNcModel.hh 90228 2015-05-21 08:49:57Z gcosmo $
//
// Geant4 Header : G4NuMuNucleusNcModel
//
// Author : V.Grichine 27.2.19
//  
// Modified:
//
// Class Description
// Default model for muon neutrino-nucleus neutral current scattering; 
// Class Description - End

#ifndef G4NuMuNucleusNcModel_h
#define G4NuMuNucleusNcModel_h 1
 
#include "globals.hh"
#include "G4NeutrinoNucleusModel.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"
#include "G4LorentzVector.hh"
#include "G4Threading.hh"

class G4ParticleDefinition;
class G4VPreCompoundModel;
class G4CascadeInterface;
class G4BinaryCascade;
class G4TheoFSGenerator;
class G4LundStringFragmentation;
class G4ExcitedStringDecay;
class G4INCLXXInterface;
class G4Nucleus;

class G4NuMuNucleusNcModel : public G4NeutrinoNucleusModel
{
public:

  G4NuMuNucleusNcModel(const G4String& name = "NuMuNuclNcModel");

  virtual ~G4NuMuNucleusNcModel();

  virtual void InitialiseModel();

  virtual G4bool IsApplicable(const G4HadProjectile & aTrack, 
  			      G4Nucleus & targetNucleus);

  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
					  G4Nucleus & targetNucleus);


  ////////// KR excitation kinematics ////////////////////

  void SampleLVkr(const G4HadProjectile & aTrack, G4Nucleus & targetNucleus);

  G4double GetMinNuMuEnergy(){ return fMnumu + 0.5*fMnumu*fMnumu/fM1 + 4.*CLHEP::keV; }; // kinematics + accuracy for sqrts

  G4double ThresholdEnergy(G4double mI, G4double mF, G4double mP) // for cluster decay
  { 
    G4double w = std::sqrt(fW2);
    return w + 0.5*( (mP+mF)*(mP+mF)-(w+mI)*(w+mI) )/mI;
  };
 
  virtual void ModelDescription(std::ostream&) const;

private:

  G4ParticleDefinition* theNuMu;
  G4ParticleDefinition* theANuMu;

  G4double  fMnumu; // = 0 for <f|-state
 
  G4bool fData, fMaster; // for one initialisation only

#ifdef G4MULTITHREADED
  static G4Mutex numuNucleusModel;
#endif

};



#endif

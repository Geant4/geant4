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
// $Id: G4DiscreteGammaTransition.hh 85841 2014-11-05 15:35:06Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4DiscreteGammaTransition
//
//      Author:        Maria Grazia Pia   (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//
//        8 March 2002,  Fan Lei (flei@space.qinetiq.com)
//              added 
//              1) SetRDM () to switch on/off IC process
//              2) GetOrbitNumber () to return the CE oribit number
//              3) GetBondEnergy() for the converted e-
//              4) IsaGamma() to separate CE from gamma
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              i) added G4int _nucleusZ initialise it through the constructor
//              ii) modified SelectGamma() to allow the generation of conversion electrons
//              iii) added #include G4AtomicShells.hh
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//
//      
//        19 April 2010, J. M. Quesada. 
//              Corrections added for taking into account mismatch between tabulated 
//              gamma energies and level energy differences (fake photons eliminated) 
//   
//	  6 October 2010, M. Kelsey
//		Store NuclearLevel as const-reference, not as copied value.
//
//	  5 May 2011, V.Ivanchenko removed unused constructor
//
// -------------------------------------------------------------------

#ifndef G4DiscreteGammaTransition_hh
#define G4DiscreteGammaTransition_hh

#include "globals.hh"
#include "G4VGammaTransition.hh"

class G4NuclearLevel;
class G4NuclearLevelManager;


class G4DiscreteGammaTransition : public G4VGammaTransition
{
public:

  G4DiscreteGammaTransition(const G4NuclearLevel* level, 
			    G4int Z, G4int ver);

  virtual ~G4DiscreteGammaTransition();

  virtual void SetEnergyFrom(G4double energy);
  virtual G4double GetGammaEnergy();
  virtual G4double GetGammaCreationTime();
  virtual void SelectGamma();

  inline void SetICM(G4bool ic)    { icm = ic; };
  inline G4bool GetICM() const     { return icm;};
  inline G4double GetBondEnergy () { return bondE;};
  inline G4int GetOrbitNumber ()   { return orbitE;};
  inline G4bool IsAGamma()         { return aGamma;};

  inline void Update(const G4NuclearLevel* level, G4int Z)
  {
    aLevel = level; 
    nucleusZ = Z; 
  }

private:
  
  G4int nucleusZ;
  G4int orbitE;
  G4int verbose;

  G4double bondE;
  G4double gammaEnergy;
  G4double excitation;
  G4double gammaCreationTime;

  G4bool aGamma;
  G4bool icm;

  const G4NuclearLevel* aLevel;
};

#endif


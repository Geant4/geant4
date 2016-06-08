//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4ContinuumGammaTransition
//
//      Authors:       Carlo Dallapiccola (dallapiccola@umdhep.umd.edu)
//                     Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//      
// -------------------------------------------------------------------
//
//  Header file for G4ContinuumGammaTransition
//

#ifndef G4ContinuumGammaTransition_hh
#define G4ContinuumGammaTransition_hh

#include "globals.hh"
#include "G4VGammaTransition.hh"
#include "G4NuclearLevelManager.hh"
#include "G4VLevelDensityParameter.hh"

class G4ContinuumGammaTransition : public G4VGammaTransition
{
public:

  // Constructor
  G4ContinuumGammaTransition(const G4NuclearLevelManager& levelManager,
			     G4int Z, G4int A, G4double excitation,
			     G4int verbose);

  // Destructor
  ~G4ContinuumGammaTransition();

  // Functions

//--  virtual G4double GammaEnergy();
//--  virtual G4double GetEnergyTo() const;
  virtual void SetEnergyFrom(const G4double energy);
  virtual G4double GetGammaEnergy();
  virtual G4double GetGammaCreationTime();
  virtual void SelectGamma();

private:

  G4double E1Pdf(G4double energy);
  G4double GammaTime();

  G4int _nucleusA;
  G4int _nucleusZ;
  G4double _eMin;
  G4double _eMax;
  G4double _maxLevelE;
  G4double _minLevelE;
  G4double _excitation;
  G4double _eGamma;
  G4NuclearLevelManager _levelManager;
  G4double _gammaCreationTime;

};

#endif

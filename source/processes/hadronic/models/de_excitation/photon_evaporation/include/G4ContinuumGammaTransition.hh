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
// $Id$
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
//        02 May 2003,   Vladimir Ivanchenko change interface to G4NuclearlevelManager
//      
// -------------------------------------------------------------------
//
//  Header file for G4ContinuumGammaTransition
//

#ifndef G4ContinuumGammaTransition_hh
#define G4ContinuumGammaTransition_hh 1

#include "globals.hh"
#include "G4VGammaTransition.hh"
#include "G4NuclearLevelManager.hh"
#include "G4VLevelDensityParameter.hh"

class G4ContinuumGammaTransition : public G4VGammaTransition
{
public:

  // Constructor
  G4ContinuumGammaTransition(const G4NuclearLevelManager* levelManager,
			     G4int Z, G4int A, G4double excitation,
			     G4int verbose);

  // Destructor
  virtual ~G4ContinuumGammaTransition();

  // Functions

  //--  virtual G4double GammaEnergy();
  //--  virtual G4double GetEnergyTo() const;
  virtual void SetEnergyFrom(G4double energy);
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
  const G4NuclearLevelManager* _levelManager;
  G4double _gammaCreationTime;

};

#endif

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
// $Id: G4ContinuumGammaDeexcitation.hh 85841 2014-11-05 15:35:06Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4ContinuumGammaDeexcitation
//
//      Authors:       Carlo Dallapiccola (dallapiccola@umdhep.umd.edu)
//                     Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------
//
//
// Class G4ContinuumGammaDeexcitation.hh
//
#ifndef G4ContinuumGammaDeexcitation_hh
#define G4ContinuumGammaDeexcitation_hh 1

#include "G4VGammaDeexcitation.hh"

#include "globals.hh"

#include "G4ContinuumGammaTransition.hh"
#include "G4Fragment.hh"

class G4NuclearLevelManager;
class G4NuclearLevelStore;
class G4DiscreteGammaTransition;

class G4ContinuumGammaDeexcitation : public G4VGammaDeexcitation
{
public:

  G4ContinuumGammaDeexcitation();

  virtual ~G4ContinuumGammaDeexcitation();

  virtual G4bool CanDoTransition(G4Fragment* aNucleus);

private:

  G4int nucleusZ;
  G4int nucleusA;  
  G4NuclearLevelStore*   store;
  G4NuclearLevelManager* levelManager;
  G4ContinuumGammaTransition* ctransition; 
};

#endif

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
// $Id: G4ContinuumGammaDeexcitation.hh,v 1.5 2010-11-17 19:17:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

class G4ContinuumGammaDeexcitation : public G4VGammaDeexcitation
{
public:

  // Constructor
  G4ContinuumGammaDeexcitation();

  // Destructor
  virtual ~G4ContinuumGammaDeexcitation();

  // Functions

public:

  virtual G4VGammaTransition* CreateTransition();

  virtual G4bool CanDoTransition();

private:

  G4int _nucleusZ;
  G4int _nucleusA;  
  G4NuclearLevelManager * _levelManager;
};

#endif

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
// $Id: G4PolarizedGammaTransition.hh 85659 2014-11-03 10:59:10Z vnivanch $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PolarizedGammaTransition
//
//      Author:        Vladimir Ivanchenko
// 
//      Creation date: 6 November 2015
//
//      Modifications: 
//
//      
// -------------------------------------------------------------------

#ifndef G4POLARIZEDGAMMATRANSITION_HH
#define G4POLARIZEDGAMMATRANSITION_HH 1

#include "G4GammaTransition.hh"

class G4PolarizationTransition;

class G4PolarizedGammaTransition : public G4GammaTransition
{
public:

  explicit G4PolarizedGammaTransition();

  virtual ~G4PolarizedGammaTransition();
  
  virtual void SampleDirection(G4Fragment* nuc, G4double ratio,
			       G4int twoJ1, G4int twoJ2, G4int mp);

private:  

  G4PolarizedGammaTransition(const G4PolarizedGammaTransition &right) = delete;
  const G4PolarizedGammaTransition& operator=(const G4PolarizedGammaTransition &right) = delete;
  G4bool operator==(const G4PolarizedGammaTransition &right) const = delete;
  G4bool operator!=(const G4PolarizedGammaTransition &right) const = delete;

  G4PolarizationTransition* fPolarization; 
};


#endif






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
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4GammaTransition
//
//      Author:        Vladimir Ivanchenko
// 
//      Creation date: 27 February 2015
//
//      Modifications: 
//
//      
// -------------------------------------------------------------------

#ifndef G4GAMMATRANSITION_HH
#define G4GAMMATRANSITION_HH 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4GammaTransition 
{
public:

  G4GammaTransition();

  virtual ~G4GammaTransition();
  
  virtual G4Fragment* SampleTransition(G4Fragment* nucleus,
				       G4double newExcEnergy,
                                       G4int  deltaS,
                                       size_t shell,
                                       G4bool isGamma,
				       G4bool isLongLived);

private:  

  G4GammaTransition(const G4GammaTransition &right);
  
  const G4GammaTransition& operator=(const G4GammaTransition &right);
  G4bool operator==(const G4GammaTransition &right) const;
  G4bool operator!=(const G4GammaTransition &right) const;
 
};


#endif






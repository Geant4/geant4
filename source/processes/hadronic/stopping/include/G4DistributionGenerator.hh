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
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4DistributionGenerator.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 12 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#ifndef G4DISTRIBUTIONGENERATOR_HH
#define G4DISTRIBUTIONGENERATOR_HH

#include <vector> 

#include "globals.hh"
#include "Randomize.hh"

class G4DistributionGenerator
{  

public:

  // Constructor
  G4DistributionGenerator(std::vector<G4double>& x,
  			  std::vector<G4double>& values);
  G4DistributionGenerator();

  // Destructor
  ~G4DistributionGenerator();

  G4double Generate(G4double ranflat);

private:

  // Hide assignment operator as private 
  //  G4DistributionGenerator& operator=(const G4DistributionGenerator &right);

  // Copy constructor
  //  G4DistributionGenerator(const G4DistributionGenerator& );

  std::vector<G4double> _x;
  std::vector<G4double> _cumProb;

};
 
#endif

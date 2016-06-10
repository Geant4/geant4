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
// $Id: G4HyperNucleiProperties.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
// ------------------------------------------------------------
// Hyper Nuclei properties based on CHIPS model (Mikhail KOSOV)
// Migrate into particles category by H.Kurashige (Sep. 2007)
// 

#ifndef G4HyperNucleiProperties_h
#define G4HyperNucleiProperties_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"

class G4HyperNucleiProperties
{
 // Class Description
 //   G4HyperNucleiProperties is an utility class to provide 
 //   mass formula of hyper nuclei based on CHIPS model
 //   (i.e. it has static member function only)


public: 

  // Destructor
  ~G4HyperNucleiProperties() { };

  // Default constructor ()
  G4HyperNucleiProperties(){};


public:  // With Description

  // Calculate Mass Excess of nucleus A, Z, L(number of Lambda)
  static G4double GetAtomicMass(G4int A, G4int Z, G4int L);
  
  static G4double GetNuclearMass(G4int A, G4int Z, G4int L);

};

#endif

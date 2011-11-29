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
//
// GEANT4 Class file
//
//
// File name:  G4RDModifiedTsai
//
// Author:     Andreia Trindade (andreia@lip.pt)
//             Pedro Rodrigues  (psilva@lip.pt)
//             Luis Peralta     (luis@lip.pt)
// 
// Creation date: 21 March 2003
//
// Modifications: 
// 21 Mar 2003       A. Trindade    First implementation acording with new design
// 24 Mar 2003                      & Fix in Tsai generator in order to prevent theta generation above pi
//               
//
// Class Description: 
//
// Concrete class for Bremsstrahlung Angular Distribution Generation - Tsai Model
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4RDModifiedTsai_h
#define G4RDModifiedTsai_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4RDVBremAngularDistribution.hh"

class G4RDModifiedTsai : public G4RDVBremAngularDistribution
{

public:

  G4RDModifiedTsai(const G4String& name);

  ~G4RDModifiedTsai();

  G4double PolarAngle(const G4double initial_energy,
		      const G4double final_energy,
		      const G4int Z);

  void PrintGeneratorInformation() const;

protected:

private:

  // hide assignment operator 
     G4RDModifiedTsai & operator=(const  G4RDModifiedTsai &right);
     G4RDModifiedTsai(const  G4RDModifiedTsai&);

};

#endif


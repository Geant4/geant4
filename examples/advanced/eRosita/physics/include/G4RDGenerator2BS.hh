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
// File name:  G4RDGenerator2BS
//
// Author:     Andreia Trindade (andreia@lip.pt)
//             Pedro Rodrigues  (psilva@lip.pt)
//             Luis Peralta     (luis@lip.pt)
// 
// Creation date: 2 June 2003
//
// Modifications: 
// 02 Jun 2003                           First implementation acording with new design
//               
//
// Class Description: 
//
// Concrete class for Bremsstrahlung Angular Distribution Generation - 2BS Distribution
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4RDGenerator2BS_h
#define G4RDGenerator2BS_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4RDVBremAngularDistribution.hh"

class G4RDGenerator2BS : public G4RDVBremAngularDistribution
{

public:

  G4RDGenerator2BS(const G4String& name);

  ~G4RDGenerator2BS();

  G4double PolarAngle(const G4double initial_energy,
		      const G4double final_energy,
		      const G4int Z);

  void PrintGeneratorInformation() const;

protected:

  G4double RejectionFunction(G4double value) const;

private:
  G4double z;
  G4double rejection_argument1, rejection_argument2, rejection_argument3;
  G4double EnergyRatio;

  // hide assignment operator 
     G4RDGenerator2BS & operator=(const  G4RDGenerator2BS &right);
     G4RDGenerator2BS(const  G4RDGenerator2BS&);

};

#endif


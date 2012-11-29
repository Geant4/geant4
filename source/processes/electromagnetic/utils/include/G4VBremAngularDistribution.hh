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
//
// GEANT4 Class file
//
//
// File name:  G4VBremAngularDistribution
//
// Author:     Pedro Rodrigues  (psilva@lip.pt)
//             Andreia Trindade (andreia@lip.pt)
//             Luis Peralta     (luis@lip.pt)
//             Design from Maria Grazia Pia (MariaGrazia.Pia@ge.infn.it) 
//            
// 
// Creation date: 21 March 2003
//
// Modifications: 
// 21 Mar 2003  A.Trindade    First implementation acording with new design
// 13 Oct 2010  V.Ivanchenko  Moved to utils  
//
// Class Description: 
//
// Abstract class for Bremsstrahlung Angular Distribution Generation

// -------------------------------------------------------------------
//

#ifndef G4VBremAngularDistribution_h
#define G4VBremAngularDistribution_h 1

#include "globals.hh"
#include "G4VEmAngularDistribution.hh"

class G4VBremAngularDistribution : public G4VEmAngularDistribution
{
public:

  G4VBremAngularDistribution(const G4String& name);

  virtual ~G4VBremAngularDistribution();

  virtual G4double PolarAngle(const G4double initial_energy,
			      const G4double final_energy,
                              const G4int Z) = 0;

  virtual void PrintGeneratorInformation() const;

private:

  // hide assignment operator 
  G4VBremAngularDistribution & operator=(const  G4VBremAngularDistribution &right);
  G4VBremAngularDistribution(const  G4VBremAngularDistribution&);

};

#endif


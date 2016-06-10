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
// $Id: G4IInterpolator.hh 70904 2013-06-07 10:34:25Z gcosmo $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 31 Jul 2008   MGP        Revised and renamed to G4IInterpolator
//
// -------------------------------------------------------------------
// Class description:
// Base class for a strategy pattern dealing with interpolation

// -------------------------------------------------------------------

#ifndef G4IINTERPOLATOR_HH
#define G4IINTERPOLATOR_HH 1

#include "globals.hh"
#include "G4DataVector.hh"

class G4IInterpolator {
 
public:

  G4IInterpolator() { }

  virtual ~G4IInterpolator() { }
 

  virtual G4double Calculate(G4double point, G4int bin, 
			     const G4DataVector& energies, 
			     const G4DataVector& data) const = 0;

  virtual G4IInterpolator* Clone() const = 0;

private:
  
  // Hide copy constructor and assignment operator
  G4IInterpolator(const G4IInterpolator&);
  G4IInterpolator& operator=(const G4IInterpolator& right);

};
 
#endif
 











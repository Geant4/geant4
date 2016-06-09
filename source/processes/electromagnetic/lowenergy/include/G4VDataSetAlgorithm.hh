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
// $Id: G4VDataSetAlgorithm.hh,v 1.7 2006/06/29 19:37:07 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Base class for a strategy pattern to encapsulate algorithms for interpolation of a data set
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4VDATASETALGORITHM_HH
#define G4VDATASETALGORITHM_HH 1

#include "globals.hh"
#include "G4DataVector.hh"

class G4VDataSetAlgorithm {
 
public:

  G4VDataSetAlgorithm() { }

  virtual ~G4VDataSetAlgorithm() { }
 

  virtual G4double Calculate(G4double point, G4int bin, 
			     const G4DataVector& energies, 
			     const G4DataVector& data) const = 0;

  virtual G4VDataSetAlgorithm* Clone() const = 0;

private:
  
  // Hide copy constructor and assignment operator
  G4VDataSetAlgorithm(const G4VDataSetAlgorithm&);
  G4VDataSetAlgorithm& operator=(const G4VDataSetAlgorithm& right);

};
 
#endif
 











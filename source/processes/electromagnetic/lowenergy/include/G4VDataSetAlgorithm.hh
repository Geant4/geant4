//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VDataSetAlgorithm.hh,v 1.4 2001-10-09 11:23:26 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  G4VDataSetAlgorithm & operator=(const G4VDataSetAlgorithm &right);

};
 
#endif
 











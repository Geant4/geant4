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
// $Id: G4RDSemiLogInterpolation.hh 107396 2017-11-10 08:28:08Z gcosmo $
// GEANT4 tag $Name: geant4-09-01-ref-00 $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Log-Log interpolation of a data set
// Part of a strategy pattern to encapsulate algorithms for interpolation of data sets
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4RDSEMILOGINTERPOLATION_HH
#define G4RDSEMILOGINTERPOLATION_HH 1

#include "globals.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4DataVector.hh"

class G4RDSemiLogInterpolation : public G4RDVDataSetAlgorithm {
 
public:

  G4RDSemiLogInterpolation();

  ~G4RDSemiLogInterpolation();
 
  G4double Calculate(G4double point, G4int bin, 
		     const G4DataVector& energies, 
		     const G4DataVector& data) const;

  virtual G4RDVDataSetAlgorithm* Clone() const;

private:

  
  // Hide copy constructor and assignment operator

};
 
#endif
 











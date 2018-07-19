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
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries, 
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class: 
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//

#ifndef G4LOGLOGINTERPOLATION_HH
#define G4LOGLOGINTERPOLATION_HH 1

#include "globals.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4DataVector.hh"

class G4DNACPA100LogLogInterpolation : public G4VDataSetAlgorithm {
 
public:

  G4DNACPA100LogLogInterpolation();

  ~G4DNACPA100LogLogInterpolation();
 
  G4double Calculate(G4double point, G4int bin, 
       const G4DataVector& energies, 
       const G4DataVector& data) const;

  G4double Calculate(G4double point, G4int bin,
       const G4DataVector& energies,
       const G4DataVector& data, 
       const G4DataVector& log_energies, 
       const G4DataVector& log_data) const;

  virtual G4VDataSetAlgorithm* Clone() const; 

private:

  
  // Hide copy constructor and assignment operator

};
 
#endif
 











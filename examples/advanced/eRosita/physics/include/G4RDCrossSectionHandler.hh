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
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
//  1 Aug 2001   MGP        Created
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Data set manager for an electromagnetic physics process
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDCROSSSECTIONHANDLER_HH
#define G4RDCROSSSECTIONHANDLER_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include <map>
#include <vector>
#include "G4RDVCrossSectionHandler.hh"

class G4RDVDataSetAlgorithm;
class G4RDVEMDataSet;
class G4Material;
class G4Element;

class G4RDCrossSectionHandler : public G4RDVCrossSectionHandler {
 
public:

  G4RDCrossSectionHandler();

  ~G4RDCrossSectionHandler();
	
   
protected: 
   
  virtual std::vector<G4RDVEMDataSet*>* BuildCrossSectionsForMaterials(const G4DataVector& energyVector, 
								       const G4DataVector* energyCuts = 0);
 
private:
 
  // Hide copy constructor and assignment operator
  G4RDCrossSectionHandler(const G4RDCrossSectionHandler&);
  G4RDCrossSectionHandler & operator=(const G4RDCrossSectionHandler &right);

};
 
#endif












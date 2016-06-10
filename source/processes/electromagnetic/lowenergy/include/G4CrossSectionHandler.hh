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
// $Id: G4CrossSectionHandler.hh 66241 2012-12-13 18:34:42Z gunter $
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

// -------------------------------------------------------------------

#ifndef G4CROSSSECTIONHANDLER_HH
#define G4CROSSSECTIONHANDLER_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include <map>
#include <vector>
#include "G4VCrossSectionHandler.hh"

class G4VDataSetAlgorithm;
class G4VEMDataSet;
class G4Material;
class G4Element;

class G4CrossSectionHandler : public G4VCrossSectionHandler {
 
public:

  G4CrossSectionHandler();

  ~G4CrossSectionHandler();
	
   
protected: 
   
  virtual std::vector<G4VEMDataSet*>* BuildCrossSectionsForMaterials(const G4DataVector& energyVector, 
								       const G4DataVector* energyCuts = 0);
 
private:
 
  // Hide copy constructor and assignment operator
  G4CrossSectionHandler(const G4CrossSectionHandler&);
  G4CrossSectionHandler & operator=(const G4CrossSectionHandler &right);

};
 
#endif












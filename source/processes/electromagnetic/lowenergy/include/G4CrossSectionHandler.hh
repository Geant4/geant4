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
// $Id: G4CrossSectionHandler.hh,v 1.7 2002-05-28 09:15:26 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#ifndef G4CROSSSECTIONHANDLER_HH
#define G4CROSSSECTIONHANDLER_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "g4std/map"
#include "g4std/vector"
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
   
  virtual G4std::vector<G4VEMDataSet*>* BuildCrossSectionsForMaterials(const G4DataVector& energyVector, 
								       const G4DataVector* energyCuts = 0);
 
private:
 
  // Hide copy constructor and assignment operator
  G4CrossSectionHandler(const G4CrossSectionHandler&);
  G4CrossSectionHandler & operator=(const G4CrossSectionHandler &right);

};
 
#endif












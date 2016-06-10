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
// $Id: G4eCrossSectionHandler.cc,v 1.21 2009-09-27 10:47:42 sincerti Exp $
//
// Author: Vladimir Ivanchenko
//
// History:
// -----------
// 1  Jun 2011   V.Ivanchenko  Created
//
// -------------------------------------------------------------------

// Class description:
// Management of electron ionisation cross sections per shell

// -------------------------------------------------------------------

#ifndef G4ECROSSSECTIONHANDLER_HH
#define G4ECROSSSECTIONHANDLER_HH 1

#include "G4VCrossSectionHandler.hh"
#include "globals.hh"
#include "G4VEMDataSet.hh"
#include <vector>

class G4VDataSetAlgorithm;


class G4eCrossSectionHandler : public G4VCrossSectionHandler 
{
 
public:

  G4eCrossSectionHandler(G4VDataSetAlgorithm* alg,
			 G4double emin, G4double emax, G4int nbin);

  virtual ~G4eCrossSectionHandler();
	
protected: 
   
  virtual std::vector<G4VEMDataSet*>* 
  BuildCrossSectionsForMaterials(const G4DataVector& energyVector, 
				 const G4DataVector* energyCuts = 0);
 
private:
 
  // Hide copy constructor and assignment operator
  G4eCrossSectionHandler(const G4eCrossSectionHandler&);
  G4eCrossSectionHandler & operator=(const G4eCrossSectionHandler &right);

};
 
#endif












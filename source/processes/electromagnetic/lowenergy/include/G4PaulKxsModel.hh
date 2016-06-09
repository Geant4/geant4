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
// History:
// -----------
//  21 Apr 2008   H. Abdelohauwed - 1st implementation
//  29 Apr 2009   ALF  Major Design Revision
//  29 Oct 2011   ALF Name changed fromn G4PaulKCrossSection to G4PaulKxsModel
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, Cross section, p ionisation, K shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------


#ifndef G4PaulKxsModel_hh
#define G4PaulKxsModel_hh 1

//#include "G4VDataSetAlgorithm.hh"
#include "globals.hh"
#include <map>

class G4VDataSetAlgorithm;
class G4VEMDataSet;

class G4PaulKxsModel 

{
public:

  G4PaulKxsModel();

  virtual ~G4PaulKxsModel();
			     
  G4double CalculateKCrossSection(G4int zTarget,G4double massIncident, G4double energyIncident);
				    
 
private:


  G4PaulKxsModel(const G4PaulKxsModel&);
  G4PaulKxsModel & operator = (const G4PaulKxsModel &right);

  G4VDataSetAlgorithm* interpolation;

  std::map< G4int , G4VEMDataSet* > protonDataSetMap;

  std::map< G4int , G4VEMDataSet* > alphaDataSetMap;
  

};

#endif

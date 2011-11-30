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
// History:
// -----------
//  01 Oct 2011   A.M., S.I. - 1st implementation
// 
// Class description
// ----------------
//  Computation of K, L & M shell ECPSSR ionisation cross sections for protons and alphas
//  Based on the work of A. Taborda et al. 
//  X-Ray Spectrom. 2011, 40, 127-134
// ---------------------------------------------------------------------------------------

#ifndef G4ecpssrFormFactorKxsModel_HH
#define G4ecpssrFormFactorKxsModel_HH 1

#include "G4VecpssrKModel.hh"
#include "globals.hh"
#include <map>

class G4VDataSetAlgorithm;
class G4VEMDataSet;

class G4ecpssrFormFactorKxsModel : public G4VecpssrKModel

{
public:

  G4ecpssrFormFactorKxsModel();

  virtual ~G4ecpssrFormFactorKxsModel();
			     
  G4double CalculateCrossSection (G4int zTarget, G4double massIncident, G4double energyIncident);
				     
private:

  G4ecpssrFormFactorKxsModel(const G4ecpssrFormFactorKxsModel&);
  G4ecpssrFormFactorKxsModel & operator = (const G4ecpssrFormFactorKxsModel &right);

  G4VDataSetAlgorithm* interpolation;

  std::map< G4int , G4VEMDataSet* > protonDataSetMap;
  std::map< G4int , G4VEMDataSet* > alphaDataSetMap;

};

#endif

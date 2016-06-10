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
// $Id: G4VecpssrLiModel.hh 66241 2012-12-13 18:34:42Z gunter $
//
// Author: Haifa Ben Abdelouahed
//         
//
// History:
// -----------
//  23 Apr 2008   H. Ben Abdelouahed   1st implementation
//  28 Apr 2008   MGP        Major revision according to a design iteration
//  29 Apr 2009   ALF Updated Desing for Integration
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, Cross section, p and alpha ionisation, L shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------


#ifndef G4VECPSSRLIMODEL_HH
#define G4VECPSSRLIMODEL_HH 1

#include "globals.hh"

class G4VecpssrLiModel 
{

public:

  G4VecpssrLiModel();

  virtual ~G4VecpssrLiModel();
			     
  virtual G4double CalculateL1CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident) = 0;//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)

  virtual G4double CalculateL2CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident) = 0;//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)

  virtual G4double CalculateL3CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident) = 0;//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)
				    

   

private:


  G4VecpssrLiModel(const G4VecpssrLiModel&);
  G4VecpssrLiModel & operator = (const G4VecpssrLiModel &right);

};

#endif

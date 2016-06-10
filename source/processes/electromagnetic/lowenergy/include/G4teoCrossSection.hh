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
// $Id: G4teoCrossSection.hh 66241 2012-12-13 18:34:42Z gunter $
//
//         
//
// History:
// -----------
//  21 Apr 2008   ALF  1st implementation
//  29 Apr 2009   ALF Updated Desing for Integration
//  15 Mar 2011   ALF introduced the usage of G4AtomicShellEnumerator
//  09 Mar 2012   LP  Changed signature of methods
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, Cross section, p ionisation, K shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------


#ifndef G4TEOCROSSSECTION_HH
#define G4TEOCROSSSECTION_HH 1

#include "globals.hh"
#include "G4VhShellCrossSection.hh"

class G4VecpssrKModel;
class G4VecpssrLiModel;
class G4VecpssrMiModel;

class G4teoCrossSection : public G4VhShellCrossSection 
{
public:

  G4teoCrossSection(const G4String& name);

  virtual ~G4teoCrossSection();
			     
  std::vector<G4double> GetCrossSection(G4int Z,
					G4double incidentEnergy,
					G4double mass,
					G4double deltaEnergy = 0,
					const G4Material* mat=0);

  G4double CrossSection(G4int Z, G4AtomicShellEnumerator shell,
			G4double incidentEnergy,
			G4double mass,
			const G4Material* mat);

  std::vector<G4double> Probabilities(G4int Z,
				      G4double incidentEnergy,
				      G4double mass,
				      G4double deltaEnergy = 0,
				      const G4Material* mat=0);
  
  
  void SetTotalCS(G4double);
    
private:
  
  G4double totalCS;
              
  G4VecpssrKModel*  ecpssrShellK;
  G4VecpssrLiModel*  ecpssrShellLi;
  G4VecpssrMiModel*  ecpssrShellMi;
			
			
  G4teoCrossSection(const G4teoCrossSection&);
  G4teoCrossSection & operator = (const G4teoCrossSection &right);
  
};

#endif

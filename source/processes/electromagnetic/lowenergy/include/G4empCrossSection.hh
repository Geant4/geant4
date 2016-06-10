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
// $Id: G4empCrossSection.hh 83410 2014-08-21 15:17:53Z gcosmo $
//
//         
//
// History:
// -----------
//  21 Apr 2009   ALF  1st implementation
//  15 Mar 2011   ALF introduced the usage of G4AtomicShellEnumerator
//  09 Mar 2012   LP   changed signature of methods
//
// -------------------------------------------------------------------

// Class description:

// -------------------------------------------------------------------


#ifndef G4EMPCROSSSECTION_HH
#define G4EMPCROSSSECTION_HH 1

#include "globals.hh"
#include "G4VhShellCrossSection.hh"

#include "G4PaulKxsModel.hh"

#include "G4OrlicLiXsModel.hh"

class G4empCrossSection : public G4VhShellCrossSection 
{
public:

  G4empCrossSection(const G4String& nam = "");

  virtual ~G4empCrossSection();
			     
  std::vector<G4double> GetCrossSection(G4int Z,
					G4double incidentEnergy,
					G4double mass,
					G4double deltaEnergy,
					const G4Material* mat);

  G4double CrossSection(G4int Z, G4AtomicShellEnumerator shell,
			G4double incidentEnergy,
			G4double mass,
			const G4Material* mat);

  std::vector<G4double> Probabilities(G4int Z,
				      G4double incidentEnergy,
				      G4double mass,
				      G4double deltaEnergy,
				      const G4Material* mat);
  
  
  void SetTotalCS(G4double);
  
  
  
private:
  
  G4double totalCS;
  G4int flag; // Flag to select Li XS set (orlic or other)
              
  G4PaulKxsModel*  paulShellK;
  G4OrlicLiXsModel* orlicShellLi;  
						
  G4empCrossSection(const G4empCrossSection&);
  G4empCrossSection & operator = (const G4empCrossSection &right);
  
};

#endif

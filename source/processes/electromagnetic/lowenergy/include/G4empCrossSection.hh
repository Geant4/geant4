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
// $Id: G4empCrossSection.hh,v 1.2 2011-01-03 19:35:11 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//         
//
// History:
// -----------
//  21 Apr 2009   ALF  1st implementation
//
// -------------------------------------------------------------------

// Class description:

// -------------------------------------------------------------------


#ifndef G4EMPCROSSSECTION_HH
#define G4EMPCROSSSECTION_HH 1

#include "globals.hh"
#include "G4VhShellCrossSection.hh"
#include "G4PaulKCrossSection.hh"
#include "G4OrlicLiCrossSection.hh"

class G4empCrossSection : public G4VhShellCrossSection 
{
public:

  G4empCrossSection(const G4String& nam = "");

  virtual ~G4empCrossSection();
			     
  std::vector<G4double> GetCrossSection(G4int Z,
					G4double incidentEnergy,
					G4double mass,
					G4double deltaEnergy,
					G4bool testFlag = false) const;

  G4double CrossSection(G4int Z, G4int shell,
			G4double incidentEnergy,
			G4double mass) const;

  std::vector<G4double> Probabilities(G4int Z,
				      G4double incidentEnergy,
				      G4double mass,
				      G4double deltaEnergy) const;
  
  
  void SetTotalCS(G4double);
  
  
  
private:
  
  G4double totalCS;
              
  G4PaulKCrossSection*  paulShellK;
  G4OrlicLiCrossSection* orlicShellLi;  
						
  G4empCrossSection(const G4empCrossSection&);
  G4empCrossSection & operator = (const G4empCrossSection &right);
  
};

#endif

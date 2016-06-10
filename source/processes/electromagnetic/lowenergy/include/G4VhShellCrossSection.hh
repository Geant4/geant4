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
// $Id: G4UAtomicDeexcitation.cc,v 1.11 
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4VhShellCrossSection
//
// Author:      S. Dussoni and A. Mantero (Alfonso.Mantero@ge.infn.it)
// 
// History:
// -----------
// 23 Oct 2001 A. Mantero   1st implementation
// 24 Oct 2001 MGP          Cleaned up
// 29 Oct 2001 VI           Add delta energy
// 15 Mar 2011 ALF          Introduced the usage of G4AtomicShellEnumerator
// 09 Mar 2012 LP           Added const G4Material* to the signature of virtual 
//                          methods.
//
// -------------------------------------------------------------------

// Class Description: 
//
// Abstract class for models of shell cross sections in proton ionisation
// -------------------------------------------------------------------
//

#ifndef G4VHSHELLCROSSSECTION_HH
#define G4VHSHELLCROSSSECTION_HH 1

#include "globals.hh"
#include <vector>
#include "G4AtomicShellEnumerator.hh"
#include "G4Material.hh"


class G4VhShellCrossSection 
{

public:

  G4VhShellCrossSection(const G4String& xname = "");

  virtual ~G4VhShellCrossSection();

  G4int SelectRandomShell(G4int Z, 
                          G4double incidentEnergy,
			  G4double mass, 
			  G4double deltaEnergy,
			  const G4Material* mat);

  virtual std::vector<G4double> GetCrossSection(G4int Z,
						G4double incidentEnergy,
						G4double mass,
						G4double deltaEnergy,
						const G4Material* mat) = 0;


  virtual G4double CrossSection(G4int Z,
                                G4AtomicShellEnumerator shell,
				G4double incidentEnergy,
				G4double mass,
				const G4Material* mat) =0;

  //protected:

  virtual std::vector<G4double> Probabilities(G4int Z,
					      G4double incidentEnergy,
					      G4double mass,
					      G4double deltaEnergy,
					      const G4Material* mat) = 0;


  virtual void SetTotalCS(G4double);

  inline const G4String& GetName() const; 

private:

  // Hide copy constructor and assignment operator 
  G4VhShellCrossSection(const  G4VhShellCrossSection&);
  G4VhShellCrossSection & operator=(const  G4VhShellCrossSection &right);

  G4String name;

};

inline const G4String& G4VhShellCrossSection::GetName() const
{
  return name;
} 

#endif


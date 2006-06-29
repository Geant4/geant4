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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4hShellCrossSectionDoubleExp   
//
// Author:     Simona Saliceti (simona.saliceti@ge.infn.it)
// 
// History:
// -----------
// 31/05/2004 Simona Saliceti 1st implementation
// -------------------------------------------------------------------
// Class Description: 
// Empiric Model for K shell cross sections in proton ionisation
// 2nd iteration, refined model at low energy
// Based on Fitted empirical Reference Cross Sections for K-Shell Ionization 
// by Protons by H.Paul and J.Sacher, At. Data Nucl. Data Tables 42 (1989) 105
// 
// -------------------------------------------------------------------
// $Id: G4hShellCrossSectionDoubleExp.hh,v 1.3 2006-06-29 19:38:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4HSHELLCROSSSECTIONDOUBLEEXP_HH
#define G4HSHELLCROSSSECTIONDOUBLEEXP_HH 1

#include "globals.hh"
#include "G4VhShellCrossSection.hh" 
#include "G4hShellCrossSectionDoubleExpData.hh"

class G4hShellCrossSectionDoubleExp : public G4VhShellCrossSection
{
public:

  G4hShellCrossSectionDoubleExp();

  virtual ~G4hShellCrossSectionDoubleExp();

  virtual std::vector<G4double> GetCrossSection(G4int Z,
						G4double incidentEnergy,
						G4double mass,
						G4double deltaEnergy,
						G4bool testFlag = false) const;
 
  void SetTotalCS(G4double);
 
protected:

  virtual std::vector<G4double> Probabilities(G4int Z,
					      G4double incidentEnergy,
					      G4double mass,
					      G4double deltaEnergy) const;


private:

  G4double GetCrossSectionDoubleExp(G4int Z,
			            G4double incidentEnergy) const;

  // Hide copy constructor and assignment operator 
  G4hShellCrossSectionDoubleExp(const G4hShellCrossSectionDoubleExp&);
  G4hShellCrossSectionDoubleExp & operator = (const G4hShellCrossSectionDoubleExp &right);
  G4hShellCrossSectionDoubleExpData* kShellData;

  G4double atomTotalCrossSection;
};

#endif

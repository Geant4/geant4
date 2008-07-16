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
// $Id: G4CrossSectionPsCreationChampionPartial.hh,v 1.1 2008-07-16 19:01:07 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------


#ifndef G4CROSSSECTIONPSCREATIONCHAMPIONPARTIAL_HH
#define G4CROSSSECTIONPSCREATIONCHAMPIONPARTIAL_HH 1
 
#include <map>
#include "G4DNACrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Positron.hh"
#include "Randomize.hh"

class G4CrossSectionPsCreationChampionPartial
{
public:
  
  G4CrossSectionPsCreationChampionPartial();
  
  virtual ~G4CrossSectionPsCreationChampionPartial();
  
  G4double CrossSection(G4double energy, G4int level, const G4ParticleDefinition* particle);

  G4double Sum(G4double energy, const G4ParticleDefinition* particle);

  G4int RandomSelectState(G4double energy, const G4ParticleDefinition* particle);

  G4int RandomSelectShell(G4double energy, const G4ParticleDefinition* particle, G4int state);
   
private:
   
  G4int numberOfPartialCrossSections;

  typedef std::map<G4String,G4String,std::less<G4String> > MapFile;
  MapFile tableFile;

  typedef std::map<G4String,G4DNACrossSectionDataSet*,std::less<G4String> > MapData;
  MapData tableData;
};
#endif

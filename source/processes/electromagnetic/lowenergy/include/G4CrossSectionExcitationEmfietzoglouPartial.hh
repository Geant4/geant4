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
// $Id: G4CrossSectionExcitationEmfietzoglouPartial.hh,v 1.2 2008-07-14 20:47:34 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4CROSSSECTIONEXCITATIONEMFIETZOGLOUPARTIAL_HH
#define G4CROSSSECTIONEXCITATIONEMFIETZOGLOUPARTIAL_HH 1
 
#include "G4WaterExcitationStructure.hh"
#include <deque>
#include "Randomize.hh"

class G4Track;
 
class G4CrossSectionExcitationEmfietzoglouPartial
{
public:
  
  G4CrossSectionExcitationEmfietzoglouPartial();
  
  virtual ~G4CrossSectionExcitationEmfietzoglouPartial();
  
  G4double CrossSection(G4double energy,G4int level);

  G4double Sum(G4double energy);

  G4int RandomSelect(G4double energy);
  
private:
   
  G4int nLevels;

  G4WaterExcitationStructure waterExcitation;
  
};

#endif

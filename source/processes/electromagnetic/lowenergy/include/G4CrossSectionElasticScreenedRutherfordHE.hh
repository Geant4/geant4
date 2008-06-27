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
// $Id: G4CrossSectionElasticScreenedRutherfordHE.hh,v 1.1 2008-06-27 20:09:54 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4CROSSSECTIONELASTICSCREENEDRUTHERFORDHE_HH
#define G4CROSSSECTIONELASTICSCREENEDRUTHERFORDHE_HH 1
 
#include "G4Track.hh"

class G4CrossSectionElasticScreenedRutherfordHE
{
public:
  
  G4CrossSectionElasticScreenedRutherfordHE();
  
  virtual ~G4CrossSectionElasticScreenedRutherfordHE();
  
  G4double CrossSection(const G4Track&);
    
private:
  
  G4double RutherfordCrossSection(G4double energy, G4double z);
  
  G4double ScreeningFactor(G4double energy, G4double z);
  
  G4double lowEnergyLimit;
  G4double highEnergyLimit;
};

#endif

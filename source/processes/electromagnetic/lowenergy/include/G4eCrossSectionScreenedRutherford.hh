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
//
// $Id: G4eCrossSectionScreenedRutherford.hh,v 1.2 2007-10-12 12:26:34 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Geant4-DNA Cross total cross section for electron elastic scattering in water
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#ifndef G4ECROSSSECTIONSCREENEDRUTHERFORD_HH
#define G4ECROSSSECTIONSCREENEDRUTHERFORD_HH 1
 
#include "globals.hh"
//#include "G4Track.hh"

class G4Track;
 
class G4eCrossSectionScreenedRutherford
{
public:
  
  G4eCrossSectionScreenedRutherford();
  
  virtual ~G4eCrossSectionScreenedRutherford();
  
  G4double CrossSection(const G4Track&);
  
  // Copy constructor and assignment operator to be added here
    
private:
  
  G4double RutherfordCrossSection(G4double energy, G4double z);
  
  G4double ScreeningFactor(G4double energy, G4double z);
  
  G4String name;  
  G4double lowEnergyLimit;
  G4double highEnergyLimit;
  
};

#endif

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
// $Id: G4OrlicLCrossSection.hh,v 1.3 2008/06/03 07:29:15 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Author: Haifa Ben Abdelouahed
//         
//
// History:
// -----------
//  23 Apr 2008   H. Ben Abdelouahed   1st implementation
//  28 Apr 2008   MGP        Major revision according to a design iteration
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, Cross section, proton ionisation, L shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------


#ifndef G4ORLICLCROSSSECTION_HH
#define G4ORLICLCROSSSECTION_HH 1

#include "globals.hh"


class G4OrlicLCrossSection 

{
public:

  G4OrlicLCrossSection();

  virtual ~G4OrlicLCrossSection();


//according to I.ORLIC, C.H.SOW and S.M.TANG,International Journal of PIXE.Vol.4(1997) 217-230	
		     
  G4double CalculateL1CrossSection(G4int zTarget, G4double energyIncident);

  G4double CalculateL2CrossSection(G4int zTarget, G4double energyIncident);				    
 
  G4double CalculateL3CrossSection(G4int zTarget, G4double energyIncident);

private:


  G4OrlicLCrossSection(const G4OrlicLCrossSection&);
  G4OrlicLCrossSection & operator = (const G4OrlicLCrossSection &right);

};

#endif

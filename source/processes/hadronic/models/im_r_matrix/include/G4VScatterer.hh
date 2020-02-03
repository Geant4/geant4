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
//
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:      G4VScatterer.hh
//
//      Author:         First implementation by an unknown person
//                      Completely rewritten by Maria Grazia Pia
// 
//      Creation date:  Unknown
//
//      Modifications: 
//      27-10-99 MGP    Removed old Dubna stuff 
//
// Class Description: 
//
// Abstract base class for the Scatterer
//
// Class Description: End 
//      
// -------------------------------------------------------------------

#ifndef G4VSCATTERER_HH
#define G4VSCATTERER_HH

#include <vector>

#include "globals.hh"
#include "G4KineticTrackVector.hh"

class G4KineticTrack;

class G4VScatterer 
{
public:
  
  G4VScatterer() {}
  
  virtual ~G4VScatterer() {}
  
  virtual G4double GetTimeToInteraction(const G4KineticTrack& trk1, 
					const G4KineticTrack& trk2) const = 0;
  
  virtual G4KineticTrackVector* Scatter(const G4KineticTrack& trk1, 	
					   const G4KineticTrack& trk2) const = 0;
  
};

#endif 



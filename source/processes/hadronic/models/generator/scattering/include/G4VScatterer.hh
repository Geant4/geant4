//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4VScatterer.hh,v 1.3 2003/06/16 17:08:47 gunter Exp $ //
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

#include "globals.hh"
#include <vector>
#include "G4KineticTrackVector.hh"

class G4KineticTrack;

class G4VScatterer 
{
public:
  
  G4VScatterer() {}
  
  virtual ~G4VScatterer() {}
  
  virtual G4double GetTimeToInteraction(const G4KineticTrack& trk1, 
					const G4KineticTrack& trk2) = 0;
  
  virtual G4KineticTrackVector* Scatter(const G4KineticTrack& trk1, 	
					   const G4KineticTrack& trk2) = 0;
  
};

#endif 



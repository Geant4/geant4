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
//
// $Id: G4PSTermination3D.hh,v 1.2 2007-08-28 10:11:29 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4PSTermination3D_h
#define G4PSTermination3D_h 1

#include "G4PSTermination.hh"
//////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring Number of tracks
//  which stops(terminated) in the cell.
// 
// Created: 2007-08-14  Tsukasa ASO
// 
///////////////////////////////////////////////////////////////////////////////


class G4PSTermination3D : public G4PSTermination
{
 
 public: // with description
      G4PSTermination3D(G4String name, 
			G4int ni=1,G4int nj=1, G4int nk=1,
			G4int di=2, G4int dj=1, G4int dk=0);
      virtual ~G4PSTermination3D();

 protected: // with description
      virtual G4int GetIndex(G4Step*);

 private:
      G4int fDepthi, fDepthj, fDepthk;

};



#endif

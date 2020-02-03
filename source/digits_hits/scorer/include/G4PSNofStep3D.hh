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

#ifndef G4PSNofStep3D_h
#define G4PSNofStep3D_h 1

#include "G4PSNofStep.hh"
///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring Number of Steps in the cell.
//
// Created: 2007-08-14  Tsukasa ASO
// 2010-07-22   Introduce Unit specification.
//
///////////////////////////////////////////////////////////////////////////////


class G4PSNofStep3D : public G4PSNofStep
{
   public: // with description
      G4PSNofStep3D(G4String name,
		    G4int ni=1,G4int nj=1, G4int nk=1,
		    G4int depi=2, G4int depj=1, G4int depk=0);
      virtual ~G4PSNofStep3D();

  protected: // with description
      virtual G4int GetIndex(G4Step*);

  private:
      G4int fDepthi, fDepthj, fDepthk;
};
#endif


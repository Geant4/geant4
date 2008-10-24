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
// $Id: G4ErrorPropagationNavigator.hh,v 1.2 2008-10-24 14:00:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
//      GEANT 4 class header file 
// --------------------------------------------------------------------
//
// Class Description:
//
// Class for performing double navigation in the detector geometry and
// on the target surface for error propagation. It overloads ComputeStep()
// and ComputeSafety() methods.

// History:
// - Created. P. Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorPropagationNavigator_hh
#define G4ErrorPropagationNavigator_hh

#include "G4Navigator.hh"
#include "G4ThreeVector.hh"

class G4ErrorPropagationNavigator : public G4Navigator
{
  public:  // with description

    G4ErrorPropagationNavigator();
   ~G4ErrorPropagationNavigator();
  
    G4double ComputeStep (const G4ThreeVector &pGlobalPoint,
                          const G4ThreeVector &pDirection,
                          const G4double pCurrentProposedStepLength,
                          G4double &pNewSafety);
      // Calls the navigation in the detector geometry and then checks
      // if the distance to surface is smaller than the proposed step

    G4double ComputeSafety(const G4ThreeVector &globalpoint,
                           const G4double pProposedMaxLength = DBL_MAX,
                           const G4bool keepState = false);
      // Calls the navigation in the detector geometry and then checks
      // if the distance to surface is smaller than the proposed safety

};

#endif

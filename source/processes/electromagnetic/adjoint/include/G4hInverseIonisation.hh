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
// $Id: G4hInverseIonisation.hh 66892 2013-01-17 10:57:59Z gunter $
//
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4hInverseIonisation.hh
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	15 February 2009 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint/reverse discrete ionisation for proton
//

#ifndef G4hInverseIonisation_h
#define G4hInverseIonisation_h 1

#include "G4VAdjointReverseReaction.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
#include "G4AdjointhIonisationModel.hh"
class G4hInverseIonisation: public G4VAdjointReverseReaction

{
public:

  G4hInverseIonisation(G4bool whichScatCase, G4String process_name, G4AdjointhIonisationModel* aEmAdjointModel);
  ~G4hInverseIonisation();
  
private:
    
};

#endif

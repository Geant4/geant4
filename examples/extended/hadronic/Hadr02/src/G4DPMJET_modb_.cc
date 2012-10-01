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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/src/G4DPMJET_modb_.cc
/// \brief Implementation of the G4DPMJET_modb_ class
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4DPMJET_modb_.cc
//
// Version:             0.B
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//

#ifdef G4_USE_DPMJET

#include "globals.hh"
#include "G4GlaubAADataSetHandler.hh"

extern "C" void g4dpmjet_modb_ (G4double &v, G4double &b)
{
  G4GlaubAADataSetHandler *theGlauberDataSetHandler = 
    G4GlaubAADataSetHandler::getInstance();
  G4ParamType1GlaubAADataSet *theCurrentGlauberDataSet = 
    (G4ParamType1GlaubAADataSet*)
    theGlauberDataSetHandler->GetCurrentGlauberDataSet();

  G4double i   = theGlauberDataSetHandler->GetValueN(v);
  b            = i * theCurrentGlauberDataSet->bstep;
//
//
// The following lines if you require information about the
// impact parameter sampled by theGlauberDataSetHandler.
//
#ifdef G4DPMJET25GDSHDEBUG
  if (theGlauberDataSetHandler->GetVerboseLevel()>=2) {
    G4cout <<"****************************************"
           <<"****************************************"
           <<G4endl;
    G4cout <<"IN G4DPMJET_modb_" <<G4endl;
    G4cout <<"PROJECTILE A     = " <<theCurrentGlauberDataSet->GetAP()
           <<" & TARGET A      = " <<theCurrentGlauberDataSet->GetAT() <<G4endl;
    G4cout <<"RANDOM NUMBER    = " <<v <<G4endl;
    G4cout <<"IMPACT PARAMETER = " <<b <<G4endl;
    G4cout <<"****************************************"
           <<"****************************************"
           <<G4endl;
  }
#endif
  
}
#endif

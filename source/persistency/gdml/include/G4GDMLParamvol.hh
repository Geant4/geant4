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
// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLPARAMVOL_INCLUDED_
#define _G4GDMLPARAMVOL_INCLUDED_

#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"

#include "G4GDMLParameterisation.hh"
#include "G4GDMLSetup.hh"

class G4GDMLParamvol : public G4GDMLSetup {

   G4GDMLParameterisation* parameterisation;

   void box_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void trd_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void trap_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void tube_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void cone_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void sphere_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void orb_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void torus_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void para_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
   void hype_dimensionsRead(const xercesc::DOMElement* const,G4GDMLParameterisation::PARAMETER&);
  
   void parametersRead(const xercesc::DOMElement* const);
   void paramvol_contentRead(const xercesc::DOMElement* const);
protected:
   void paramvolRead(const xercesc::DOMElement* const,G4LogicalVolume*);
};

#endif

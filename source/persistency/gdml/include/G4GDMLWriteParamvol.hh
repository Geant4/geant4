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
// $Id: G4GDMLWriteParamvol.hh,v 1.7 2008-07-11 07:50:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4GDMLWriteParamvol
//
// Class description:
//
// GDML class for writing parameterised entities dimensions.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLWRITEPARAMVOL_INCLUDED_
#define _G4GDMLWRITEPARAMVOL_INCLUDED_

#include "G4PVParameterised.hh"
#include "G4VPhysicalVolume.hh"

#include "G4GDMLWriteSetup.hh"

class G4GDMLWriteParamvol : public G4GDMLWriteSetup
{
 private:

   void box_dimensionsWrite(xercesc::DOMElement*, const G4Box* const);
   void trd_dimensionsWrite(xercesc::DOMElement*, const G4Trd* const);
   void trap_dimensionsWrite(xercesc::DOMElement*, const G4Trap* const);
   void tube_dimensionsWrite(xercesc::DOMElement*, const G4Tubs* const);
   void cone_dimensionsWrite(xercesc::DOMElement*, const G4Cons* const);
   void sphere_dimensionsWrite(xercesc::DOMElement*, const G4Sphere* const);
   void orb_dimensionsWrite(xercesc::DOMElement*, const G4Orb* const);
   void torus_dimensionsWrite(xercesc::DOMElement*, const G4Torus* const);
   void para_dimensionsWrite(xercesc::DOMElement*, const G4Para* const);
   void hype_dimensionsWrite(xercesc::DOMElement*, const G4Hype* const);
   void parametersWrite(xercesc::DOMElement*,
                        const G4VPhysicalVolume* const, const G4int&);

 protected:

   void paramvolWrite(xercesc::DOMElement*, const G4VPhysicalVolume* const);
};

#endif

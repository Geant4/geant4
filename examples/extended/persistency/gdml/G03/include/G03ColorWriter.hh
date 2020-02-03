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
/// \file persistency/gdml/G03/include/G03ColorWriter.hh
/// \brief Definition of the G03ColorWriter class
//
//
//
//
// class G03ColorWriter
//
// Custom writer for handling "color" tags extensions in GDML.
// -------------------------------------------------------------------------

#ifndef G03ColorWriter_H
#define G03ColorWriter_H 1

#include <vector>
#include "G4GDMLWriteStructure.hh"

class G4LogicalVolume;
class G4VisAttributes;

/// GDML writer for the color attributes

class G03ColorWriter : public G4GDMLWriteStructure
{

 public:

   G03ColorWriter();
  ~G03ColorWriter();

   void AddExtension(xercesc::DOMElement* volumeElement,
                     const G4LogicalVolume* const vol);
   void ExtensionWrite(xercesc::DOMElement* element);
   void ColorWrite(xercesc::DOMElement* volumeElement,
                   const G4VisAttributes* const att);

   G4bool BookAttribute(const G4VisAttributes* const att);

 private:

   std::vector<const G4VisAttributes*> fAttribs;
};

#endif

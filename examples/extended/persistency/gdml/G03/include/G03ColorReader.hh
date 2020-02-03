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
/// \file persistency/gdml/G03/include/G03ColorReader.hh
/// \brief Definition of the G03ColorReader class
//
//
//
//
// class G03ColorReader
//
// Custom reader for handling "color" tags extensions in GDML.
// -------------------------------------------------------------------------

#ifndef G03ColorReader_H
#define G03ColorReader_H 1

#include <map>
#include "G4GDMLReadStructure.hh"

class G4VisAttributes;

/// GDML reader for the color attributes

class G03ColorReader : public G4GDMLReadStructure
{

 public:

   G03ColorReader();
  ~G03ColorReader();

   void ExtensionRead(const xercesc::DOMElement* const element);
   void ColorRead(const xercesc::DOMElement* const element);

   G4VisAttributes* GetVisAttribute(const G4String& ref);

 protected:

   virtual void VolumeRead(const xercesc::DOMElement* const);

 private:

   std::map<G4String, G4VisAttributes*> fAttribs;
};

#endif

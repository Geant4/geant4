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
// $Id: G4GDMLRead.hh,v 1.23 2008-08-22 09:35:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4GDMLRead
//
// Class description:
//
// GDML reader.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLBASE_INCLUDED_
#define _G4GDMLBASE_INCLUDED_

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/dom/DOM.hpp>

#include "G4GDMLEvaluator.hh"

#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"

#include <sstream>

class G4GDMLErrorHandler : public xercesc::ErrorHandler
{
   G4bool Suppress;

 public:

   G4GDMLErrorHandler(const G4bool set) { Suppress = set; }

   void warning(const xercesc::SAXParseException& exception)
   {
      if (Suppress)  { return; }
      char* message = xercesc::XMLString::transcode(exception.getMessage());
      G4cout << "G4GDML: VALIDATION WARNING! " << message
             << " at line: " << exception.getLineNumber() << G4endl;
      xercesc::XMLString::release(&message);
   }

   void error(const xercesc::SAXParseException& exception)
   {
      if (Suppress)  { return; }
      char* message = xercesc::XMLString::transcode(exception.getMessage());
      G4cout << "G4GDML: VALIDATION ERROR! " << message
             << " at line: " << exception.getLineNumber() << G4endl;
      xercesc::XMLString::release(&message);
   }

   void fatalError(const xercesc::SAXParseException& exception)
   {
      error(exception);
   }
   void resetErrors() {}
};

class G4GDMLRead
{

 public:

   virtual void DefineRead(const xercesc::DOMElement* const)=0;
   virtual void MaterialsRead(const xercesc::DOMElement* const)=0;
   virtual void SetupRead(const xercesc::DOMElement* const)=0;
   virtual void SolidsRead(const xercesc::DOMElement* const)=0;
   virtual void Paramvol_contentRead(const xercesc::DOMElement* const)=0;
   virtual void Volume_contentRead(const xercesc::DOMElement* const)=0;
   virtual void StructureRead(const xercesc::DOMElement* const)=0;
   virtual G4LogicalVolume* GetVolume(const G4String&) const=0;
   virtual G4String GetSetup(const G4String&)=0;
   void Read(const G4String&, G4bool SetValidate, G4bool IsModule);

 protected:

   G4GDMLEvaluator eval;
   G4bool Validate;

   G4String Transcode(const XMLCh* const);
   G4String GenerateName(const G4String& name, G4bool strip=false);
   G4String GenerateMatName(const G4String& name, G4bool strip=false);
   G4String Strip(const G4String&) const;
   void StripName(G4String&) const;
   void GeneratePhysvolName(const G4String&,G4VPhysicalVolume*);
   void LoopRead(const xercesc::DOMElement* const,
                 void(G4GDMLRead::*)(const xercesc::DOMElement* const));
 private:

   void StripNames() const;

 private:

   G4String ModuleName;
   G4int InLoop;

};

#endif

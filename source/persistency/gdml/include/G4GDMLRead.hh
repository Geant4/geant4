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
// $Id: G4GDMLRead.hh,v 1.9 2008-06-10 09:12:45 ztorzsok Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4GDMLBase
//
// Class description:
//
// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLBASE_INCLUDED_
#define _G4GDMLBASE_INCLUDED_

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/dom/DOM.hpp>

#include <sstream>

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

#include "G4GDMLEvaluator.hh"

class G4GDMLRead {
private:
   G4String prename;
   G4int InLoop;
protected:
   G4GDMLEvaluator eval;

   G4String GenerateName(const G4String&);
   void loopRead(const xercesc::DOMElement* const,void(G4GDMLRead::*)(const xercesc::DOMElement* const));
public:
   virtual void defineRead(const xercesc::DOMElement* const)=0;
   virtual void materialsRead(const xercesc::DOMElement* const)=0;
   virtual void setupRead(const xercesc::DOMElement* const)=0;
   virtual void solidsRead(const xercesc::DOMElement* const)=0;
   virtual void paramvol_contentRead(const xercesc::DOMElement* const)=0;
   virtual void volume_contentRead(const xercesc::DOMElement* const)=0;
   virtual void structureRead(const xercesc::DOMElement* const)=0;
   virtual G4LogicalVolume* getVolume(const G4String&) const=0;
   virtual G4String getSetup(const G4String&)=0;
   void Read(const G4String&,bool external=false);
};

#endif

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
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLWriter.hh"

void G4GDMLWriter::Write(xercesc::DOMNode* node) {

   XMLCh tempStr[100];

   xercesc::XMLString::transcode("LS", tempStr, 99);
   xercesc::DOMImplementation* impl = xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
   xercesc::DOMWriter* theSerializer = ((xercesc::DOMImplementationLS*)impl)->createDOMWriter();

   xercesc::XMLFormatTarget *myFormTarget = new xercesc::StdOutFormatTarget();

   try {
      theSerializer->writeNode(myFormTarget,*node);
   } catch (const xercesc::XMLException& toCatch) {
   
      char* message = xercesc::XMLString::transcode(toCatch.getMessage());
      G4cout << "Exception message is: \n" << message << "\n" << G4endl;
      xercesc::XMLString::release(&message);
      return;
   } catch (const xercesc::DOMException& toCatch) {
   
      char* message = xercesc::XMLString::transcode(toCatch.msg);
      G4cout << "Exception message is: \n" << message << "\n" << G4endl;
      xercesc::XMLString::release(&message);
      return;
   } catch (...) {
      
      G4cout << "Unexpected Exception \n" << G4endl;
      return;
   }        

   theSerializer->release();
}

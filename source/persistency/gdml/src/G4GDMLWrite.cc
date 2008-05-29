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

#include <sstream>
#include "G4GDMLWrite.hh"

G4GDMLWrite::ModuleMapType& moduleMap() {

   static G4GDMLWrite::ModuleMapType instance;
   return instance;
}

xercesc::DOMAttr* G4GDMLWrite::newAttribute(const G4String& name,const G4String& value) {

   xercesc::XMLString::transcode(name,tempStr,99);
   xercesc::DOMAttr* att = doc->createAttribute(tempStr);
   xercesc::XMLString::transcode(value,tempStr,99);
   att->setValue(tempStr);
   return att;
}

xercesc::DOMAttr* G4GDMLWrite::newAttribute(const G4String& name,const G4double& value) {

   xercesc::XMLString::transcode(name,tempStr,99);
   xercesc::DOMAttr* att = doc->createAttribute(tempStr);
   std::ostringstream ostream;
   ostream << value;
   G4String str = ostream.str();
   xercesc::XMLString::transcode(str,tempStr,99);
   att->setValue(tempStr);
   return att;
}

xercesc::DOMElement* G4GDMLWrite::newElement(const G4String& name) {

   xercesc::XMLString::transcode(name,tempStr,99);
   return doc->createElement(tempStr);
}

G4Transform3D G4GDMLWrite::Write(const G4String& fname,const G4LogicalVolume* const logvol) {

   G4cout << "Writing '" << fname << "'..." << G4endl;

   xercesc::XMLString::transcode("LS", tempStr, 99);
   xercesc::DOMImplementation* impl = xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
   xercesc::DOMWriter* writer = ((xercesc::DOMImplementationLS*)impl)->createDOMWriter();

   xercesc::XMLString::transcode("Range", tempStr, 99);
   impl = xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);

   xercesc::XMLString::transcode("gdml", tempStr, 99);
   doc = impl->createDocument(0,tempStr,0);
   xercesc::DOMElement* gdml = doc->getDocumentElement();

   if (writer->canSetFeature(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint,true))
      writer->setFeature(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint,true);

   defineWrite(gdml);
   materialsWrite(gdml);
   solidsWrite(gdml);
   structureWrite(gdml);
   setupWrite(gdml,logvol->GetName());

   G4Transform3D R = TraverseVolumeTree(logvol);

   xercesc::XMLFormatTarget *myFormTarget = new xercesc::LocalFileFormatTarget(fname.c_str());

   try {
      writer->writeNode(myFormTarget,*doc);
   } catch (const xercesc::XMLException& toCatch) {
   
      char* message = xercesc::XMLString::transcode(toCatch.getMessage());
      G4cout << "Exception message is: \n" << message << "\n" << G4endl;
      xercesc::XMLString::release(&message);
      return G4Transform3D::Identity;
   } catch (const xercesc::DOMException& toCatch) {
   
      char* message = xercesc::XMLString::transcode(toCatch.msg);
      G4cout << "Exception message is: \n" << message << "\n" << G4endl;
      xercesc::XMLString::release(&message);
      return G4Transform3D::Identity;
   } catch (...) {
      
      G4cout << "Unexpected Exception \n" << G4endl;
      return G4Transform3D::Identity;
   }        

   delete myFormTarget;
   writer->release();

   G4cout << "Writing done!" << G4endl << G4endl;

   return R;
}

void G4GDMLWrite::SetModule(const G4VPhysicalVolume* const physvol,const G4String& name) {

   moduleMap()[physvol] = name;
}

G4String G4GDMLWrite::GetModule(const G4VPhysicalVolume* const physvol) {

   if (moduleMap().find(physvol) != moduleMap().end()) return moduleMap()[physvol];
   return G4String("");
}

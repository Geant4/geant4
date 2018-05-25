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
// $Id: G4GDMLWrite.cc 110108 2018-05-15 11:46:54Z gcosmo $
//
// class G4GDMLWrite Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include <sys/stat.h>
#include <iostream>

#include "G4GDMLWrite.hh"

#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4PVDivision.hh"

G4bool G4GDMLWrite::addPointerToName = true;

G4GDMLWrite::G4GDMLWrite() : doc(0), extElement(0)
{
}

G4GDMLWrite::~G4GDMLWrite()
{
}

G4bool G4GDMLWrite::FileExists(const G4String& fname) const
{
  struct stat FileInfo;
  return (stat(fname.c_str(),&FileInfo) == 0); 
}

G4GDMLWrite::VolumeMapType& G4GDMLWrite::VolumeMap()
{
   static VolumeMapType instance;
   return instance;
}

G4GDMLWrite::PhysVolumeMapType& G4GDMLWrite::PvolumeMap()
{
   static PhysVolumeMapType instance;
   return instance;
}

G4GDMLWrite::DepthMapType& G4GDMLWrite::DepthMap()
{
   static DepthMapType instance;
   return instance;
}

void G4GDMLWrite::AddExtension(xercesc::DOMElement*,
                               const G4LogicalVolume* const)
{
   // Empty implementation. To be overwritten by user for specific extensions
   // related to attributes associated to volumes
}

void G4GDMLWrite::ExtensionWrite(xercesc::DOMElement*)
{
   // Empty implementation. To be overwritten by user for specific extensions
}

void G4GDMLWrite::AddAuxInfo(G4GDMLAuxListType* auxInfoList,
                             xercesc::DOMElement* element)
{
  for(std::vector<G4GDMLAuxStructType>::const_iterator
      iaux = auxInfoList->begin(); iaux != auxInfoList->end(); iaux++ )
  {
    xercesc::DOMElement* auxiliaryElement = NewElement("auxiliary");
    element->appendChild(auxiliaryElement);
      
    auxiliaryElement->setAttributeNode(NewAttribute("auxtype", (*iaux).type));
    auxiliaryElement->setAttributeNode(NewAttribute("auxvalue", (*iaux).value));
    if (((*iaux).unit)!="")
    {
      auxiliaryElement->setAttributeNode(NewAttribute("auxunit", (*iaux).unit));
    }

    if (iaux->auxList)  { AddAuxInfo(iaux->auxList, auxiliaryElement); }
  }
  return;
}

void G4GDMLWrite::UserinfoWrite(xercesc::DOMElement* gdmlElement)
{
  if(auxList.size()>0)
  {
#ifdef G4VERBOSE
    G4cout << "G4GDML: Writing userinfo..." << G4endl;
#endif
    userinfoElement = NewElement("userinfo");
    gdmlElement->appendChild(userinfoElement);
    AddAuxInfo(&auxList, userinfoElement);
  }
}

G4String G4GDMLWrite::GenerateName(const G4String& name, const void* const ptr)
{
   G4String nameOut;
   std::stringstream stream; stream << name;
   if (addPointerToName) { stream << ptr; };

   nameOut=G4String(stream.str());
   if(nameOut.contains(' '))
   nameOut.erase(std::remove(nameOut.begin(),nameOut.end(),' '),nameOut.end());

   return nameOut;
}

xercesc::DOMAttr* G4GDMLWrite::NewAttribute(const G4String& name,
                                            const G4String& value)
{
   xercesc::XMLString::transcode(name,tempStr,9999);
   xercesc::DOMAttr* att = doc->createAttribute(tempStr);
   xercesc::XMLString::transcode(value,tempStr,9999);
   att->setValue(tempStr);
   return att;
}

xercesc::DOMAttr* G4GDMLWrite::NewAttribute(const G4String& name,
                                            const G4double& value)
{
   xercesc::XMLString::transcode(name,tempStr,9999);
   xercesc::DOMAttr* att = doc->createAttribute(tempStr);
   std::ostringstream ostream;
   ostream.precision(15);
   ostream << value;
   G4String str = ostream.str();
   xercesc::XMLString::transcode(str,tempStr,9999);
   att->setValue(tempStr);
   return att;
}

xercesc::DOMElement* G4GDMLWrite::NewElement(const G4String& name)
{
   xercesc::XMLString::transcode(name,tempStr,9999);
   return doc->createElement(tempStr);
}

G4Transform3D G4GDMLWrite::Write(const G4String& fname,
                                 const G4LogicalVolume* const logvol,
                                 const G4String& setSchemaLocation,
                                 const G4int depth,
                                       G4bool refs)
{
   SchemaLocation = setSchemaLocation;
   addPointerToName = refs;
#ifdef G4VERBOSE
   if (depth==0) { G4cout << "G4GDML: Writing '" << fname << "'..." << G4endl; }
   else   { G4cout << "G4GDML: Writing module '" << fname << "'..." << G4endl; }
#endif
   if (FileExists(fname))
   {
     G4String ErrorMessage = "File '"+fname+"' already exists!";
     G4Exception("G4GDMLWrite::Write()", "InvalidSetup",
                 FatalException, ErrorMessage);
   }
   
   VolumeMap().clear(); // The module map is global for all modules,
                        // so clear it only at once!

   xercesc::XMLString::transcode("LS", tempStr, 9999);
     xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
   xercesc::XMLString::transcode("Range", tempStr, 9999);
   xercesc::DOMImplementation* impl =
     xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
   xercesc::XMLString::transcode("gdml", tempStr, 9999);
   doc = impl->createDocument(0,tempStr,0);
   xercesc::DOMElement* gdml = doc->getDocumentElement();

#if XERCES_VERSION_MAJOR >= 3
                                             // DOM L3 as per Xerces 3.0 API
    xercesc::DOMLSSerializer* writer =
      ((xercesc::DOMImplementationLS*)impl)->createLSSerializer();

    xercesc::DOMConfiguration *dc = writer->getDomConfig();
    dc->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);

#else

   xercesc::DOMWriter* writer =
     ((xercesc::DOMImplementationLS*)impl)->createDOMWriter();

   if (writer->canSetFeature(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true))
       writer->setFeature(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);

#endif

   gdml->setAttributeNode(NewAttribute("xmlns:xsi",
                          "http://www.w3.org/2001/XMLSchema-instance"));
   gdml->setAttributeNode(NewAttribute("xsi:noNamespaceSchemaLocation",
                          SchemaLocation));

   ExtensionWrite(gdml);
   DefineWrite(gdml);
   MaterialsWrite(gdml);
   SolidsWrite(gdml);
   StructureWrite(gdml);
   UserinfoWrite(gdml);
   SetupWrite(gdml,logvol);

   G4Transform3D R = TraverseVolumeTree(logvol,depth);

   SurfacesWrite();
   xercesc::XMLFormatTarget *myFormTarget =
     new xercesc::LocalFileFormatTarget(fname.c_str());

   try
   {
#if XERCES_VERSION_MAJOR >= 3
                                            // DOM L3 as per Xerces 3.0 API
      xercesc::DOMLSOutput *theOutput =
        ((xercesc::DOMImplementationLS*)impl)->createLSOutput();
      theOutput->setByteStream(myFormTarget);
      writer->write(doc, theOutput);
#else
      writer->writeNode(myFormTarget, *doc);
#endif
   }
   catch (const xercesc::XMLException& toCatch)
   {
      char* message = xercesc::XMLString::transcode(toCatch.getMessage());
      G4cout << "G4GDML: Exception message is: " << message << G4endl;
      xercesc::XMLString::release(&message);
      return G4Transform3D::Identity;
   }
   catch (const xercesc::DOMException& toCatch)
   {
      char* message = xercesc::XMLString::transcode(toCatch.msg);
      G4cout << "G4GDML: Exception message is: " << message << G4endl;
      xercesc::XMLString::release(&message);
      return G4Transform3D::Identity;
   }
   catch (...)
   {   
      G4cout << "G4GDML: Unexpected Exception!" << G4endl;
      return G4Transform3D::Identity;
   }        

   delete myFormTarget;
   writer->release();

   if (depth==0)
   {
     G4cout << "G4GDML: Writing '" << fname << "' done !" << G4endl;
   }
   else
   {
#ifdef G4VERBOSE
     G4cout << "G4GDML: Writing module '" << fname << "' done !" << G4endl;
#endif
   }

   return R;
}

void G4GDMLWrite::AddModule(const G4VPhysicalVolume* const physvol)
{
   G4String fname = GenerateName(physvol->GetName(),physvol);
   G4cout << "G4GDML: Adding module '" << fname << "'..." << G4endl;

   if (physvol == 0)
   {
     G4Exception("G4GDMLWrite::AddModule()", "InvalidSetup", FatalException,
                 "Invalid NULL pointer is specified for modularization!");
     return;
   }
   if (dynamic_cast<const G4PVDivision*>(physvol))
   {
     G4Exception("G4GDMLWrite::AddModule()", "InvalidSetup", FatalException,
                 "It is not possible to modularize by divisionvol!");
     return;
   }
   if (physvol->IsParameterised())
   {
     G4Exception("G4GDMLWrite::AddModule()", "InvalidSetup", FatalException,
                 "It is not possible to modularize by parameterised volume!");
     return;
   }
   if (physvol->IsReplicated())
   {
     G4Exception("G4GDMLWrite::AddModule()", "InvalidSetup", FatalException,
                 "It is not possible to modularize by replicated volume!");
     return;
   }

   PvolumeMap()[physvol] = fname;
}

void G4GDMLWrite::AddModule(const G4int depth)
{
   if (depth<0)
   {
     G4Exception("G4GDMLWrite::AddModule()", "InvalidSetup", FatalException,
                 "Depth must be a positive number!");
   }
   if (DepthMap().find(depth) != DepthMap().end())
   {
     G4Exception("G4GDMLWrite::AddModule()", "InvalidSetup", FatalException,
                 "Adding module(s) at this depth is already requested!");
   }
   DepthMap()[depth] = 0;
}

G4String G4GDMLWrite::Modularize( const G4VPhysicalVolume* const physvol,
                                  const G4int depth )
{
   if (PvolumeMap().find(physvol) != PvolumeMap().end())
   {
     return PvolumeMap()[physvol]; // Modularize via physvol
   }

   if (DepthMap().find(depth) != DepthMap().end()) // Modularize via depth
   {
     std::stringstream stream;
     stream << "depth" << depth << "_module" << DepthMap()[depth] << ".gdml";
     DepthMap()[depth]++;           // There can be more modules at this depth!
     return G4String(stream.str());
   }

   return G4String(""); // Empty string for module name = no modularization
                        // was requested at that level/physvol!
}

void G4GDMLWrite::AddAuxiliary(G4GDMLAuxStructType myaux)
{
   auxList.push_back(myaux);
}

void G4GDMLWrite::SetAddPointerToName(G4bool set)
{
   addPointerToName = set;
}

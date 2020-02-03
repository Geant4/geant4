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
/// \file persistency/gdml/G03/src/G03ColorWriter.cc
/// \brief Implementation of the G03ColorWriter class
//
//
// --------------------------------------------------------------------

#include "G03ColorWriter.hh"

#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G03ColorWriter::G03ColorWriter()
  : G4GDMLWriteStructure()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G03ColorWriter::~G03ColorWriter()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03ColorWriter::AddExtension(xercesc::DOMElement* volumeElement,
                               const G4LogicalVolume* const vol)
{
   const G4VisAttributes* vis = vol->GetVisAttributes();
   if (vis)  { ColorWrite(volumeElement, vis); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03ColorWriter::ExtensionWrite(xercesc::DOMElement* element)
{
   G4cout << "G4GDML: Writing GDML extension..." << G4endl;

   // Mandatory calls the -first- time an extension is created
   //
   extElement = NewElement("extension");
   element->appendChild(extElement);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03ColorWriter::ColorWrite(xercesc::DOMElement* volumeElement,
                             const G4VisAttributes* const att)
{
   G4bool book = BookAttribute(att);
   G4Color color = att->GetColor();

   const G4String& name = GenerateName("test_color", att);
   G4double r=color.GetRed(), g=color.GetGreen(),
            b=color.GetBlue(), a=color.GetAlpha();

   if (book)
   {
     xercesc::DOMElement* colElement = NewElement("color");
     colElement->setAttributeNode(NewAttribute("name",name));
     colElement->setAttributeNode(NewAttribute("R",r));
     colElement->setAttributeNode(NewAttribute("G",g));
     colElement->setAttributeNode(NewAttribute("B",b));
     colElement->setAttributeNode(NewAttribute("A",a));
     extElement->appendChild(colElement);
   }

   xercesc::DOMElement* colorrefElement = NewElement("colorref");
   colorrefElement->setAttributeNode(NewAttribute("ref",name));
   volumeElement->appendChild(colorrefElement);

   G4cout << "Written color attribute (R,G,B,A) is: "
          << r << ", " << g << ", " << b << ", " << a << " !" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G03ColorWriter::BookAttribute(const G4VisAttributes* const ref)
{
  G4bool booking = true;
  std::vector<const G4VisAttributes*>::const_iterator pos =
     std::find(fAttribs.begin(), fAttribs.end(), ref);

  if (pos != fAttribs.end())  { booking = false; }
  else                       { fAttribs.push_back(ref); }

  return booking;
}

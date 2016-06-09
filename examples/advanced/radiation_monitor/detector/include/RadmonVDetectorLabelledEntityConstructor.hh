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
// File name:     RadmonVDetectorLabelledEntityConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLabelledEntityConstructor.hh,v 1.8 2006/06/29 16:13:16 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   Abstract class of a detector-entity constructor with label
//                and attributes
//

#ifndef   RADMONVDETECTORLABELLEDENTITYCONSTRUCTOR_HH
 #define  RADMONVDETECTORLABELLEDENTITYCONSTRUCTOR_HH
 
 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"
 #include "RadmonVDetectorEntityConstructor.hh"
 
 // Forward declaration
 class G4Material;
 class G4VisAttributes;
 class G4VSensitiveDetector;
 
 class RadmonVDetectorLabelledEntityConstructor : public RadmonVDetectorEntityConstructor, public RadmonLayoutEntityWithAttributes
 {
  public:
   inline virtual                               ~RadmonVDetectorLabelledEntityConstructor();

   inline const G4String &                      GetLabel(void) const;
   inline virtual void                          SetEntityAttribute(const G4String & attributeName, const G4String & value);
   virtual RadmonVDetectorLabelledEntityConstructor * New(void) const = 0;

   G4double                                     GetWidth(void) const;
   G4double                                     GetHeight(void) const;
   G4double                                     GetThickness(void) const;

   G4Material *                                 GetMaterial(const G4String & attributeName) const;
   G4VisAttributes *                            AllocateVisAttributes(const G4String & attributeName, const G4Material * material) const;
   
   G4VSensitiveDetector *                       AllocateSensitiveDetector(const G4String & attributeName, const G4String & defaultAttrbuteValue) const;

  protected:
   inline                                       RadmonVDetectorLabelledEntityConstructor(const G4String & label);
   
  private:
  // Hidden constructors and operators
                                                RadmonVDetectorLabelledEntityConstructor(const RadmonVDetectorLabelledEntityConstructor & copy);
   RadmonVDetectorLabelledEntityConstructor &   operator=(const RadmonVDetectorLabelledEntityConstructor & copy);

  // Private attributes
   G4String                                     entityLabel;
 };
 
 #include "RadmonVDetectorLabelledEntityConstructor.icc"
#endif /* RADMONVDETECTORLABELLEDENTITYCONSTRUCTOR_HH */

//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonVDetectorLabelledEntityConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLabelledEntityConstructor.hh,v 1.7 2006-06-28 13:50:01 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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

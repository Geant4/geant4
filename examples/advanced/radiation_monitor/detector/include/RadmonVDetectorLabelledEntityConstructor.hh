//
// File name:     RadmonVDetectorLabelledEntityConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLabelledEntityConstructor.hh,v 1.4 2005-10-24 14:52:03 capra Exp $
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
 
 class RadmonVDetectorLabelledEntityConstructor : public RadmonVDetectorEntityConstructor, public RadmonLayoutEntityWithAttributes
 {
  public:
   inline virtual                               ~RadmonVDetectorLabelledEntityConstructor();

   inline const G4String &                      GetLabel(void) const;
   inline virtual void                          SetEntityAttribute(const G4String & attributeName, const G4String & value);
   virtual RadmonVDetectorLabelledEntityConstructor * New(void) const = 0;

   G4double                                     GetAttributeAsDouble(const G4String & attributeName, double defaultValue) const;
   G4double                                     GetAttributeAsMeasure(const G4String & attributeName, const char * category, double defaultValue) const;
   G4int                                        GetAttributeAsInteger(const G4String & attributeName, G4int defaultValue) const;
   
   G4double                                     GetWidth(void) const;
   G4double                                     GetHeight(void) const;
   G4double                                     GetThickness(void) const;

   G4Material *                                 GetMaterial(const G4String & attributeName) const;
   G4VisAttributes *                            AllocateVisAttributes(const G4String & attributeName, const G4Material * material) const;

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

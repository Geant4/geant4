//
// File name:     RadmonVDetectorLabelledEntityConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLabelledEntityConstructor.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a detector-entity constructor with label
//                and attributes
//

#ifndef   RADMONVDETECTORLABELLEDENTITYCONSTRUCTOR_HH
 #define  RADMONVDETECTORLABELLEDENTITYCONSTRUCTOR_HH
 
 // Include files
 #include "RadmonDetectorLayoutEntityWithAttributes.hh"
 #include "RadmonVDetectorEntityConstructor.hh"
 
 class RadmonVDetectorLabelledEntityConstructor : public RadmonVDetectorEntityConstructor, public RadmonDetectorLayoutEntityWithAttributes
 {
  public:
   inline virtual                               ~RadmonVDetectorLabelledEntityConstructor();

   inline const G4String &                      GetLabel(void) const;
   inline virtual void                          SetEntityAttribute(const G4String & attributeName, const G4String & value);
   virtual RadmonVDetectorLabelledEntityConstructor * New(void) = 0;

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

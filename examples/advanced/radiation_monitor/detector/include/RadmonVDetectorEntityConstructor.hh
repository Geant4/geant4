//
// File name:     RadmonVDetectorEntityConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorEntityConstructor.hh,v 1.2 2005-09-19 19:42:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a detector-entity constructor
//

#ifndef   RADMONVDETECTORENTITYCONSTRUCTOR_HH
 #define  RADMONVDETECTORENTITYCONSTRUCTOR_HH
 
 // Forward declarations
 class G4LogicalVolume;
 class G4String;
 
 class RadmonVDetectorEntityConstructor
 {
  public:
   inline virtual                              ~RadmonVDetectorEntityConstructor();
    
   virtual void                                 SetEntityAttribute(const G4String & attributeName, const G4String &value) = 0;
   virtual G4LogicalVolume *                    ConstructLogicalVolume() = 0;

  protected:
   inline                                       RadmonVDetectorEntityConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonVDetectorEntityConstructor(const RadmonVDetectorEntityConstructor & copy);
   RadmonVDetectorEntityConstructor &           operator=(const RadmonVDetectorEntityConstructor & copy);

 };
 
 // Inline implementations
 #include "RadmonVDetectorEntityConstructor.icc"
#endif /* RADMONVDETECTORENTITYCONSTRUCTOR_HH */

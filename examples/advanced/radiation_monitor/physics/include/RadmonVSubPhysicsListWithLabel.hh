//
// File name:     RadmonVSubPhysicsListWithLabel.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVSubPhysicsListWithLabel.hh,v 1.2 2005-11-10 08:14:10 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a physics list piece with label
//                and attributes
//

#ifndef   RADMONVSUBPHYSICSLISTWITHLABEL_HH
 #define  RADMONVSUBPHYSICSLISTWITHLABEL_HH
 
 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"
 #include "RadmonVSubPhysicsList.hh"
 
 class RadmonVSubPhysicsListWithLabel : public RadmonVSubPhysicsList, public RadmonLayoutEntityWithAttributes
 {
  public:
   inline virtual                              ~RadmonVSubPhysicsListWithLabel();

   inline const G4String &                      GetLabel(void) const;
   inline virtual void                          SetPhysicsListAttribute(const G4String & attributeName, const G4String & value);

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const = 0;

  protected:
   inline                                       RadmonVSubPhysicsListWithLabel(const G4String & label);
   
  private:
  // Hidden constructors and operators
                                                RadmonVSubPhysicsListWithLabel();
                                                RadmonVSubPhysicsListWithLabel(const RadmonVSubPhysicsListWithLabel & copy);
   RadmonVSubPhysicsListWithLabel &             operator=(const RadmonVSubPhysicsListWithLabel & copy);

  // Private attributes
   G4String                                     physiscListLabel;
 };
 
 #include "RadmonVSubPhysicsListWithLabel.icc"
#endif /* RADMONVSUBPHYSICSLISTWITHLABEL_HH */

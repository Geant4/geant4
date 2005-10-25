//
// File name:     RadmonVGeneratorWithLabel.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGeneratorWithLabel.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Primary generators with a label abstract class
//

#ifndef   RADMONVGENERATORWITHLABEL_HH
 #define  RADMONVGENERATORWITHLABEL_HH

 // Include files
 #include "RadmonVGenerator.hh"
 #include "RadmonLayoutEntityWithAttributes.hh"

 class RadmonVGeneratorWithLabel : public RadmonVGenerator, public RadmonLayoutEntityWithAttributes
 {
  public:
   inline virtual                              ~RadmonVGeneratorWithLabel();

   inline const G4String &                      GetLabel(void) const;

   inline virtual void                          SetGeneratorAttribute(const G4String & attribute, const G4String & value);

   virtual RadmonVGeneratorWithLabel *          New(void) const = 0;

  protected:
   inline                                       RadmonVGeneratorWithLabel(const G4String & label);

  // Hidden constructors and operators
  private:
                                                RadmonVGeneratorWithLabel(const RadmonVGeneratorWithLabel & copy);
   RadmonVGeneratorWithLabel &                   operator=(const RadmonVGeneratorWithLabel & copy);

  // Private attributes
   G4String                                     generatorLabel;
 };
 
 #include "RadmonVGeneratorWithLabel.icc"   
#endif /* RADMONVGENERATORWITHLABEL_HH */

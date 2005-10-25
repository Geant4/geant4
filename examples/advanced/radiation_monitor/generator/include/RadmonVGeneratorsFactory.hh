//
// File name:     RadmonVGeneratorsFactory.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGeneratorsFactory.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract factory for a primary generators
//

#ifndef   RADMONVGENERATORSFACTORY_HH
 #define  RADMONVGENERATORSFACTORY_HH

 // Forward declarations
 class RadmonVGenerator;
 class G4String;
 
 class RadmonVGeneratorsFactory
 {
  public:
   inline                                     RadmonVGeneratorsFactory();
   inline virtual                            ~RadmonVGeneratorsFactory();

   virtual RadmonVGenerator *                 GetGenerator(const G4String & generatorType) = 0;

  // Hidden constructors and operators
  private:
                                              RadmonVGeneratorsFactory(const RadmonVGeneratorsFactory & copy);
   RadmonVGeneratorsFactory &                 operator=(const RadmonVGeneratorsFactory & copy);
 };
 
 // Inline implementations
 #include "RadmonVGeneratorsFactory.icc"
#endif /* RADMONVGENERATORSFACTORY_HH */

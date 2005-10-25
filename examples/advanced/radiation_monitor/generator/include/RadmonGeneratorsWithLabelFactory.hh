//
// File name:     RadmonGeneratorsWithLabelFactory.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelFactory.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Factory of primary generators with a label
//

#ifndef   RADMONGENERATORSWITHLABELFACTORY_HH
 #define  RADMONGENERATORSWITHLABELFACTORY_HH
 
 // Include files
 #include "RadmonVGeneratorsFactory.hh"
 #include <list>
 
 // Forward declarations
 class RadmonVGeneratorWithLabel;

 class RadmonGeneratorsWithLabelFactory : public RadmonVGeneratorsFactory
 {
  public:
                                                RadmonGeneratorsWithLabelFactory();
   virtual                                     ~RadmonGeneratorsWithLabelFactory();

   virtual RadmonVGenerator *                   GetGenerator(const G4String & generatorType);

   void                                         AppendGenerator(RadmonVGeneratorWithLabel * generator);

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorsWithLabelFactory(const RadmonGeneratorsWithLabelFactory & copy);
   RadmonGeneratorsWithLabelFactory &           operator=(const RadmonGeneratorsWithLabelFactory & copy);

  // Private data types
   typedef std::list<RadmonVGeneratorWithLabel *> GeneratorsList;

  // Private attributes
   GeneratorsList                                generatorsList;
 };
#endif /* RADMONGENERATORSWITHLABELFACTORY_HH */

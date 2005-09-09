//
// File name:     RadmonDetectorLabelledEntitiesConstructorsFactory.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLabelledEntitiesConstructorsFactory.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Concrete factory that constructs 
//                RadmonVDetectorLabelledEntityConstructor objects
//

#ifndef   RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSFACTORY_HH
 #define  RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSFACTORY_HH

 // Include files
 #include "RadmonVDetectorEntitiesConstructorsFactory.hh"
 #include <list>

 // Forward declarations
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonVDetectorEntityConstructor;

 class RadmonDetectorLabelledEntitiesConstructorsFactory : public RadmonVDetectorEntitiesConstructorsFactory
 {
  public:
   inline                                       RadmonDetectorLabelledEntitiesConstructorsFactory();
   virtual                                     ~RadmonDetectorLabelledEntitiesConstructorsFactory();

   virtual RadmonVDetectorEntityConstructor *   GetEntityConstructor(const G4String & entityName);

   void                                         AppendLabelledEntityConstructor(RadmonVDetectorLabelledEntityConstructor * constructor);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorLabelledEntitiesConstructorsFactory(const RadmonDetectorLabelledEntitiesConstructorsFactory & copy);
   RadmonDetectorLabelledEntitiesConstructorsFactory & operator=(const RadmonDetectorLabelledEntitiesConstructorsFactory & copy);

  // Private attributes
   typedef std::list<RadmonVDetectorLabelledEntityConstructor *> EntityConstructorsList;
   EntityConstructorsList                       entityConstructorsList;
 };
 
 // Inline implementations
 #include "RadmonDetectorLabelledEntitiesConstructorsFactory.icc"
#endif /* RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSFACTORY_HH */

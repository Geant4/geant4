//
// File name:     RadmonVDetectorEntitiesConstructorsFactory.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorEntitiesConstructorsFactory.hh,v 1.2 2005-09-19 19:42:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a factory of detector-entity constructor
//

#ifndef   RADMONVDETECTORENTITIESCONSTRUCTORSFACTORY_HH
 #define  RADMONVDETECTORENTITIESCONSTRUCTORSFACTORY_HH

 // Forward declaration
 class RadmonVDetectorEntityConstructor;
 class G4String;

 class RadmonVDetectorEntitiesConstructorsFactory
 {
  public:
   inline virtual                              ~RadmonVDetectorEntitiesConstructorsFactory();

   virtual RadmonVDetectorEntityConstructor *   CreateEntityConstructor(const G4String & entityName) = 0;

  protected:
   inline                                       RadmonVDetectorEntitiesConstructorsFactory();

  private:
  // Hidden constructors and operators
                                                RadmonVDetectorEntitiesConstructorsFactory(const RadmonVDetectorEntitiesConstructorsFactory & copy);
   RadmonVDetectorEntitiesConstructorsFactory & operator=(const RadmonVDetectorEntitiesConstructorsFactory & copy);
 };
 
 // Inline implementations
 #include "RadmonVDetectorEntitiesConstructorsFactory.icc"
#endif /* RADMONVDETECTORENTITIESCONSTRUCTORSFACTORY_HH */

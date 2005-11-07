//
// File name:     RadmonVSubPhysicsListFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVSubPhysicsListFactory.hh,v 1.1 2005-11-07 17:52:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a factory of detector-entity constructor
//

#ifndef   RADMONVSUBPHYSICSLISTFACTORY_HH
 #define  RADMONVSUBPHYSICSLISTFACTORY_HH

 // Forward declaration
 class RadmonVSubPhysicsList;
 class G4String;

 class RadmonVSubPhysicsListFactory
 {
  public:
   inline virtual                              ~RadmonVSubPhysicsListFactory();

   virtual RadmonVSubPhysicsList *              CreateSubPhysicsList(const G4String & subPhysicsListName) = 0;

  protected:
   inline                                       RadmonVSubPhysicsListFactory();

  private:
  // Hidden constructors and operators
                                                RadmonVSubPhysicsListFactory(const RadmonVSubPhysicsListFactory & copy);
   RadmonVSubPhysicsListFactory &               operator=(const RadmonVSubPhysicsListFactory & copy);
 };
 
 // Inline implementations
 #include "RadmonVSubPhysicsListFactory.icc"
#endif /* RADMONVSUBPHYSICSLISTFACTORY_HH */

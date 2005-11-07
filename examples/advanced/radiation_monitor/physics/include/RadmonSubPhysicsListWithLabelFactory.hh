//
// File name:     RadmonSubPhysicsListWithLabelFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelFactory.hh,v 1.1 2005-11-07 17:52:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Concrete factory that constructs 
//                RadmonSubPhysiscsListWithLabel objects
//

#ifndef   RADMONSUBPHYSICSLISTWITHLABELFACTORY_HH
 #define  RADMONSUBPHYSICSLISTWITHLABELFACTORY_HH

 // Include files
 #include "RadmonVSubPhysicsListFactory.hh"
 #include <list>

 // Forward declarations
 class RadmonVSubPhysicsList;
 class RadmonVSubPhysicsListWithLabel;

 class RadmonSubPhysicsListWithLabelFactory : public RadmonVSubPhysicsListFactory
 {
  public:
   inline                                       RadmonSubPhysicsListWithLabelFactory();
   virtual                                     ~RadmonSubPhysicsListWithLabelFactory();

   virtual RadmonVSubPhysicsList *              CreateSubPhysicsList(const G4String & subPhysicsListName);

   void                                         AppendSubPhysicsListWithLabel(RadmonVSubPhysicsListWithLabel * physicsList);

  private:
  // Hidden constructors and operators
                                                RadmonSubPhysicsListWithLabelFactory(const RadmonSubPhysicsListWithLabelFactory & copy);
   RadmonSubPhysicsListWithLabelFactory &       operator=(const RadmonSubPhysicsListWithLabelFactory & copy);

  // Private attributes
   typedef std::list<RadmonVSubPhysicsListWithLabel *> SubPhysicsLists;
   SubPhysicsLists                              subPhysicsLists;
 };
 
 // Inline implementations
 #include "RadmonSubPhysicsListWithLabelFactory.icc"
#endif /* RADMONSUBPHYSICSLISTWITHLABELFACTORY_HH */

//
// File name:     RadmonPhysicsInfoList.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfoList.hh,v 1.1 2005-11-10 08:14:10 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Provides informations about a process
//

#ifndef   RADMONPHYSICSINFOLIST_HH
 #define  RADMONPHYSICSINFOLIST_HH
 
 // Include files
 #include "RadmonPhysicsInfo.hh"
 #include "globals.hh"
 #include <vector>
 
 class RadmonPhysicsInfoList
 {
  public:
   inline                                       RadmonPhysicsInfoList();
                                                RadmonPhysicsInfoList(const RadmonPhysicsInfoList & copy);
   inline                                      ~RadmonPhysicsInfoList();
   
   RadmonPhysicsInfoList &                      operator=(const RadmonPhysicsInfoList & copy);
   
   G4bool                                       CollidesWith(const RadmonPhysicsInfoList & other) const;
   
   void                                         InsertPhysicsInfo(const RadmonPhysicsInfo & info);
   
   inline G4int                                 GetNPhysicsInfos(void) const;
   inline const RadmonPhysicsInfo &             GetPhysicsInfo(G4int index) const;
   
  private:
  // Private data types
   typedef std::vector<RadmonPhysicsInfo>       InfoVector;

  // Private attributes
   InfoVector                                   infoVector;
 };
 
 std::ostream &                                 operator<<(std::ostream & out, const RadmonPhysicsInfoList & infoList);
 
 // Inline implementations
 #include "RadmonPhysicsInfoList.icc"
#endif /* RADMONPHYSICSINFOLIST_HH */

//
// File name:     RadmonVSubPhysicsList.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVSubPhysicsList.hh,v 1.2 2005-11-10 08:14:10 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a physics list piece
//

#ifndef   RADMONVSUBPHYSICSLIST_HH
 #define  RADMONVSUBPHYSICSLIST_HH
 
 // Forward declarations
 class G4LogicalVolume;
 class RadmonPhysicsInfoList;
 class G4String;
 
 class RadmonVSubPhysicsList
 {
  public:
   inline virtual                              ~RadmonVSubPhysicsList();
    
   virtual void                                 SetPhysicsListAttribute(const G4String & attributeName, const G4String &value) = 0;

   virtual void                                 ConstructParticle(void) = 0;
   virtual void                                 ConstructProcess(void) = 0;
   virtual void                                 SetCuts(void) = 0;
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const = 0 ;

  protected:
   inline                                       RadmonVSubPhysicsList();

  private:
  // Hidden constructors and operators
                                                RadmonVSubPhysicsList(const RadmonVSubPhysicsList & copy);
   RadmonVSubPhysicsList &                      operator=(const RadmonVSubPhysicsList & copy);
 };
 
 // Inline implementations
 #include "RadmonVSubPhysicsList.icc"
#endif /* RADMONVSUBPHYSICSLIST_HH */

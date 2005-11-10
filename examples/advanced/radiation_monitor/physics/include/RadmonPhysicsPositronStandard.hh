//
// File name:     RadmonPhysicsPositronStandard.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsPositronStandard.hh,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Standard processes for electrons
//

#ifndef   RADMONPHYSICSPOSITRONSTANDARD_HH
 #define  RADMONPHYSICSPOSITRONSTANDARD_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsPositronStandard : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsPositronStandard();
   inline virtual                              ~RadmonPhysicsPositronStandard();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsPositronStandard(const RadmonPhysicsPositronStandard & copy);
   RadmonPhysicsPositronStandard &              operator=(const RadmonPhysicsPositronStandard & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsPositronStandard.icc"
#endif /* RADMONPHYSICSPOSITRONSTANDARD_HH */

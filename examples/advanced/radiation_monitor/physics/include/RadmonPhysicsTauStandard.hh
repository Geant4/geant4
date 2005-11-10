//
// File name:     RadmonPhysicsTauStandard.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsTauStandard.hh,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   EEDL processes for electrons
//

#ifndef   RADMONPHYSICSTAUSTANDARD_HH
 #define  RADMONPHYSICSTAUSTANDARD_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsTauStandard : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsTauStandard();
   inline virtual                              ~RadmonPhysicsTauStandard();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsTauStandard(const RadmonPhysicsTauStandard & copy);
   RadmonPhysicsTauStandard &                   operator=(const RadmonPhysicsTauStandard & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsTauStandard.icc"
#endif /* RADMONPHYSICSTAUSTANDARD_HH */

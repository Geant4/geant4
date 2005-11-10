//
// File name:     RadmonPhysicsPhotonStandard.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsPhotonStandard.hh,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   EEDL processes for electrons
//

#ifndef   RADMONPHYSICSPHOTONSTANDARD_HH
 #define  RADMONPHYSICSPHOTONSTANDARD_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsPhotonStandard : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsPhotonStandard();
   inline virtual                              ~RadmonPhysicsPhotonStandard();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsPhotonStandard(const RadmonPhysicsPhotonStandard & copy);
   RadmonPhysicsPhotonStandard &                operator=(const RadmonPhysicsPhotonStandard & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsPhotonStandard.icc"
#endif /* RADMONPHYSICSPHOTONSTANDARD_HH */

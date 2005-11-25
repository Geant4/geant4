//
// File name:     RadmonPhysicsICRUIonization.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsICRUIonization.hh,v 1.1 2005-11-25 01:52:01 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   ICRU ionization processe
//

#ifndef   RADMONPHYSICSICRUIONIZATION_HH
 #define  RADMONPHYSICSICRUIONIZATION_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsICRUIonization : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsICRUIonization();
   inline virtual                              ~RadmonPhysicsICRUIonization();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsICRUIonization(const RadmonPhysicsICRUIonization & copy);
   RadmonPhysicsICRUIonization &                         operator=(const RadmonPhysicsICRUIonization & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsICRUIonization.icc"
#endif /* RADMONPHYSICSICRUIONIZATION_HH */

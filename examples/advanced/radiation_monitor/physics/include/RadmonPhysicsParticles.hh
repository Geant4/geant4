//
// File name:     RadmonPhysicsParticles.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsParticles.hh,v 1.1 2006-01-06 12:52:32 guatelli Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Particles production cuts
//

#ifndef   RADMONPHYSICSPARTICLES_HH
 #define  RADMONPHYSICSPARTICLES_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsParticles : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsParticles();
   inline virtual                              ~RadmonPhysicsParticles();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  private:
  // Hidden constructors and operators
                                                RadmonPhysicsParticles(const RadmonPhysicsParticles & copy);
   RadmonPhysicsParticles &                     operator=(const RadmonPhysicsParticles & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsParticles.icc"
#endif /* RADMONPHYSICSPARTICLES_HH */

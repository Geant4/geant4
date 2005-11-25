//
// File name:     RadmonPhysicsNeutronBertini.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsNeutronBertini.hh,v 1.1 2005-11-25 01:52:01 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Neutron processes based on bertini cascade model
//

#ifndef   RADMONPHYSICSNEUTRONBERTINI_HH
 #define  RADMONPHYSICSNEUTRONBERTINI_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsNeutronBertini : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsNeutronBertini();
   inline virtual                              ~RadmonPhysicsNeutronBertini();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsNeutronBertini(const RadmonPhysicsNeutronBertini & copy);
   RadmonPhysicsNeutronBertini &                operator=(const RadmonPhysicsNeutronBertini & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsNeutronBertini.icc"
#endif /* RADMONPHYSICSNEUTRONBERTINI_HH */

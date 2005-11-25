//
// File name:     RadmonPhysicsNeutronBinary.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsNeutronBinary.hh,v 1.2 2005-11-25 01:52:01 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Neutron processes based on binary cascade model
//

#ifndef   RADMONPHYSICSNEUTRONBINARY_HH
 #define  RADMONPHYSICSNEUTRONBINARY_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsNeutronBinary : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsNeutronBinary();
   inline virtual                              ~RadmonPhysicsNeutronBinary();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsNeutronBinary(const RadmonPhysicsNeutronBinary & copy);
   RadmonPhysicsNeutronBinary &                 operator=(const RadmonPhysicsNeutronBinary & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsNeutronBinary.icc"
#endif /* RADMONPHYSICSNEUTRONBINARY_HH */

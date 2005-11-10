//
// File name:     RadmonPhysicsNuclear.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsNuclear.hh,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   EEDL processes for electrons
//

#ifndef   RADMONPHYSICSNUCLEAR_HH
 #define  RADMONPHYSICSNUCLEAR_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsNuclear : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsNuclear();
   inline virtual                              ~RadmonPhysicsNuclear();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsNuclear(const RadmonPhysicsNuclear & copy);
   RadmonPhysicsNuclear &                       operator=(const RadmonPhysicsNuclear & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsNuclear.icc"
#endif /* RADMONPHYSICSNUCLEAR_HH */

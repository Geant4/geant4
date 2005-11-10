//
// File name:     RadmonPhysicsElectronStandard.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsElectronStandard.hh,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Standard processes for electrons
//

#ifndef   RADMONPHYSICSELECTRONSTANDARD_HH
 #define  RADMONPHYSICSELECTRONSTANDARD_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsElectronStandard : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsElectronStandard();
   inline virtual                              ~RadmonPhysicsElectronStandard();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsElectronStandard(const RadmonPhysicsElectronStandard & copy);
   RadmonPhysicsElectronStandard &              operator=(const RadmonPhysicsElectronStandard & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsElectronStandard.icc"
#endif /* RADMONPHYSICSELECTRONSTANDARD_HH */

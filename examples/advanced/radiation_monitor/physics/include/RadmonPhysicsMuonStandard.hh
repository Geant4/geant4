//
// File name:     RadmonPhysicsMuonStandard.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsMuonStandard.hh,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   EEDL processes for electrons
//

#ifndef   RADMONPHYSICSMUONSTANDARD_HH
 #define  RADMONPHYSICSMUONSTANDARD_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsMuonStandard : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsMuonStandard();
   inline virtual                              ~RadmonPhysicsMuonStandard();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsMuonStandard(const RadmonPhysicsMuonStandard & copy);
   RadmonPhysicsMuonStandard &                  operator=(const RadmonPhysicsMuonStandard & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsMuonStandard.icc"
#endif /* RADMONPHYSICSMUONSTANDARD_HH */

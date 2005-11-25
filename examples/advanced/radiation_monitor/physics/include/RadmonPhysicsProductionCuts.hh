//
// File name:     RadmonPhysicsProductionCuts.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsProductionCuts.hh,v 1.1 2005-11-25 11:52:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Particles production cuts
//

#ifndef   RADMONPHYSICSPRODUCTIONCUTS_HH
 #define  RADMONPHYSICSPRODUCTIONCUTS_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsProductionCuts : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsProductionCuts();
   inline virtual                              ~RadmonPhysicsProductionCuts();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  private:
   void                                         SetProductionCut(G4double cut, const G4String& particleName) const;

  // Hidden constructors and operators
                                                RadmonPhysicsProductionCuts(const RadmonPhysicsProductionCuts & copy);
   RadmonPhysicsProductionCuts &                         operator=(const RadmonPhysicsProductionCuts & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsProductionCuts.icc"
#endif /* RADMONPHYSICSPRODUCTIONCUTS_HH */

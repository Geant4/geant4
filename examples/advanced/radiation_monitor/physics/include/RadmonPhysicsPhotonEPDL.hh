//
// File name:     RadmonPhysicsPhotonEPDL.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsPhotonEPDL.hh,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   EEDL processes for electrons
//

#ifndef   RADMONPHYSICSPHOTONEPDL_HH
 #define  RADMONPHYSICSPHOTONEPDL_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsPhotonEPDL : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsPhotonEPDL();
   inline virtual                              ~RadmonPhysicsPhotonEPDL();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsPhotonEPDL(const RadmonPhysicsPhotonEPDL & copy);
   RadmonPhysicsPhotonEPDL &                    operator=(const RadmonPhysicsPhotonEPDL & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsPhotonEPDL.icc"
#endif /* RADMONPHYSICSPHOTONEPDL_HH */

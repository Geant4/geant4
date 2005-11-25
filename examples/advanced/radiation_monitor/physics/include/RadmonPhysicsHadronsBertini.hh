//
// File name:     RadmonPhysicsHadronsBertini.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsHadronsBertini.hh,v 1.1 2005-11-25 01:52:01 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Hadrons (except neutron) processes based on bertini cascade model
//

#ifndef   RADMONPHYSICSHADRONSBERTINI_HH
 #define  RADMONPHYSICSHADRONSBERTINI_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsHadronsBertini : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsHadronsBertini();
   inline virtual                              ~RadmonPhysicsHadronsBertini();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsHadronsBertini(const RadmonPhysicsHadronsBertini & copy);
   RadmonPhysicsHadronsBertini &                operator=(const RadmonPhysicsHadronsBertini & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsHadronsBertini.icc"
#endif /* RADMONPHYSICSHADRONSBERTINI_HH */

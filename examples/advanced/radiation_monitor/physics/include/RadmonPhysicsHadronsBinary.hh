//
// File name:     RadmonPhysicsHadronsBinary.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsHadronsBinary.hh,v 1.1 2005-11-25 01:52:01 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Hadrons (except neutron) processes based on binary cascade model
//

#ifndef   RADMONPHYSICSHADRONSBINARY_HH
 #define  RADMONPHYSICSHADRONSBINARY_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsHadronsBinary : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsHadronsBinary();
   inline virtual                              ~RadmonPhysicsHadronsBinary();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsHadronsBinary(const RadmonPhysicsHadronsBinary & copy);
   RadmonPhysicsHadronsBinary &                 operator=(const RadmonPhysicsHadronsBinary & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsHadronsBinary.icc"
#endif /* RADMONPHYSICSHADRONSBINARY_HH */

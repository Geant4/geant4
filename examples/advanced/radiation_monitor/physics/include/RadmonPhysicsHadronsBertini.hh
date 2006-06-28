//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonPhysicsHadronsBertini.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsHadronsBertini.hh,v 1.2 2006-06-28 13:54:33 gunter Exp $
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

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
// File name:     RadmonPhysicsProductionCuts.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsProductionCuts.hh,v 1.2 2006-06-28 13:55:27 gunter Exp $
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

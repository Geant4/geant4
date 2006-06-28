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
// File name:     RadmonPhysicsSteppingAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsSteppingAction.hh,v 1.2 2006-06-28 13:55:31 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class for the Stepping Action
//

#ifndef   RADMONPHYSICSSTEPPINGACTION_HH
 #define  RADMONPHYSICSSTEPPINGACTION_HH

 #include "RadmonSteppingActionObserver.hh"
 
 class RadmonPhysicsSteppingAction : public RadmonSteppingActionObserver
 {
  public:
   inline                                       RadmonPhysicsSteppingAction();
   inline virtual                              ~RadmonPhysicsSteppingAction();
   
   virtual void                                 OnUserStepping(const G4Step * step);

  private:
                                                RadmonPhysicsSteppingAction(const RadmonPhysicsSteppingAction & copy);
   RadmonPhysicsSteppingAction &                operator=(const RadmonPhysicsSteppingAction & copy);
 };
 
 // Inline implementations
 #include "RadmonPhysicsSteppingAction.icc"
#endif /* RADMONPHYSICSSTEPPINGACTION_HH */

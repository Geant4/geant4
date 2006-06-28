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
// File name:     RadmonSteppingActionObserver.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSteppingActionObserver.hh,v 1.2 2006-06-28 13:57:35 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class for the Stepping Action
//

#ifndef   RADMONSTEPPINGACTIONOBSERVER_HH
 #define  RADMONSTEPPINGACTIONOBSERVER_HH

 class G4Step;
 
 class RadmonSteppingActionObserver
 {
  public:
   inline                                       RadmonSteppingActionObserver();
   inline virtual                              ~RadmonSteppingActionObserver();
   
   virtual void                                 OnUserStepping(const G4Step * step) = 0;

  private:
                                                RadmonSteppingActionObserver(const RadmonSteppingActionObserver & copy);
   RadmonSteppingActionObserver &               operator=(const RadmonSteppingActionObserver & copy);
 };
 
 // Inline implementations
 #include "RadmonSteppingActionObserver.icc"
#endif /* RADMONSTEPPINGACTIONOBSERVER_HH */

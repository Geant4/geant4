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
// File name:     RadmonSteppingAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSteppingAction.hh,v 1.2 2006-06-28 13:57:31 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon stepping action class
//

#ifndef   RADMONSTEPPINGACTION_HH
 #define  RADMONSTEPPINGACTION_HH
 
 // Include files
 #include "G4UserSteppingAction.hh"
 #include "globals.hh"
 #include <set>
 
 // Forward declarations
 class RadmonSteppingActionObserver;
 
 class RadmonSteppingAction : public G4UserSteppingAction
 {
  public:
   inline virtual                              ~RadmonSteppingAction();
   
   inline static RadmonSteppingAction *         Instance(void);

   void                                         AttachObserver(RadmonSteppingActionObserver * observer);
   void                                         DetachObserver(RadmonSteppingActionObserver * observer);
   
   virtual void                                 UserSteppingAction(const G4Step * step);
   
  private:
                                                RadmonSteppingAction();
  // Hidden constructors and operators
                                                RadmonSteppingAction(const RadmonSteppingAction & copy);
   RadmonSteppingAction &                       operator=(const RadmonSteppingAction & copy);
   
  // Private attributes
   typedef std::set<RadmonSteppingActionObserver *> ObserversSet;
   ObserversSet                                 observersSet;
   
   static RadmonSteppingAction *                instance;
 };
 
 // Inline implementations
 #include "RadmonSteppingAction.icc"
#endif /* RADMONSTEPPINGACTION_HH */

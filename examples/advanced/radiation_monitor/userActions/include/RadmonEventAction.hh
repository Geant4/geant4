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
// File name:     RadmonEventAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonEventAction.hh,v 1.2 2006-06-28 13:57:15 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon event action class
//

#ifndef   RADMONEVENTACTION_HH
 #define  RADMONEVENTACTION_HH
 
 // Include files
 #include "G4UserEventAction.hh"
 #include "globals.hh"
 #include <set>
 
 // Forward declarations
 class RadmonEventActionObserver;
 
 class RadmonEventAction : public G4UserEventAction
 {
  public:
   inline virtual                              ~RadmonEventAction();
   
   inline static RadmonEventAction *            Instance(void);

   void                                         AttachObserver(RadmonEventActionObserver * observer);
   void                                         DetachObserver(RadmonEventActionObserver * observer);
   
   virtual void                                 BeginOfEventAction(const G4Event * event);
   virtual void                                 EndOfEventAction(const G4Event * event);
   
  private:
                                                RadmonEventAction();
  // Hidden constructors and operators
                                                RadmonEventAction(const RadmonEventAction & copy);
   RadmonEventAction &                          operator=(const RadmonEventAction & copy);
   
  // Private attributes
   typedef std::set<RadmonEventActionObserver *> ObserversSet;
   ObserversSet                                 observersSet;
   
   static RadmonEventAction *                   instance;
 };
 
 // Inline implementations
 #include "RadmonEventAction.icc"
#endif /* RADMONEVENTACTION_HH */

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
// File name:     RadmonRunAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonRunAction.hh,v 1.2 2006-06-28 13:57:23 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon run action class
//

#ifndef   RADMONRUNACTION_HH
 #define  RADMONRUNACTION_HH
 
 // Include files
 #include "G4UserRunAction.hh"
 #include "globals.hh"
 #include <set>
 
 // Forward declarations
 class RadmonRunActionObserver;
 
 class RadmonRunAction : public G4UserRunAction
 {
  public:
   inline virtual                              ~RadmonRunAction();
   
   inline static RadmonRunAction *              Instance(void);

   void                                         AttachObserver(RadmonRunActionObserver * observer);
   void                                         DetachObserver(RadmonRunActionObserver * observer);
   
   virtual void                                 BeginOfRunAction(const G4Run * run);
   virtual void                                 EndOfRunAction(const G4Run * run);
   
  private:
                                                RadmonRunAction();
  // Hidden constructors and operators
                                                RadmonRunAction(const RadmonRunAction & copy);
   RadmonRunAction &                            operator=(const RadmonRunAction & copy);
   
  // Private attributes
   typedef std::set<RadmonRunActionObserver *>  ObserversSet;
   ObserversSet                                 observersSet;
   
   static RadmonRunAction *                     instance;
 };
 
 // Inline implementations
 #include "RadmonRunAction.icc"
#endif /* RADMONRUNACTION_HH */

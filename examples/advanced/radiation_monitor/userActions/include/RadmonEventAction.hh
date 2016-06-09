//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonEventAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonEventAction.hh,v 1.3 2006/06/29 16:20:36 gunter Exp $
// Tag:           $Name: geant4-09-01 $
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

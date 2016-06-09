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
// File name:     RadmonSteppingAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSteppingAction.hh,v 1.3 2006/06/29 16:21:48 gunter Exp $
// Tag:           $Name: geant4-08-02 $
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

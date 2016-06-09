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
// File name:     RadmonRunAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonRunAction.hh,v 1.3 2006/06/29 16:21:00 gunter Exp $
// Tag:           $Name: geant4-09-01 $
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

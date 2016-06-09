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
// File name:     RadmonHitsManager.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHitsManager.hh,v 1.3 2006/06/29 16:20:18 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//
// Description:   Global features for hits
//

#ifndef   RADMONHITSMANAGER_HH
 #define  RADMONHITSMANAGER_HH
 
 // Include files
 #include "G4Allocator.hh"
 #include "globals.hh"
 #include <map>
 
 class RadmonHit;
 
 class RadmonHitsManager
 {
  public:
   typedef G4Allocator<RadmonHit>                Allocator;

   inline static RadmonHitsManager *             Instance();

   inline Allocator &                            HitsAllocator();
   G4int                                         ReserveIdByLabel(const G4String & label);
   
  private:
                                                RadmonHitsManager();
   inline                                      ~RadmonHitsManager();

  // Hidden constructors and operators
                                                RadmonHitsManager(const RadmonHitsManager & copy);
   RadmonHitsManager &                          operator=(const RadmonHitsManager & copy);
   
   typedef std::map<G4String, G4int>            IdLabels;
   
   static RadmonHitsManager *                   instance;

   IdLabels                                     idLabels;
   G4int                                        idCounter;
   Allocator *                                  allocator;
 };
 
 // Inline implementations
 #include "RadmonHitsManager.icc"
#endif /* RADMONHITSMANAGER_HH */

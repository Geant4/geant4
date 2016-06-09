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
// File name:     RadmonHit.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHit.hh,v 1.3 2006/06/29 16:20:14 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//
// Description:   Radmon hit
//

#ifndef   RADMONHIT_HH
 #define  RADMONHIT_HH
 
 // Include files
 #include "G4VHit.hh"
 #include "G4THitsCollection.hh"
 #include "RadmonHitsManager.hh"
 #include "globals.hh"
 #include <map>
 
 class RadmonHit : public G4VHit
 {
  public:
   inline                                       RadmonHit();
   inline virtual                              ~RadmonHit();

   void *                                       operator new(size_t size);
   void                                         operator delete(void * hit);
   
   G4double                                     RetrieveById(G4int id) const;
   void                                         StoreById(G4int id, G4double value);

   inline G4double                              RetrieveByLabel(const G4String & label) const;
   inline void                                  StoreByLabel(const G4String & label, G4double value);
   
  private:
  // Hidden constructors and operators
                                                RadmonHit(const RadmonHit & copy);
   RadmonHit &                                  operator=(const RadmonHit & copy);
   
   typedef std::map<G4int, G4double>            HitData;
   
   HitData                                      hitData;
 };
 
 typedef G4THitsCollection<RadmonHit>           RadmonHitsCollection;

 // Inline implementations
 #include "RadmonHit.icc"
#endif /* RADMONHIT_HH */

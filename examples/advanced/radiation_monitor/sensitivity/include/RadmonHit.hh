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
// File name:     RadmonHit.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHit.hh,v 1.2 2006-06-28 13:56:53 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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

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
// File name:     RadmonSensitiveDetector.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSensitiveDetector.hh,v 1.3 2006-06-28 13:57:01 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon sensitive detector class
//

#ifndef   RADMONSENSITIVEDETECTOR_HH
 #define  RADMONSENSITIVEDETECTOR_HH
 
 // Include files
 #include "G4VSensitiveDetector.hh"
 #include "RadmonHit.hh"
 #include "globals.hh"
 #include <set>
 
 class RadmonSensitiveDetectorDataStorer;
 
 class RadmonSensitiveDetector : public G4VSensitiveDetector
 {
  public:
                                                RadmonSensitiveDetector(const G4String & name);
   virtual                                     ~RadmonSensitiveDetector();
   
   void                                         ClearDataStorersList();
   void                                         AttachDataStorer(RadmonSensitiveDetectorDataStorer * observer);
   
   virtual void                                 Initialize(G4HCofThisEvent * hitsCollections);
   inline virtual void                          EndOfEvent(G4HCofThisEvent * hitsCollections);
   virtual G4bool                               ProcessHits(G4Step * step, G4TouchableHistory * touchableHistory);
   
   RadmonHitsCollection *                       GetDetectorCollection(void) const;
   
  private:
  // Hidden constructors and operators
                                                RadmonSensitiveDetector();
                                                RadmonSensitiveDetector(const RadmonSensitiveDetector & copy);
   RadmonSensitiveDetector &                    operator=(const RadmonSensitiveDetector & copy);
   
  // Private attributes
   typedef std::set<RadmonSensitiveDetectorDataStorer *> DataStorersSet;
   DataStorersSet                               dataStorersSet;

   RadmonHitsCollection *                       hitsCollection;
   G4String                                     collName;
 };
 
 // Inline implementations
 #include "RadmonSensitiveDetector.icc"
#endif /* RADMONSENSITIVEDETECTOR_HH */

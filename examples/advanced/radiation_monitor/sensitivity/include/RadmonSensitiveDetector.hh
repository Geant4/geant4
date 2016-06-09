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
// File name:     RadmonSensitiveDetector.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSensitiveDetector.hh,v 1.4 2006/06/29 16:20:22 gunter Exp $
// Tag:           $Name: geant4-09-00 $
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

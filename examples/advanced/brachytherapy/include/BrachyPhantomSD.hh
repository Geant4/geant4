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
// $Id: BrachyPhantomSD.hh,v 1.3 2002-11-18 15:18:36 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    ********************************
//    *                              *  
//    *     BrachyPhantomSD.hh      *
//    *                              *
//    ********************************

#ifndef BrachyPhantomSD_h
#define BrachyPhantomSD_h 1

#include "G4VSensitiveDetector.hh"
#include "BrachyPhantomHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class BrachyPhantomSD : public G4VSensitiveDetector
{
 public:
      BrachyPhantomSD(G4String name, G4int NumVoxelX, G4int NumVoxelZ);
      ~BrachyPhantomSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();
 private:
      BrachyPhantomHitsCollection *m_pPhantomHitsCollection;
	
      const G4int m_NumVoxelX;
      const G4int m_NumVoxelZ;
      G4int *m_pVoxelID;
};
#endif



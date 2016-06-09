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
// $Id: ExN07CalorimeterSD.hh,v 1.1 2003/03/10 01:43:35 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#ifndef ExN07CalorimeterSD_h
#define ExN07CalorimeterSD_h 1

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;
#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "ExN07CalorHit.hh"

class ExN07CalorimeterSD : public G4VSensitiveDetector
{
  public:
      ExN07CalorimeterSD(G4String);
      virtual ~ExN07CalorimeterSD();

      virtual void Initialize(G4HCofThisEvent*);
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();
      virtual void DrawAll();
      virtual void PrintAll();

  private:
      static G4int numberOfLayers;
      ExN07CalorHitsCollection*  AbsCollection;      
      ExN07CalorHitsCollection*  GapCollection;      
      G4int AbsCollID;
      G4int GapCollID;

  public:
      static void SetNumberOfLayers(G4int nl)
      { numberOfLayers = nl; }
};

#endif


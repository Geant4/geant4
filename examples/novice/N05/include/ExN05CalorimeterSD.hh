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
// $Id: ExN05CalorimeterSD.hh,v 1.4 2002-01-09 17:24:18 ranjard Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN05CalorimeterSD_h
#define ExN05CalorimeterSD_h 1

#include "ExN05CalorimeterHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class ExN05CalorimeterSD : public G4VSensitiveDetector
{

  public:
      ExN05CalorimeterSD(G4String name, G4int nCells, G4String colName);
      ~ExN05CalorimeterSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      ExN05CalorimeterHitsCollection *CalCollection;

      int* CellID;
      int numberOfCells;
      int HCID;
};




#endif


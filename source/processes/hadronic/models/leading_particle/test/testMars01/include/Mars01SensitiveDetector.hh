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
// $Id: Mars01SensitiveDetector.hh,v 1.1 2001-12-13 14:58:43 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Mars01SensitiveDetector_h
#define Mars01SensitiveDetector_h 1

#include "Mars01Hit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class Mars01SensitiveDetector : public G4VSensitiveDetector
{

  public:
      Mars01SensitiveDetector(G4String name);
      ~Mars01SensitiveDetector();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      Mars01HitsCollection *HitCollection;
      int HCID;
};




#endif


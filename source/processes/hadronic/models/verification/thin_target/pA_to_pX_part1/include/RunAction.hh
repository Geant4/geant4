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
// $Id: RunAction.hh,v 1.2 2003-05-29 15:32:39 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "EnergyHists.hh"

class G4Run;
class G4Material;
class TargetConstruction;
class EnergyAngleCrossSection;


class RunAction : public G4UserRunAction
{
  public:
    RunAction(TargetConstruction* targ);
   ~RunAction();

    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

    void FillHists(G4double charge, G4double ke, G4double costheta);

  private:

    G4double loBinEdge;
    G4double hiBinEdge;
    G4double binSize;

    EnergyHists ehists;

    TargetConstruction* target;

    EnergyAngleCrossSection* theData; 

    G4std::vector<G4std::pair<G4double,G4double> > binEdges;

    G4Material* material;
    G4double twopi; 
};

#endif



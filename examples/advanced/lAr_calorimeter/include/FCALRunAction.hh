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
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALRunAction.hh,v 1.5 2003-02-14 15:54:43 pmendez Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



//#ifndef  G4ANALYSIS_USE

#ifndef FCALRunAction_h
#define FCALRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

//namespace AIDA {
//  class ITree;
//  class IHistogram1D;
//  class IHistogram2D;
//  class ITuple;
//}

class FCALRunAction : public G4UserRunAction
{
  public:
    FCALRunAction();
   ~FCALRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
//
//    AIDA::IHistogram1D* GetHisto(G4int id) {return histo[id];}
//    AIDA::ITuple* GetTuple(G4int id){return ntuple[id];}
//
//  private:  
//    void bookHisto();
//    void cleanHisto();
//
//private:      
//    AIDA::ITree* tree;
//    AIDA::IHistogram1D* histo[4];
//    AIDA::ITuple* ntuple[3];
//
};

//#endif

#endif

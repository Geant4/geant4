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
// $Id: Em4RunAction.hh,v 1.7 2001-10-17 14:04:14 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em4RunAction_h
#define Em4RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class Em4RunActionMessenger;

#ifndef G4NOHIST
  class HepTupleManager;
  class HepHistogram;
#endif

class Em4RunAction : public G4UserRunAction
{
  public:
    Em4RunAction();
   ~Em4RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void  SetRndmFreq(G4int    val) {saveRndm = val;}
    G4int GetRndmFreq()             {return saveRndm;}
    
#ifndef G4NOHIST    
    HepHistogram* GetHisto(G4int id) {return histo[id];}
#endif
        
  private:
    void bookHisto();
    
  private:  
    Em4RunActionMessenger* runMessenger;
    G4int saveRndm;
    
#ifndef G4NOHIST         
    HepTupleManager* hbookManager;    
    HepHistogram* histo[1];
#endif         
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


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
// $Id: EventAction.hh,v 1.5 2009-08-29 08:48:30 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"

#include <vector>

class RunAction;
class PrimaryGeneratorAction;
class EventActionMessenger;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:  
    EventAction(DetectorConstruction*, RunAction*, PrimaryGeneratorAction*,
                HistoManager*);
   ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag    = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    
    void SumVisibleEnergy(G4int pixel, G4double de)
                               {visibleEnergy[pixel] += de;};
			         	    
    void SumTotalEnergy  (G4int pixel, G4double de)
                                 {totalEnergy[pixel] += de;};
				 
    void SumNbRadLength  (G4double dn)  {nbRadLen += dn;};
    
    void SetWriteFile(G4bool);    
    void WritePixels(const G4Event*);
    				         
  private:  
    DetectorConstruction*   detector;
    RunAction*              runAct;
    PrimaryGeneratorAction* primary;
          
    std::vector<G4double>   visibleEnergy;
    std::vector<G4double>     totalEnergy;
    G4double                nbRadLen;
    
    G4bool                trigger;
    G4double              Eseuil;   
    
    G4bool                writeFile;
                    
    G4String              drawFlag; 
    G4int                 printModulo;         
    EventActionMessenger* eventMessenger;
    HistoManager*         histoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    

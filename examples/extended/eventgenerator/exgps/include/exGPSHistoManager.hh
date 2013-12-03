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
/// \file eventgenerator/exgps/include/exGPSHistoManager.hh
/// \brief Definition of the exGPSHistoManager class
//
//
// $Id: exGPSHistoManager.hh 74272 2013-10-02 14:48:50Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef exGPSHistoManager_h
#define exGPSHistoManager_h 1

#include "globals.hh"
#include "g4root.hh"
////#include "g4xml.hh"
////#include "g4hbook.hh"
class exGPSHistoMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exGPSHistoManager
{
  public:

    exGPSHistoManager();
   ~exGPSHistoManager();

    void book();
    void save();
    
    void Fill(G4int PDGid, G4double e,
                                    G4double x, G4double y, G4double z,
                                    G4double t, G4double p, G4double w);

    void PrintStatistic();        

  private:

    G4String      fFileName[2];
    G4bool        fFactoryOn;
    G4double fMinpos,fMaxpos,fMineng,fMaxeng;
    G4AnaH1 *fEnerHisto;
    G4AnaH2 *fPosiXY, *fPosiZX, *fPosiYZ, *fAnglCTP, *fAnglTP;
    exGPSHistoMessenger* fMessenger;

  public:
    inline void SetEMin(G4double emin) {fMineng = emin;}
    inline void SetEMax(G4double emax) {fMaxeng = emax;}
    inline void SetPosMin(G4double posmin) {fMinpos = posmin;}
    inline void SetPosMax(G4double posmax) {fMaxpos = posmax;}
    inline void SetFileName(G4String name) {fFileName[0] = name;}

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


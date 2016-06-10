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
// $Id: exGPSEventAction.hh 76468 2013-11-11 10:27:19Z gcosmo $
//
/// \file eventgenerator/exgps/include/exGPSEventAction.hh
/// \brief Definition of the exGPSEventAction class
//

#ifndef exGPSEventAction_h
#define exGPSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class exGPSEventActionMessenger;
class exGPSHistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exGPSEventAction : public G4UserEventAction
{
  public:
    exGPSEventAction(exGPSHistoManager* histoManager) ;
    virtual ~exGPSEventAction();

    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {fDrawFlag = val;};
    void SetPrintModulo(G4int    val)  {fPrintModulo = val;};
    
  private:
    G4String                    fDrawFlag;
    G4int                       fPrintModulo;                         
    exGPSEventActionMessenger*  fEventMessenger;
    exGPSHistoManager* fexGPSHistoManager;
};

#endif

    

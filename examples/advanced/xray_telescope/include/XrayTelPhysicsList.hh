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
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelPhysicsList.hh                           *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
//
//
// **********************************************************************

#ifndef XrayTelPhysicsList_h
#define XrayTelPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class XrayTelPhysicsList: public G4VUserPhysicsList
{
public:
  XrayTelPhysicsList();
  ~XrayTelPhysicsList();

protected:
  // Construct particle and physics process
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();

public:
  // Set/Get cut values 
  void      SetCutForGamma(G4double);
  void      SetCutForElectron(G4double);
  G4double  GetCutForGamma() const;
  G4double  GetCutForElectron() const;
    
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructIons();
  void ConstructAllShortLiveds();

protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();

private:
  G4double cutForGamma;
  G4double cutForElectron; 
};

#endif








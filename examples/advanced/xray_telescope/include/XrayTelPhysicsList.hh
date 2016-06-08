// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
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
  void      SetCutForProton(G4double);           
  G4double  GetCutForGamma() const;
  G4double  GetCutForElectron() const;
  G4double  GetCutForProton() const;
    
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructAllShortLiveds();

protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();

private:
  G4double cutForGamma;
  G4double cutForElectron; 
  G4double cutForProton;
};

#endif








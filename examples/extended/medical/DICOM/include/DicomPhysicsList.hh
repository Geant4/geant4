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
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************//

#ifndef DicomPhysicsList_h
#define DicomPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;

class DicomPhysicsList: public G4VUserPhysicsList
{
public:
  DicomPhysicsList();
  ~DicomPhysicsList();

  // Set Cuts
  void SetGammaCut(G4double);
  void SetElectronCut(G4double);
  void SetPositronCut(G4double);

  void SetGammaLowLimit(G4double);
  void SetElectronLowLimit(G4double);
  void SetGELowLimit(G4double);
  void SetLowEnSecPhotCut(G4double);
  void SetLowEnSecElecCut(G4double);

protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();

  void SetCuts();

  // these methods Construct particles
  void ConstructBosons();
  void ConstructLeptons();

  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();

private:

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  G4LowEnergyIonisation*  loweIon;
  G4LowEnergyPhotoElectric* lowePhot;
  G4LowEnergyBremsstrahlung* loweBrem;

};

#endif




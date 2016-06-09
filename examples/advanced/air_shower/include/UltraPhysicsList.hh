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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//   **********************************************
//   *        UltraPhysicsList.hh
//   **********************************************
//
//    Ultra Physics List class; Standard and Low Energy EM processes are defined for
//    the relevant particles. Optical processes are declared.
//
#ifndef UltraPhysicsList_H
#define UltraPhysicsList_H 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class UltraPhysicsList : public G4VUserPhysicsList
{
  public:
    UltraPhysicsList();
    ~UltraPhysicsList();
 
  protected:
    // Construct particles and processes
    void ConstructParticle();
    void ConstructProcess();

    //
    void SetCuts();

  protected:
    // these methods Construct particles
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();

  protected:
    // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    void ConstructOp();

};

#endif


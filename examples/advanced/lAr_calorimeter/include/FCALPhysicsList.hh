// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALPhysicsList.hh,v 1.4 2002-10-02 12:10:27 araujo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// $Id: FCALPhysicsList.hh,v 1.4 2002-10-02 12:10:27 araujo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALPhysicsList_h
#define FCALPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALPhysicsList: public G4VUserPhysicsList
{
  public:
    FCALPhysicsList();
   ~FCALPhysicsList();

  protected:
    // Construct particle and physics
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

  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();

  private:
    G4double cutForGamma;
    G4double cutForElectron; 
    G4double cutForProton;
    G4double currentDefaultCut;
};

#endif


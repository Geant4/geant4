// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05PhysicsList.hh,v 1.1 1999-01-07 16:06:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#ifndef ExN05PhysicsList_h
#define ExN05PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class ExN05PhysicsList: public G4VUserPhysicsList
{
public:
  ExN05PhysicsList();
  virtual ~ExN05PhysicsList();
  
protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
  // 
  virtual void SetCuts(G4double aCut);
  
protected:
  // these methods Construct particles 
  virtual void ConstructBosons();
  virtual void ConstructLeptons();
  virtual void ConstructMesons();
  virtual void ConstructBarions();
  virtual void ConstructIons();
  
protected:
  // these methods Construct physics processes and register them
  void AddParameterisation();
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructHad();
  virtual void ConstructLeptHad();
  virtual void AddTransportation();
};

#endif




// -*- C++ -*-
// $Id: g4rwPhysicsList.hh,v 1.1 2004-12-06 11:01:14 radoone Exp $

#ifndef gogdmlPhysicsList_h
#define gogdmlPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class gogdmlPhysicsList: public G4VUserPhysicsList
{
  public:
    gogdmlPhysicsList();
    ~gogdmlPhysicsList();

  protected:
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

};

#endif








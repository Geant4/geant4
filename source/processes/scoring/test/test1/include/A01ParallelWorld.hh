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
// $Id: A01ParallelWorld.hh,v 1.1 2006-07-14 14:43:01 asaim Exp $
// --------------------------------------------------------------
//

#ifndef A01ParallelWorld_h
#define A01ParallelWorld_h 1

#include "globals.hh"
#include "G4VUserParallelWorld.hh"

class A01ParallelWorld : public G4VUserParallelWorld
{
public:
    A01ParallelWorld(G4String worldName);
    virtual ~A01ParallelWorld();

public:
    virtual void Construct();

};

#endif


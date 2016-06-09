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
// $Id: ExN07PrimaryGeneratorAction.hh,v 1.1 2003/03/10 01:43:36 asaim Exp $
// GEANT4 tag $Name: geant4-06-00 $
//

#ifndef ExN07PrimaryGeneratorAction_h
#define ExN07PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class ExN07PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExN07PrimaryGeneratorAction();    
    virtual ~ExN07PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun*                particleGun;
    G4bool                        serial;

  public:
    inline void SetSerial(G4bool ser)
    { serial = ser; }
};


#endif



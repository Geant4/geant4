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
// $Id: Tst01ParticleGun.hh,v 1.4 2006-06-29 19:31:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst01ParticleGun_h
#define Tst01ParticleGun_h 1


#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4PrimaryParticle.hh"
#include <vector>

class G4Event;
class G4DecayProducts;

class Tst01ParticleGunMessenger;

class Tst01ParticleGun : public G4ParticleGun
{
  public:
    Tst01ParticleGun();
    Tst01ParticleGun(G4int numberofparticles);
    Tst01ParticleGun(G4ParticleDefinition * particleDef, 
                   G4int numberofparticles = 1);

  public:
    virtual ~Tst01ParticleGun();
    Tst01ParticleGun(const Tst01ParticleGun &right);

    const Tst01ParticleGun & operator=(const Tst01ParticleGun &right);

  public:
    void GeneratePrimaryVertex(G4Event* evt);
    void SetDecayProperTime(G4double t);
    void AddDaughter(const G4PrimaryParticle *daughter);

  public:
    G4double GetDecayProperTime() const;

  protected:  
    virtual              void SetInitialValues();
    G4double             particle_decay_time;

    std::vector<const G4PrimaryParticle*>  decay_products;  

  private:
    Tst01ParticleGunMessenger* theMessenger;
};

 
inline
 void Tst01ParticleGun::SetDecayProperTime(G4double t)
{
  particle_decay_time = t;
}

inline
 G4double  Tst01ParticleGun::GetDecayProperTime() const
{
  return particle_decay_time;
}

inline 
 void Tst01ParticleGun::AddDaughter(const G4PrimaryParticle *daughter)
{
    decay_products.push_back(daughter);
}

#endif








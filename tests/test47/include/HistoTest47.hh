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
#ifndef HistoTest47_H
#define HistoTest47_H

#include "G4ParticleDefinition.hh"
#include "G4LorentzVector.hh"

#include <string>
#include <map>
#include <vector>

class G4VParticleChange;

class HistoTest47 {

public:

  HistoTest47(std::string namePart, std::string nameMat, G4double momentum,
	      std::string nameGen);
  virtual ~HistoTest47();

  void setParticle(std::string namePart);
  void setTarget(std::string nameMat)    {target = nameMat; unInitialized =true;}
  void setMomentum(double momentum)      {energy = momentum; unInitialized =true;}
  void setGenerator(std::string nameGen) {generator = nameGen; unInitialized =true;}
  void setJobID( int id ) { jobID = id ; return ; }
  void setClusterID( int id ) { clusterID = id ; return ; }
  virtual void fill(G4VParticleChange*, G4LorentzVector)=0;
  virtual void write(G4double cross_sec, G4int nevt)=0;

protected:

  G4int particleType(G4ParticleDefinition*);
  bool                                   unInitialized;
  std::string                            particle, target, generator;
  double                                 energy;
  int                                    jobID;
  int                                    clusterID;

private:

  std::map<G4ParticleDefinition*, G4int> mapParticle;
  
};

#endif

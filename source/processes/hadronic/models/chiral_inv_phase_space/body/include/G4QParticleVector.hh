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
// $Id: G4QParticleVector.hh,v 1.20 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for Decay Channel Vector in CHIPS model
// ---------------------------------------------------------------
// Short description: The G4QParticle is a part of the CHIPS World. It is
// characterized by the quark content, spin, mass, width and a vector of
// the decay channels (G4QDecayCannelVector).
// -----------------------------------------------------------------------

#ifndef G4QParticleVector_h
#define G4QParticleVector_h 1

#include "G4QParticle.hh"
#include <vector>

//typedef std::vector<G4QParticle *> G4QParticleVector;

struct DeleteQParticle
{
  void operator()(G4QParticle *aN)
  {
    if(aN)
    {
      //G4cout<<"G4QParticleVector::DeleteQParticle: aN="<<aN<<",PDG="<<aN->GetPDGCode()
      //      <<",Q="<<aN->GetQCode()<<G4endl; // TMP
      delete aN;
    }
    //else G4cout<<"***G4QParticleVector::DeleteQParticle: aN="<<aN<<G4endl; // TMP
  }
};

class G4QParticleVector : public  std::vector<G4QParticle *>
{
public:
  ~G4QParticleVector() 
  {
    std::for_each(begin(),end(),DeleteQParticle());
    clear();
  }
};

#endif



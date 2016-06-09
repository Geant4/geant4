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
// $Id: G4QParticleVector.hh,v 1.18 2005/06/04 13:08:23 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for Decay Channel Vector in CHIPS model
// ---------------------------------------------------------------

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



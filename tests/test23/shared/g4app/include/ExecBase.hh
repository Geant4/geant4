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
#ifndef ExecBase_H
#define ExecBase_H 1

// fwd declaration
//
class TstReader;
class G4VParticleChange;

class ExecBase
{

   public:
   
      // ctor & dtor
      ExecBase( const TstReader* pset ) { Init(pset); }
      virtual ~ExecBase() {};
      
      
      virtual G4VParticleChange* DoEvent() = 0;

   protected:
   
      ExecBase() {};
      virtual void InitSetup( const TstReader* ) = 0;
      virtual void InitBeam(  const TstReader* ) = 0;
      
      void InitParticles();

   private:
   
      void Init( const TstReader* );

};

#endif

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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name: G4QMDFermiG.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 18 November 2008
// -----------------------------------------------------------------------------
//

#ifndef G4QMDFermiG_hh
#define G4QMDFermiG_hh

#include "G4TheoFSGenerator.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4VPartonStringModel.hh"
#include "G4QGSModel.hh"

#include "G4QGSFGParticipants.hh" 
#include "G4HadronicInteraction.hh"

class G4QMDFermiG : public G4HadronicInteraction
{
   public:
      G4QMDFermiG();
      ~G4QMDFermiG();

      G4HadFinalState* ApplyYourself( const G4HadProjectile &aTrack , G4Nucleus & targetNucleus );

   private:

      G4VIntraNuclearTransportModel* theTransport;
      G4VPartonStringModel* theHEG;
};

#endif

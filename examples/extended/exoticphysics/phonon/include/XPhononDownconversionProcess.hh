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
/// \file exoticphysics/phonon/include/XPhononDownconversionProcess.hh
/// \brief Definition of the XPhononDownconversionProcess class
//
// $Id$
//
#ifndef XPhononDownconversionProcess_h
#define XPhononDownconversionProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"

#include "XPhysicalLattice.hh"

class XPhononDownconversionProcess : public G4VDiscreteProcess 
{
  public:


     XPhononDownconversionProcess(const G4String& processName ="XPhononDownconversionProcess" );

     virtual ~XPhononDownconversionProcess();

     virtual G4VParticleChange* PostStepDoIt(
                const G4Track&, const G4Step& );
 
     virtual G4bool IsApplicable(const G4ParticleDefinition&);
                           
  protected:

     virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );



 
  private:
    double fBeta, fGamma, fLambda, fMu;

  // hide assignment operator as private 
     XPhononDownconversionProcess(XPhononDownconversionProcess&);
     XPhononDownconversionProcess& operator=(const XPhononDownconversionProcess& right);

  // relative probability that anharmonic decay occurs L->L'+T'
  inline double GetLTDecayProb(G4double, G4double);
  inline double GetTTDecayProb(G4double, G4double);
  inline double MakeLDeviation(G4double, G4double);
  inline double MakeTTDeviation(G4double, G4double);
  inline double MakeTDeviation(G4double, G4double);
  void MakeTTSecondaries(const G4Track&);
  void MakeLTSecondaries(const G4Track&);



  XPhysicalLattice* Lattice;


};

#endif











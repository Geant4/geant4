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
// G4IntegrHadrNucleus.hh

#ifndef  G4IntegrHadrNucleus_h
#define  G4IntegrHadrNucleus_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4Nucleus.hh"
#include "G4HadronValues.hh"

         class G4IntegrHadrNucleus : public G4HadronValues
 {

   public:

       G4IntegrHadrNucleus() : G4HadronValues() {;}
      ~G4IntegrHadrNucleus() {;}

       G4double GetElasticCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 

       G4double GetTotalCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 

       G4double GetProductionCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 


       G4double GetInelasticCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 


       G4double GetQuasyElasticCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 

   private:  

       void GetIntegralCrSec(G4Nucleus * aNucleus);

       G4double  TotalCrSec,  InelCrSec,  ProdCrSec,  ElasticCrSec, 
                 QuasyElasticCrSec,  HadrEnergy; 

 };

#endif

 

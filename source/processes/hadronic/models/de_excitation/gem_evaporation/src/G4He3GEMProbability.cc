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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//
// J.M. Quesada (July 2009) C's and k's  values according to Furihata's paper 
// (based on notes added on proof in Dostrovskii's paper)

#include "G4He3GEMProbability.hh"

G4He3GEMProbability::G4He3GEMProbability() :
    G4GEMProbability(3,2,1.0/2.0) // A,Z,Gamma
{}

G4He3GEMProbability::~G4He3GEMProbability()
{}

//G4double G4He3GEMProbability::CCoeficient(G4int aZ) const
//{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    //	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};
    // C for He3 is equal to C for alpha times 4/3
//    G4double C = 0.0;
	
	
//    if (aZ <= 30) {
//        C = 0.10;
//    } else if (aZ <= 50) {
//        C = 0.1 + -((aZ-50.)/20.)*0.02;
//   } else if (aZ < 70) {
//        C = 0.08 + -((aZ-70.)/20.)*0.02;
//    } else {
//        C = 0.06;
//    }
//    return C*(4.0/3.0);
//}


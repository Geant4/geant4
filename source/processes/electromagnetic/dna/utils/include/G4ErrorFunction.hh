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
/*
 * G4ErrorFunction.hh
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
*/

/*
 Extracted from http://ab-initio.mit.edu/Faddeeva
 Steven G. Johnson, October 2012.

 Copyright © 2012 Massachusetts Institute of Technology

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge,
 publish, distribute, sublicense, and/or sell copies of the Software,
 and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
*/

#ifndef G4ERRORFUNCTION_HH_
#define G4ERRORFUNCTION_HH_

#include "globals.hh"
#include <vector>

class G4ErrorFunction {
public:
    G4ErrorFunction();
    virtual ~G4ErrorFunction();

    static G4double NormQuantile(G4double x);
    static G4double erfcx_y100(G4double x);
    static G4double erfcx(G4double x);
    static G4double erfc(G4double x);
    static G4double erfcInv(G4double x);
    static G4double erfcWxy(G4double c, G4double x, G4double y);

    static G4double Lambda(G4double x, G4double beta, G4double alpha);
};

#endif /* G4ERRORFUNCTION_HH_ */

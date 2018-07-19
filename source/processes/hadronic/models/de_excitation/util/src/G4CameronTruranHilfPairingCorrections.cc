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
// $Id: G4CameronTruranHilfPairingCorrections.cc 96634 2016-04-27 09:31:49Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up

#include "G4CameronTruranHilfPairingCorrections.hh"
#include <CLHEP/Units/SystemOfUnits.h>

// Data comes from:
// J.W. Truran, A.G.W. Cameron, and E. Hilf, 
// Proc. Int. Conf. on the Properties of Nuclei Far From the Beta-Stability,
// Leysin, Switzerland, August 31 - September 4, 1970, Vol.1, p. 275
// S(Z)
G4double G4CameronTruranHilfPairingCorrections::PairingZTable[] = 
{ // 93 from Z = 10 to Z = 42
-2.200, 0.   ,-2.120, 0.   ,-1.981, 0.   ,-1.491, 0.   ,-1.450, 0.,
-1.701, 0.   ,-1.344, 0.   ,-1.349, 0.   ,-1.397, 0.   ,-1.311, 0.,
-1.161, 0.   ,-1.201, 0.   ,-1.449, 0.   ,-1.331, 0.   ,-1.272, 0.,
-1.198, 0.   ,-1.340, 0.   ,-1.407, 0.   ,-1.287, 0.   ,-1.334, 0.,
-1.307, 0.   ,-1.128, 0.   ,-1.152, 0.   ,-1.139, 0.   ,-1.138, 0.,
-1.115, 0.   ,-1.070, 0.   ,-1.096, 0.   ,-1.123, 0.   ,-0.901, 0.,
-0.933, 0.   ,-0.714, 0.   ,-0.799, 0.   ,-0.840, 0.   ,-0.726, 0.,
-0.815, 0.   ,-0.715, 0.   ,-0.788, 0.   ,-0.793, 0.   ,-0.663, 0.,
-0.705, 0.   ,-0.711, 0.   ,-0.561, 0.   ,-0.694, 0.   ,-0.683, 0.,
-0.501, 0.   ,-0.491
};
// S(N)
G4double G4CameronTruranHilfPairingCorrections::PairingNTable[] = 
{ // 145 from N = 10 to N = 154
-2.400, 0.   ,-2.358, 0.   ,-2.057, 0.   ,-1.462, 0.   ,-1.592, 0.,
-1.528, 0.   ,-1.470, 0.   ,-1.310, 0.   ,-1.316, 0.   ,-1.265, 0.,
-1.279, 0.   ,-1.256, 0.   ,-1.285, 0.   ,-1.440, 0.   ,-1.517, 0.,
-1.486, 0.   ,-1.456, 0.   ,-1.471, 0.   ,-1.336, 0.   ,-1.341, 0.,
-1.278, 0.   ,-0.821, 0.   ,-0.814, 0.   ,-1.095, 0.   ,-1.147, 0.,
-1.295, 0.   ,-1.281, 0.   ,-1.245, 0.   ,-1.197, 0.   ,-1.227, 0.,
-1.291, 0.   ,-1.254, 0.   ,-1.310, 0.   ,-1.171, 0.   ,-1.092, 0.,
-1.062, 0.   ,-0.713, 0.   ,-0.822, 0.   ,-0.843, 0.   ,-0.968, 0.,
-1.117, 0.   ,-0.999, 0.   ,-0.877, 0.   ,-0.844, 0.   ,-0.889, 0.,
-0.729, 0.   ,-0.706, 0.   ,-0.623, 0.   ,-0.511, 0.   ,-0.773, 0.,
-0.662, 0.   ,-0.808, 0.   ,-0.889, 0.   ,-0.930, 0.   ,-0.771, 0.,
-0.751, 0.   ,-0.835, 0.   ,-0.658, 0.   ,-0.607, 0.   ,-0.657, 0.,
-0.695, 0.   ,-0.457, 0.   ,-0.345, 0.   ,-0.452, 0.   ,-0.648, 0.,
-0.681, 0.   ,-0.416, 0.   ,-0.545, 0.   ,-0.482, 0.   ,-0.481, 0.,
-0.611, 0.   ,-0.654, 0.   ,-0.557
};

G4CameronTruranHilfPairingCorrections::G4CameronTruranHilfPairingCorrections()
{
  for(size_t i=0; i<ZTableSize; ++i) { PairingZTable[i] *= CLHEP::MeV; }
  for(size_t i=0; i<NTableSize; ++i) { PairingNTable[i] *= CLHEP::MeV; }
}



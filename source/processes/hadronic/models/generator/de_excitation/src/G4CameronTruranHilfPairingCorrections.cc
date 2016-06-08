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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CameronTruranHilfPairingCorrections.cc,v 1.5 2001/10/05 16:13:42 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4CameronTruranHilfPairingCorrections.hh"



// Data comes from:
// J.W. Truran, A.G.W. Cameron, and E. Hilf, 
// Proc. Int. Conf. on the Properties of Nuclei Far From the Beta-Stability,
// Leysin, Switzerland, August 31 - September 4, 1970, Vol.1, p. 275
// S(Z)
const G4double G4CameronTruranHilfPairingCorrections::PairingZTable
[G4CameronTruranHilfPairingCorrections::ZTableSize] = { // 102
 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 
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
const G4double G4CameronTruranHilfPairingCorrections::PairingNTable
[G4CameronTruranHilfPairingCorrections::NTableSize] = { 
 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
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
-0.611, 0.   ,-0.654, 0.   ,-0.557, 0.
};


G4CameronTruranHilfPairingCorrections  G4CameronTruranHilfPairingCorrections::theInstance(10.0);

G4CameronTruranHilfPairingCorrections::G4CameronTruranHilfPairingCorrections(G4double dummy)
{
    G4double even_more_dummy = dummy;
    even_more_dummy/=2.;
}

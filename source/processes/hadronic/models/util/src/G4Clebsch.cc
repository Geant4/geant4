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

#include "G4ios.hh"
#include "G4Clebsch.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "Randomize.hh"

const G4int G4POWLOGFACTMAX = 512;

using namespace std;

G4double G4Clebsch::ClebschGordanCoeff(G4int twoJ1, G4int twoM1, 
                                       G4int twoJ2, G4int twoM2, 
                                       G4int twoJ)
{
  if(twoJ1 < 0 || twoJ2 < 0 || twoJ < 0 ||
     ((twoJ1-twoM1) % 2) || ((twoJ2-twoM2) % 2)) { return 0; }

  G4int twoM = twoM1 + twoM2;
  if(twoM1 > twoJ1 || twoM1 < -twoJ1 ||
     twoM2 > twoJ2 || twoM2 < -twoJ2 ||
     twoM > twoJ || twoM < -twoJ) { return 0; }

  // Checks limits on J1, J2, J3
  G4double triangle = TriangleCoeff(twoJ1, twoJ2, twoJ);
  if(triangle == 0) { return 0; }

  G4Pow* g4pow = G4Pow::GetInstance();
  G4double factor = g4pow->logfactorial((twoJ1 + twoM1)/2) +
                    g4pow->logfactorial((twoJ1 - twoM1)/2);
  factor += g4pow->logfactorial((twoJ2 + twoM2)/2) +
            g4pow->logfactorial((twoJ2 - twoM2)/2);
  factor += g4pow->logfactorial((twoJ + twoM)/2) +
            g4pow->logfactorial((twoJ - twoM)/2);
  factor *= 0.5;

  G4int kMin = 0;
  G4int sum1 = (twoJ1 - twoM1)/2;
  G4int kMax = sum1;
  G4int sum2 = (twoJ - twoJ2 + twoM1)/2;
  if(-sum2 > kMin) kMin = -sum2;
  G4int sum3 = (twoJ2 + twoM2)/2;
  if(sum3 < kMax) kMax = sum3;
  G4int sum4 = (twoJ - twoJ1 - twoM2)/2;
  if(-sum4 > kMin) kMin = -sum4;
  G4int sum5 = (twoJ1 + twoJ2 - twoJ)/2;
  if(sum5 < kMax) kMax = sum5;

  // sanity / boundary checks
  if(kMin < 0) {
    G4Exception("G4Clebsch::ClebschGordanCoeff()", "Clebsch001", 
		JustWarning, "kMin < 0");
    return 0;
  }
  if(kMax < kMin) {
    G4Exception("G4Clebsch::ClebschGordanCoeff()", "Clebsch002", 
		JustWarning, "kMax < kMin");
    return 0;
  }
  if(kMax >= G4POWLOGFACTMAX) {
    G4Exception("G4Clebsch::ClebschGordanCoeff()", "Clebsch003", 
		JustWarning, "kMax too big for G4Pow");
    return 0;
  }

  // Now do the sum over k
  G4double kSum = 0.;
  for(G4int k = kMin; k <= kMax; k++) {
    G4double sign = (k % 2) ? -1 : 1;
    kSum += sign * G4Exp(factor - g4pow->logfactorial(sum1-k) -
			 g4pow->logfactorial(sum2+k) -
			 g4pow->logfactorial(sum3-k) -
			 g4pow->logfactorial(sum4+k) -
			 g4pow->logfactorial(k) -
			 g4pow->logfactorial(sum5-k));
  }

  return triangle*sqrt(twoJ+1)*kSum;
}

G4double G4Clebsch::ClebschGordan(G4int twoJ1, G4int twoM1, 
				  G4int twoJ2, G4int twoM2, 
				  G4int twoJ)
{
  // ClebschGordanCoeff() will do all input checking
  G4double clebsch = ClebschGordanCoeff(twoJ1, twoM1, twoJ2, twoM2, twoJ);
  return clebsch*clebsch;
}

std::vector<G4double> 
G4Clebsch::GenerateIso3(G4int twoJ1, G4int twoM1, 
			G4int twoJ2, G4int twoM2, 
			G4int twoJOut1, G4int twoJOut2)
{
  std::vector<G4double> temp;

  // ---- Special cases first ----

  // Special case, both Jin are zero
  if (twoJ1 == 0 && twoJ2 == 0) {
    G4Exception("G4Clebsch::GenerateIso3()", "Clebsch010", 
		JustWarning, "both twoJ are zero");
    temp.push_back(0.);
    temp.push_back(0.);
    return temp;
  }

  G4int twoM3 = twoM1 + twoM2;

  // Special case, either Jout is zero
  if (twoJOut1 == 0) {  
    temp.push_back(0.);
    temp.push_back(twoM3);
    return temp;
  }
  if (twoJOut2 == 0) {
    temp.push_back(twoM3);
    temp.push_back(0.);
    return temp;
  }
  
  // Number of possible states, in 
  G4int twoJMinIn = std::max(std::abs(twoJ1 - twoJ2), std::abs(twoM3));
  G4int twoJMaxIn = twoJ1 + twoJ2;

  // Number of possible states, out
  G4int twoJMinOut = 9999;
  for(G4int i=-1; i<=1; i+=2) {
    for(G4int j=-1; j<=1; j+=2) {
      G4int twoJTmp= std::abs(i*twoJOut1 + j*twoJOut2);
      if(twoJTmp < twoJMinOut) twoJMinOut = twoJTmp;
    }
  }
  twoJMinOut = std::max(twoJMinOut, std::abs(twoM3));
  G4int twoJMaxOut = twoJOut1 + twoJOut2;

  // Possible in and out common states 
  G4int twoJMin  =  std::max(twoJMinIn, twoJMinOut);
  G4int twoJMax  =  std::min(twoJMaxIn, twoJMaxOut);
  if (twoJMin > twoJMax) {
    G4Exception("G4Clebsch::GenerateIso3()", "Clebsch020", 
		JustWarning, "twoJMin > twoJMax");
    return temp;
  }
  
  // Number of possible isospins
  G4int nJ = (twoJMax - twoJMin) / 2 + 1;

  // A few consistency checks
  
  if ( (twoJ1 == 0 || twoJ2 == 0) && twoJMin != twoJMax ) {
    G4Exception("G4Clebsch::GenerateIso3()", "Clebsch021", 
		JustWarning, "twoJ1 or twoJ2 = 0, but twoJMin != JMax");
    return temp;
  }

  // MGP ---- Shall it be a warning or an exception?
  if (nJ == 0) {
    G4Exception("G4Clebsch::GenerateIso3()", "Clebsch022", 
		JustWarning, "nJ is zero, no overlap between in and out");
    return temp;
  }

  // Loop over all possible combinations of twoJ1, twoJ2, twoM11, twoM2, twoJTot
  // to get the probability of each of the in-channel couplings

  std::vector<G4double> clebsch;
  G4double sum = 0.0;
  for(G4int twoJ=twoJMin; twoJ<=twoJMax; twoJ+=2) {
    sum += ClebschGordan(twoJ1, twoM1, twoJ2, twoM2, twoJ);
    clebsch.push_back(sum);
  }     

  // Consistency check
  if (static_cast<G4int>(clebsch.size()) != nJ) {
    G4Exception("G4Clebsch::GenerateIso3()", "Clebsch023", 
		JustWarning, "nJ inconsistency");
    return temp;
  }

  // Consistency check
  if (sum <= 0.) {
    G4Exception("G4Clebsch::GenerateIso3()", "Clebsch024", 
		JustWarning, "Sum of Clebsch-Gordan probabilities <=0");
    return temp;
  }

  // Generate a random twoJTot according to the Clebsch-Gordan pdf
  sum *= G4UniformRand();
  G4int twoJTot = twoJMin;
  for (G4int i=0; i<nJ; ++i) {
    if (sum < clebsch[i]) {
      twoJTot += 2*i;
      break;
    }
  }

  // Generate twoM3Out

  std::vector<G4double> mMin;
  mMin.push_back(-twoJOut1);
  mMin.push_back(-twoJOut2);

  std::vector<G4double> mMax;
  mMax.push_back(twoJOut1);
  mMax.push_back(twoJOut2);

  // Calculate the possible |J_i M_i> combinations and their probability

  std::vector<G4double> m1Out;
  std::vector<G4double> m2Out;

  const G4int size = 20;
  G4double prbout[size][size];

  G4int m1pos(0), m2pos(0);
  G4int j12;
  G4int m1pr(0), m2pr(0);

  sum = 0.;
  for(j12 = std::abs(twoJOut1-twoJOut2); j12<=(twoJOut1+twoJOut2); j12+=2)
  {
    m1pos = -1;
    for (m1pr = static_cast<G4int>(mMin[0]+.00001); m1pr <= mMax[0]; m1pr+=2)
    {
      m1pos++;
      if (m1pos >= size) {
        G4Exception("G4Clebsch::GenerateIso3()", "Clebsch025", 
		    JustWarning, "m1pos > size");
        return temp;
      }
      m1Out.push_back(m1pr);
      m2pos = -1;
      for (m2pr = static_cast<G4int>(mMin[1]+.00001); m2pr <= mMax[1]; m2pr+=2)
      {
        m2pos++;
        if (m2pos >= size)
        {
          G4Exception("G4Clebsch::GenerateIso3()", "Clebsch026", 
		      JustWarning, "m2pos > size");
          return temp;
        }
        m2Out.push_back(m2pr);

        if(m1pr + m2pr == twoM3) 
        {
          G4int m12 = m1pr + m2pr;
          G4double c12 = ClebschGordan(twoJOut1, m1pr, twoJOut2,m2pr, j12);
          G4double c34 = ClebschGordan(0,0,0,0,0);
          G4double ctot = ClebschGordan(j12, m12, 0, 0, twoJTot);
          G4double cleb = c12*c34*ctot;
          prbout[m1pos][m2pos] = cleb;
          sum += cleb;
        }
        else
        {
          prbout[m1pos][m2pos] = 0.;
        }
      }
    }
  }
  
  if (sum <= 0.) {
    G4Exception("G4Clebsch::GenerateIso3()", "Clebsch027", 
		JustWarning, "sum (out) <=0");
    return temp;
  }

  for (G4int i=0; i<size; i++) {
    for (G4int j=0; j<size; j++) {
      prbout[i][j] /= sum;
    }
  }

  G4double rand = G4UniformRand();

  G4int m1p, m2p;

  for (m1p=0; m1p<m1pos; m1p++) {
    for (m2p=0; m2p<m2pos; m2p++) {
      if (rand < prbout[m1p][m2p]) {
        temp.push_back(m1Out[m1p]);
        temp.push_back(m2Out[m2p]);
        return temp;
      }   
      else rand -= prbout[m1p][m2p];
    }     
  }   

  G4Exception("G4Clebsch::GenerateIso3()", "Clebsch028", 
	      JustWarning, "Should never get here");
  return temp;
}

G4double G4Clebsch::Weight(G4int twoJ1,  G4int twoM1, 
			   G4int twoJ2,  G4int twoM2, 
			   G4int twoJOut1, G4int twoJOut2)
{
  G4double value = 0.;
  
  G4int twoM = twoM1 + twoM2;

  G4int twoJMinIn = std::max(std::abs(twoJ1 - twoJ2), std::abs(twoM));
  G4int twoJMaxIn = twoJ1 + twoJ2;

  G4int twoJMinOut = std::max(std::abs(twoJOut1 - twoJOut2), std::abs(twoM));
  G4int twoJMaxOut = twoJOut1 + twoJOut2;

  G4int twoJMin = std::max(twoJMinIn,twoJMinOut);
  G4int twoJMax = std::min(twoJMaxIn,twoJMaxOut);

  for (G4int twoJ=twoJMin; twoJ<=twoJMax; twoJ+=2) {
    // ClebschGordan() will do all input checking
    value += ClebschGordan(twoJ1, twoM1, twoJ2, twoM2, twoJ);
  }

  return value;
}

G4double G4Clebsch::Wigner3J(G4double j1, G4double j2, G4double j3, 
			     G4double m1, G4double m2, G4double m3)
{
  //  G4Exception("G4Clebsch::Wigner3J()", "Clebsch030", JustWarning, 
  //  "G4Clebsch::Wigner3J with double arguments is deprecated. Please use G4int version.");
  G4int twoJ1 = (G4int) (2.*j1);
  G4int twoJ2 = (G4int) (2.*j2);
  G4int twoJ3 = (G4int) (2.*j3);
  G4int twoM1 = (G4int) (2.*m1);
  G4int twoM2 = (G4int) (2.*m2);
  G4int twoM3 = (G4int) (2.*m3);
  return Wigner3J(twoJ1, twoM1, twoJ2, twoM2, twoJ3, twoM3);
}

G4double G4Clebsch::Wigner3J(G4int twoJ1, G4int twoM1, 
                             G4int twoJ2, G4int twoM2, 
                             G4int twoJ3)
{
  G4double clebsch = ClebschGordanCoeff(twoJ1, twoM1, twoJ2, twoM2, twoJ3);
  if(clebsch == 0) return clebsch;
  if( (twoJ1-twoJ2+twoM1+twoM2)/2 % 2) clebsch = -clebsch;
  return clebsch / sqrt(twoJ3+1);
}

G4double G4Clebsch::Wigner3J(G4int twoJ1, G4int twoM1, 
                             G4int twoJ2, G4int twoM2, 
                             G4int twoJ3, G4int twoM3)
{
  if(twoM1 + twoM2 != -twoM3) return 0;
  G4double clebsch = ClebschGordanCoeff(twoJ1, twoM1, twoJ2, twoM2, twoJ3);
  if(clebsch == 0) return clebsch;
  if( (twoJ1-twoJ2-twoM3)/2 % 2) clebsch = -clebsch;
  return clebsch / sqrt(twoJ3+1);
}

G4double G4Clebsch::NormalizedClebschGordan(G4int twoJ, G4int twoM, 
					    G4int twoJ1, G4int twoJ2,
					    G4int twoM1, G4int twoM2)
{
  // Calculate the normalized Clebsch-Gordan coefficient, that is the prob 
  // of isospin decomposition of (J,m) into J1, J2, m1, m2

  G4double cleb = 0.;
  if(twoJ1 == 0 || twoJ2 == 0) return cleb; 
  
  // Loop over all J1,J2,Jtot,m1,m2 combinations
  G4double sum = 0.0;
  for(G4int twoM1Current=-twoJ1; twoM1Current<=twoJ1;  twoM1Current+=2) {
    G4int twoM2Current = twoM - twoM1Current;
    // ClebschGordan() will do all further input checking
    G4double prob = ClebschGordan(twoJ1, twoM1Current, twoJ2, 
				  twoM2Current, twoJ);
    sum += prob;
    if (twoM2Current == twoM2 && twoM1Current == twoM1) cleb += prob;
  }

  // Normalize probs to 1 
  if (sum > 0.) cleb /= sum; 

  return cleb;
}

G4double G4Clebsch::TriangleCoeff(G4int twoA, G4int twoB, G4int twoC)
{
  // TC(ABC) = sqrt[ (A+B-C)! (A-B+C)! (-A+B+C)! / (A+B+C+1)! ]
  // return 0 if the triad does not satisfy the triangle inequalities
  G4Pow* g4pow =  G4Pow::GetInstance();

  double val = 0;
  G4int i = twoA+twoB-twoC;
  // only have to check that i is even the first time
  if(i<0 || (i%2)) return 0;
  else val += g4pow->logfactorial(i/2);

  i = twoA-twoB+twoC;
  if(i<0) return 0;
  else val += g4pow->logfactorial(i/2);

  i = -twoA+twoB+twoC;
  if(i<0) return 0;
  else val += g4pow->logfactorial(i/2);

  i = twoA+twoB+twoC+2;
  if(i<0) return 0;
  return G4Exp(0.5*(val - g4pow->logfactorial(i/2)));
}

G4double G4Clebsch::Wigner6J(G4int twoJ1, G4int twoJ2, G4int twoJ3, 
		             G4int twoJ4, G4int twoJ5, G4int twoJ6)
{
  if(twoJ1 < 0 || twoJ2 < 0 || twoJ3 < 0 || 
     twoJ4 < 0 || twoJ5 < 0 || twoJ6 < 0) return 0;

  // There is a fast calculation (no sums or exps) when twoJ6 = 0, 
  // so permute to use it when possible
  if(twoJ6 == 0) {
    if(twoJ1 != twoJ5) return 0;
    if(twoJ2 != twoJ4) return 0;
    if(twoJ1+twoJ2 < twoJ3) return 0;
    if((twoJ1 > twoJ2) && (twoJ3 < (twoJ1-twoJ2))) return 0;
    if((twoJ2 > twoJ1) && (twoJ3 < (twoJ2-twoJ1))) return 0;
    if((twoJ1+twoJ2+twoJ3) % 2) return 0;
    return (((twoJ1+twoJ2+twoJ3)/2) % 2 ? -1. : 1.) /sqrt((twoJ1+1)*(twoJ2+1));
  }
  if(twoJ1 == 0) return Wigner6J(twoJ6, twoJ2, twoJ4, twoJ3, twoJ5, 0);
  if(twoJ2 == 0) return Wigner6J(twoJ1, twoJ6, twoJ5, twoJ4, twoJ3, 0);
  if(twoJ3 == 0) return Wigner6J(twoJ4, twoJ2, twoJ6, twoJ1, twoJ5, 0);
  if(twoJ4 == 0) return Wigner6J(twoJ3, twoJ2, twoJ1, twoJ6, twoJ5, 0);
  if(twoJ5 == 0) return Wigner6J(twoJ1, twoJ3, twoJ2, twoJ4, twoJ6, 0);

  // Check triangle inequalities and calculate triangle coefficients.
  // Also check evenness of sums
  G4Pow* g4pow =  G4Pow::GetInstance();
  double triangles = 0;
  G4int i;
  i =  twoJ1+twoJ2-twoJ3;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i =  twoJ1-twoJ2+twoJ3;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i = -twoJ1+twoJ2+twoJ3;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i =  twoJ1+twoJ2+twoJ3+2; if(i<0 || i%2) return 0; else triangles -= g4pow->logfactorial(i/2);
  i =  twoJ1+twoJ5-twoJ6;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i =  twoJ1-twoJ5+twoJ6;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i = -twoJ1+twoJ5+twoJ6;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i =  twoJ1+twoJ5+twoJ6+2; if(i<0 || i%2) return 0; else triangles -= g4pow->logfactorial(i/2);
  i =  twoJ4+twoJ2-twoJ6;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i =  twoJ4-twoJ2+twoJ6;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i = -twoJ4+twoJ2+twoJ6;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i =  twoJ4+twoJ2+twoJ6+2; if(i<0 || i%2) return 0; else triangles -= g4pow->logfactorial(i/2);
  i =  twoJ4+twoJ5-twoJ3;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i =  twoJ4-twoJ5+twoJ3;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i = -twoJ4+twoJ5+twoJ3;   if(i<0 || i%2) return 0; else triangles += g4pow->logfactorial(i/2);
  i =  twoJ4+twoJ5+twoJ3+2; if(i<0 || i%2) return 0; else triangles -= g4pow->logfactorial(i/2);
  triangles = G4Exp(0.5*triangles);
  
  // Prepare to sum over k. If we have made it this far, all of the following
  // sums must be non-negative and divisible by two

  // k must be >= all of the following sums:
  G4int sum1 = (twoJ1 + twoJ2 + twoJ3)/2;
  G4int kMin = sum1;
  G4int sum2 = (twoJ1 + twoJ5 + twoJ6)/2;
  if(sum2 > kMin) kMin = sum2;
  G4int sum3 = (twoJ4 + twoJ2 + twoJ6)/2;
  if(sum3 > kMin) kMin = sum3;
  G4int sum4 = (twoJ4 + twoJ5 + twoJ3)/2;
  if(sum4 > kMin) kMin = sum4;

  // and k must be <= all of the following sums:
  G4int sum5 = (twoJ1 + twoJ2 + twoJ4 + twoJ5)/2;
  G4int kMax = sum5;
  G4int sum6 = (twoJ2 + twoJ3 + twoJ5 + twoJ6)/2;
  if(sum6 < kMax) kMax = sum6;
  G4int sum7 = (twoJ1 + twoJ3 + twoJ4 + twoJ6)/2;
  if(sum7 < kMax) kMax = sum7;

  // sanity / boundary checks
  if(kMin < 0) {
    G4Exception("G4Clebsch::Wigner6J()", "Clebsch040", 
		JustWarning, "kMin < 0");
    return 0;
  }
  if(kMax < kMin) {
    G4Exception("G4Clebsch::Wigner6J()", "Clebsch041", 
		JustWarning, "kMax < kMin");
    return 0;
  }
  if(kMax >= G4POWLOGFACTMAX) {
    G4Exception("G4Clebsch::Wigner6J()", "Clebsch041", 
		JustWarning, "kMax too big for G4Pow");
    return 0;
  }

  // Now do the sum over k
  G4double kSum = 0.;
  G4double sign = (kMin % 2) ? -1 : 1;
  for(G4int k = kMin; k <= kMax; k++) {
    kSum += sign * G4Exp(g4pow->logfactorial(k+1) -
                            g4pow->logfactorial(k-sum1) -
                            g4pow->logfactorial(k-sum2) -
                            g4pow->logfactorial(k-sum3) -
                            g4pow->logfactorial(k-sum4) -
                            g4pow->logfactorial(sum5-k) -
                            g4pow->logfactorial(sum6-k) -
                            g4pow->logfactorial(sum7-k));
    sign *= -1;
  }
  return triangles*kSum;
}

G4double G4Clebsch::Wigner9J(G4int twoJ1, G4int twoJ2, G4int twoJ3, 
		             G4int twoJ4, G4int twoJ5, G4int twoJ6,
		             G4int twoJ7, G4int twoJ8, G4int twoJ9)
{
  if(twoJ1 < 0 || twoJ2 < 0 || twoJ3 < 0 || 
     twoJ4 < 0 || twoJ5 < 0 || twoJ6 < 0 ||
     twoJ7 < 0 || twoJ8 < 0 || twoJ9 < 0) return 0;

  if(twoJ9 == 0) {
    if(twoJ3 != twoJ6) return 0;
    if(twoJ7 != twoJ8) return 0;
    G4double sixJ = Wigner6J(twoJ1, twoJ2, twoJ3, twoJ5, twoJ4, twoJ7);
    if(sixJ == 0) return 0;
    if((twoJ2+twoJ3+twoJ4+twoJ7)/2 % 2) sixJ = -sixJ;
    return sixJ/sqrt((twoJ3+1)*(twoJ7+1));
  }
  if(twoJ1 == 0) return Wigner9J(twoJ9, twoJ6, twoJ3, twoJ8, twoJ5, twoJ2, twoJ7, twoJ4, twoJ1);
  if(twoJ2 == 0) return Wigner9J(twoJ7, twoJ9, twoJ8, twoJ4, twoJ6, twoJ5, twoJ1, twoJ3, twoJ2);
  if(twoJ4 == 0) return Wigner9J(twoJ3, twoJ2, twoJ1, twoJ9, twoJ8, twoJ7, twoJ6, twoJ5, twoJ4);
  if(twoJ5 == 0) return Wigner9J(twoJ1, twoJ3, twoJ2, twoJ7, twoJ9, twoJ8, twoJ4, twoJ6, twoJ5);
  G4int twoS = twoJ1+twoJ2+twoJ3+twoJ4+twoJ5+twoJ6+twoJ7+twoJ8+twoJ9;
  if(twoS % 2) return 0;
  G4double sign = (twoS/2 % 2) ? -1 : 1;
  if(twoJ3 == 0) return sign*Wigner9J(twoJ7, twoJ8, twoJ9, twoJ4, twoJ5, twoJ6, twoJ1, twoJ2, twoJ3);
  if(twoJ6 == 0) return sign*Wigner9J(twoJ1, twoJ2, twoJ3, twoJ7, twoJ8, twoJ9, twoJ4, twoJ5, twoJ6);
  if(twoJ7 == 0) return sign*Wigner9J(twoJ3, twoJ2, twoJ1, twoJ6, twoJ5, twoJ4, twoJ9, twoJ8, twoJ7);
  if(twoJ8 == 0) return sign*Wigner9J(twoJ1, twoJ3, twoJ2, twoJ4, twoJ6, twoJ5, twoJ7, twoJ9, twoJ8);

  // No element is zero: check triads now for speed
  G4int i;
  i =  twoJ1+twoJ2-twoJ3; if(i<0 || i%2) return 0;
  i =  twoJ1-twoJ2+twoJ3; if(i<0 || i%2) return 0;
  i = -twoJ1+twoJ2+twoJ3; if(i<0 || i%2) return 0;
  i =  twoJ4+twoJ5-twoJ6; if(i<0 || i%2) return 0;
  i =  twoJ4-twoJ5+twoJ6; if(i<0 || i%2) return 0;
  i = -twoJ4+twoJ5+twoJ6; if(i<0 || i%2) return 0;
  i =  twoJ7+twoJ8-twoJ9; if(i<0 || i%2) return 0;
  i =  twoJ7-twoJ8+twoJ9; if(i<0 || i%2) return 0;
  i = -twoJ7+twoJ8+twoJ9; if(i<0 || i%2) return 0;
  i =  twoJ1+twoJ4-twoJ7; if(i<0 || i%2) return 0;
  i =  twoJ1-twoJ4+twoJ7; if(i<0 || i%2) return 0;
  i = -twoJ1+twoJ4+twoJ7; if(i<0 || i%2) return 0;
  i =  twoJ2+twoJ5-twoJ8; if(i<0 || i%2) return 0;
  i =  twoJ2-twoJ5+twoJ8; if(i<0 || i%2) return 0;
  i = -twoJ2+twoJ5+twoJ8; if(i<0 || i%2) return 0;
  i =  twoJ3+twoJ6-twoJ9; if(i<0 || i%2) return 0;
  i =  twoJ3-twoJ6+twoJ9; if(i<0 || i%2) return 0;
  i = -twoJ3+twoJ6+twoJ9; if(i<0 || i%2) return 0;
  
  // Okay, have to do the full sum over 6J's
  // Find limits for K sum
  G4int twoKMax = twoJ1+twoJ9;
  if(twoJ4+twoJ8 < twoKMax) twoKMax = twoJ4+twoJ8;
  if(twoJ2+twoJ6 < twoKMax) twoKMax = twoJ2+twoJ6;
  G4int twoKMin = twoJ1-twoJ9;
  if(twoJ9-twoJ1 > twoKMin) twoKMin = twoJ9-twoJ1;
  if(twoJ4-twoJ8 > twoKMin) twoKMin = twoJ4-twoJ8;
  if(twoJ8-twoJ4 > twoKMin) twoKMin = twoJ8-twoJ4;
  if(twoJ2-twoJ6 > twoKMin) twoKMin = twoJ2-twoJ6;
  if(twoJ6-twoJ2 > twoKMin) twoKMin = twoJ6-twoJ2;
  if(twoKMin > twoKMax) return 0;

  G4double sum = 0;
  for(G4int twoK = twoKMin; twoK <= twoKMax; twoK += 2) {
    G4double value = Wigner6J(twoJ1, twoJ4, twoJ7, twoJ8, twoJ9, twoK);
    if(value == 0) continue;
    value *= Wigner6J(twoJ2, twoJ5, twoJ8, twoJ4, twoK, twoJ6);
    if(value == 0) continue;
    value *= Wigner6J(twoJ3, twoJ6, twoJ9, twoK, twoJ1, twoJ2);
    if(value == 0) continue;
    if(twoK % 2) value = -value;
    sum += value*G4double(twoK+1);
  }
  return sum;
}

G4double G4Clebsch::WignerLittleD(G4int twoJ, G4int twoM, G4int twoN, 
				  G4double cosTheta)
{
  if(twoM < -twoJ || twoM > twoJ || twoN < -twoJ || twoN > twoJ 
     || ((twoM % 2) != (twoJ % 2)) || ((twoN % 2) != (twoJ % 2))) 
    { return 0; }

     if(cosTheta == 1.0) { return G4double(twoM == twoN); }

  G4int kMin = 0;
  if(twoM > twoN) kMin = (twoM-twoN)/2;
  G4int kMax = (twoJ + twoM)/2;
  if((twoJ-twoN)/2 < kMax) kMax = (twoJ-twoN)/2;

  G4double lnCosHalfTheta = G4Log((cosTheta+1.)*0.5) * 0.5;
  G4double lnSinHalfTheta = G4Log((1.-cosTheta)*0.5) * 0.5;

  G4Pow* g4pow =  G4Pow::GetInstance();
  G4double d = 0;
  for(G4int k = kMin; k <= kMax; k++) {
    G4double logSum = 0.5*(g4pow->logfactorial((twoJ+twoM)/2) +
                           g4pow->logfactorial((twoJ-twoM)/2) +
                           g4pow->logfactorial((twoJ+twoN)/2) +
                           g4pow->logfactorial((twoJ-twoN)/2));
    logSum += -g4pow->logfactorial((twoJ+twoM)/2 - k) -
               g4pow->logfactorial((twoJ-twoN)/2 - k) -
               g4pow->logfactorial(k) - 
               g4pow->logfactorial(k+(twoN-twoM)/2);
    logSum += (twoJ+(twoM-twoN)/2 - 2*k)*lnCosHalfTheta + 
      (2*k + (twoN-twoM)/2)*lnSinHalfTheta;
    G4double sign = (k % 2) ? -1 : 1;
    d += sign * G4Exp(logSum);
  }
  return d;
}

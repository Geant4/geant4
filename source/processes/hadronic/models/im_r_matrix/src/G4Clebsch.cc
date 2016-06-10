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

#include "globals.hh"
#include "G4ios.hh"
#include "G4HadronicException.hh"
#include "G4Clebsch.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4HadTmpUtil.hh"

G4Clebsch::G4Clebsch()
{
  G4int nLogs = 101;
  logs.push_back(0.);
  G4int i;
  for (i=1; i<nLogs; i++)
    {
      G4double previousLog = logs.back();
      G4double value = previousLog + std::log((G4double)i);
      logs.push_back(value);
    }
}


G4Clebsch::~G4Clebsch() 
{  }


G4bool G4Clebsch::operator==(const G4Clebsch &right) const
{
  return (this == (G4Clebsch *) &right);
}


G4bool G4Clebsch::operator!=(const G4Clebsch &right) const
{
  return (this != (G4Clebsch *) &right);
}


G4double G4Clebsch::Weight(G4int isoIn1,  G4int iso3In1, 
			   G4int isoIn2,  G4int iso3In2, 
			   G4int isoOut1, G4int isoOut2) const
{
  G4double value = 0.;
  
  G4int an_m = iso3In1 + iso3In2;

  G4int jMinIn = std::max(std::abs(isoIn1 - isoIn2), std::abs(an_m));
  G4int jMaxIn = isoIn1 + isoIn2;

  G4int jMinOut = std::max(std::abs(isoOut1 - isoOut2), std::abs(an_m));
  G4int jMaxOut = isoOut1 + isoOut2;

  G4int jMin = std::max(jMinIn,jMinOut);
  G4int jMax = std::min(jMaxIn,jMaxOut);

  G4int j;
  for (j=jMin; j<=jMax; j+=2)
  {
    value += ClebschGordan(isoIn1,iso3In1, isoIn2,iso3In2, j);
  }

  return value;
}


G4double G4Clebsch::ClebschGordan(G4int isoIn1, G4int iso3In1, 
				  G4int isoIn2, G4int iso3In2, 
				  G4int jOut) const
{
  // Calculates Clebsch-Gordan coefficient

  G4double j1 = isoIn1 / 2.0;
  G4double j2 = isoIn2 / 2.0;
  G4double j3 = jOut / 2.0;

  G4double m_1 = iso3In1 / 2.0;
  G4double m_2 = iso3In2 / 2.0;
  G4double m_3 = - (m_1 + m_2);

  G4int n = G4lrint(m_3+j1+j2+.1);
  G4double argument = 2. * j3 + 1.;
  if (argument < 0.) 
    throw G4HadronicException(__FILE__, __LINE__, "G4Clebsch::ClebschGordan - sqrt of negative argument");
  G4double coeff = std::sqrt(argument) / (std::pow(-1.,n));
  G4double clebsch = coeff * Wigner3J(j1,j2,j3, m_1,m_2,m_3);
  G4double value = clebsch * clebsch;

//   G4cout << "ClebschGordan(" 
// 	 << isoIn1 << "," << iso3In1 << ","
// 	 << isoIn2 << "," << iso3In2 << "," << jOut
// 	 << ") = " << value << G4endl;

  return value;
}


G4double G4Clebsch::Wigner3J(G4double j1, G4double j2, G4double j3, 
			     G4double m_1, G4double m_2, G4double m_3) const
{
  // Calculates Wigner 3-j symbols

  G4double value = 0.;

  G4double sigma = j1 + j2 + j3;
  std::vector<G4double> n;
  n.push_back(-j1 + j2 + j3);      // n0
  n.push_back(j1 - m_1);            // n1
  n.push_back(j1 + m_1);            // n2
  n.push_back(j1 - j2 + j3);       // n3
  n.push_back(j2 - m_2);            // n4
  n.push_back(j2 + m_2);            // n5
  n.push_back(j1 + j2 - j3);       // n6
  n.push_back(j3 - m_3);            // n7
  n.push_back(j3 + m_3);            // n8

  // Some preliminary checks

  G4bool ok = true;
  size_t i;
  for(i=1; i<=3; i++)
  {
    G4double sum1 = n[i-1] + n[i+2] + n[i+5];
    G4double sum2 = n[3*i-1] + n[3*i-2] + n[3*i-3];
    if (sum1 != sigma || sum2 != sigma) ok = false;
    G4int j;
    for(j=1; j<=3; j++) 
    {
      if (n[i+3*j-4] < 0.) ok = false; 
    }
  }

  if (ok)
  {
    G4int iMin = 1;
    G4int jMin = 1;
    G4double smallest = n[0];

    // Find the smallest n
    for (i=1; i<=3; i++)
    {
      G4int j;
      for (j=1; j<=3; j++)
      {
	if (n[i+3*j-4] < smallest)
	{
	  smallest = n[i+3*j-4];
	  iMin = i;
	  jMin = j;
	}
      }
    }

    G4int sign = 1;

    if(iMin > 1)
    {
      for(G4int j=1; j<=3; ++j)
      {
	G4double tmp = n[j*3-3];
	n[j*3-3] = n[iMin+j*3-4];
	n[iMin+j*3-4] = tmp;
      }
      sign = (G4int) std::pow(-1.,sigma);
    }

    if (jMin > 1)
    {
      for(i=1; i<=3; i++)
      {
	G4double tmp = n[i-1];
	n[i-1] = n[i+jMin*3-4];
	n[i+jMin*3-4] = tmp;
      }
      sign *= (G4int) std::pow(-1.,sigma);
    }

    const std::vector<G4double>& logVector = logs;//GetLogs();
    size_t n1 = G4lrint(n[0]);

    // Some boundary checks
    G4int logEntries = logVector.size() - 1;
    for (i=0; i<n.size(); i++)
    {
      if (n[i] < 0. || n[i] > logEntries)
	 throw G4HadronicException(__FILE__, __LINE__, "G4Clebsch::Wigner3J - Outside logVector boundaries, n");
    }

    G4double r1 = n[0];
    G4double r2 = n[3];
    G4double r3 = n[6];
    G4double r4 = n[1];
    G4double r5 = n[4];
    G4double r6 = n[7];
    G4double r7 = n[2];
    G4double r8 = n[5];
    G4double r9 = n[8];

    G4double l1 = logVector[(G4int)r1];
    G4double l2 = logVector[(G4int)r2];
    G4double l3 = logVector[(G4int)r3];
    G4double l4 = logVector[(G4int)r4];
    G4double l5 = logVector[(G4int)r5];
    G4double l6 = logVector[(G4int)r6];
    G4double l7 = logVector[(G4int)r7];
    G4double l8 = logVector[(G4int)r8];
    G4double l9 = logVector[(G4int)r9];

    G4double sigma1 = sigma + 1.;
    if (sigma1 < 0. || sigma1 > logEntries)
      throw G4HadronicException(__FILE__, __LINE__, "G4Clebsch::Wigner3J - Outside logVector boundaries, sigma");

    G4double ls = logVector[static_cast<G4int>(sigma1+.00001)];
    G4double hlp1 = (l2 + l3 + l4 +l7 -ls -l1 -l5 -l9 -l6 -l8) / 2.;
    G4int expon = static_cast<G4int>(r6 + r8+.00001);
    G4double sgn = std::pow(-1., expon);
    G4double coeff = std::exp(hlp1) * sgn;

    G4int n61 = static_cast<G4int>(r6 - r1+.00001);
    if (n61 < 0. || n61 > logEntries)
      throw G4HadronicException(__FILE__, __LINE__, "G4Clebsch::Wigner3J - Outside logVector boundaries, n61");
    G4int n81 = static_cast<G4int>(r8 - r1+.00001);
    if (n81 < 0. || n81 > logEntries)
      throw G4HadronicException(__FILE__, __LINE__, "G4Clebsch::Wigner3J - Outside logVector boundaries, n81");

    G4double hlp2 = l6 - logVector[n61] + l8 - logVector[n81];
    G4double sum = std::exp(hlp2);
    std::vector<G4double> S;
    S.push_back(sum);
    n1 = (size_t)r1;
    for (i=1; i<=n1; i++)
    {
      G4double last = S.back();
      G4double den = i * (r6 - r1 + i) * (r8 - r1 + i);
      if (den == 0) 
	throw G4HadronicException(__FILE__, __LINE__, "G4Clebsch::Wigner3J - divide by zero");
      G4double data = -last * (r1 + 1.0 - i) * (r5 + 1.0 - i) * (r9 + 1. - i) / den;
      S.push_back(data);
      sum += data;
    }
    value = coeff * sum * sign;
  } // endif ok
  else
  {
  }


//  G4cout << "Wigner3j(" 
//	 << j1 << "," << j2 << "," << j3 << "," 
//	 << m1 << "," << m2 << "," << m3 << ") = " 
//	 << value
//	 << G4endl;

  return value;
}



std::vector<G4double> G4Clebsch::GenerateIso3(G4int isoIn1, G4int iso3In1, 
						G4int isoIn2, G4int iso3In2, 
						G4int isoA,   G4int isoB) const
{
  std::vector<G4double> temp;

  // ---- Special cases first ----

  // Special case, both Jin are zero
  if (isoIn1 == 0 && isoIn2 == 0)
  {
    G4cout << "WARNING: G4Clebsch::GenerateIso3 - both isoIn are zero" << G4endl;
    temp.push_back(0.);
    temp.push_back(0.);
    return temp;
  }

  G4int iso3 = iso3In1 + iso3In2;

  // Special case, either Jout is zero
  if (isoA == 0)
  {  
    temp.push_back(0.);
    temp.push_back(iso3);
    return temp;
  }
  if (isoB == 0)
  {
    temp.push_back(iso3);
    temp.push_back(0.);
    return temp;
  }
  
  // Number of possible states, in 
  G4int jMinIn = std::max(std::abs(isoIn1 - isoIn2), std::abs(iso3));
  G4int jMaxIn = isoIn1 + isoIn2;

  // Number of possible states, out
    
  G4int jMinOut = 9999;
  G4int jTmp, i, j;
 
   for(i=-1; i<=1; i+=2)
   {
     for(j=-1; j<=1; j+=2)
     {
       jTmp= std::abs(i*isoA + j*isoB);
       if(jTmp < jMinOut) jMinOut = jTmp;
     }
   }
   jMinOut = std::max(jMinOut, std::abs(iso3));
   G4int jMaxOut = isoA + isoB;

   // Possible in and out common states 
   G4int jMin  =  std::max(jMinIn, jMinOut);
   G4int jMax  =  std::min(jMaxIn, jMaxOut);
   if (jMin > jMax)
   {
     throw G4HadronicException(__FILE__, __LINE__,  "G4Clebsch::GenerateIso3 - jMin > JMax");
   }
   
   // Number of possible isospins
   G4int nJ = (jMax - jMin) / 2 + 1;

   // A few consistency checks
   
   if ( (isoIn1 == 0 || isoIn2 == 0) && jMin != jMax )
     throw G4HadronicException(__FILE__, __LINE__,  "G4Clebsch::GenerateIso3 - J1 or J2 = 0, but jMin != JMax");

   // MGP ---- Shall it be a warning or an exception?
   if (nJ == 0)
     throw G4HadronicException(__FILE__, __LINE__,  "G4Clebsch::GenerateIso3 - nJ is zero, no overlap between in and out");

   // Loop over all possible combinations of isoIn1, isoIn2, iso3In11, iso3In2, jTot
   // to get the probability of each of the in-channel couplings

   std::vector<G4double> clebsch;

   for(j=jMin; j<=jMax; j+=2)
     {
       G4double cg = ClebschGordan(isoIn1, iso3In1, isoIn2, iso3In2, j);
       clebsch.push_back(cg);
     }     

   // Consistency check
   if (static_cast<G4int>(clebsch.size()) != nJ)
     throw G4HadronicException(__FILE__, __LINE__,  "G4Clebsch::GenerateIso3 - nJ inconsistency");

   G4double sum = clebsch[0];
   
   for (j=1; j<nJ; j++)
   {
     sum += clebsch[j];
   }
   // Consistency check
   if (sum <= 0.)
     throw G4HadronicException(__FILE__, __LINE__,  "G4Clebsch::GenerateIso3 - Sum of Clebsch-Gordan probabilities <=0");

   // Generate a normalized pdf 

   std::vector<G4double> clebschPdf;
   G4double previous = clebsch[0];
   clebschPdf.push_back(previous/sum);
   for (j=1; j<nJ; j++)
   {
     previous += clebsch[j];
     G4double prob = previous / sum;
     clebschPdf.push_back(prob);
   }

   // Generate a random jTot according to the Clebsch-Gordan pdf
   G4double rand = G4UniformRand();
   G4int jTot = 0;
   for (j=0; j<nJ; j++)
   {
     G4bool found = false;
     if (rand < clebschPdf[j])
     {
       found = true;
       jTot = jMin + 2*j;
     }
     if (found) break;
   }

   // Generate iso3Out

   std::vector<G4double> mMin;
   mMin.push_back(-isoA);
   mMin.push_back(-isoB);

   std::vector<G4double> mMax;
   mMax.push_back(isoA);
   mMax.push_back(isoB);

   // Calculate the possible |J_i M_i> combinations and their probability

   std::vector<G4double> m1Out;
   std::vector<G4double> m2Out;

   const G4int size = 20;
   G4double prbout[size][size];

   G4int m1pos(0), m2pos(0);
   G4int j12;
   G4int m1pr(0), m2pr(0);

   sum = 0.;
   for(j12 = std::abs(isoA-isoB); j12<=(isoA+isoB); j12+=2)
   {
     m1pos = -1;
     for (m1pr = static_cast<G4int>(mMin[0]+.00001); m1pr <= mMax[0]; m1pr+=2)
     {
       m1pos++;
       if (m1pos >= size)
	 throw G4HadronicException(__FILE__, __LINE__,  "G4Clebsch::GenerateIso3 - m1pos > size");
       m1Out.push_back(m1pr);
       m2pos = -1;
       for (m2pr = static_cast<G4int>(mMin[1]+.00001); m2pr <= mMax[1]; m2pr+=2)
       {
	 m2pos++;
	 if (m2pos >= size)
	 {
	   throw G4HadronicException(__FILE__, __LINE__,  "G4Clebsch::GenerateIso3 - m2pos > size");
	 }
	 m2Out.push_back(m2pr);

	 if(m1pr + m2pr == iso3) 
	 {
	   G4int m12 = m1pr + m2pr;
	   G4double c12 = ClebschGordan(isoA, m1pr, isoB,m2pr, j12);
	   G4double c34 = ClebschGordan(0,0,0,0,0);
           G4double ctot = ClebschGordan(j12, m12, 0, 0, jTot);
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
   
   if (sum <= 0.)
     throw G4HadronicException(__FILE__, __LINE__, "G4Clebsch::GenerateIso3 - sum (out) <=0");

   for (i=0; i<size; i++)
   {
     for (j=0; j<size; j++)
     {
       prbout[i][j] /= sum;
     }
   }

   rand = G4UniformRand();

   G4int m1p, m2p;

   for (m1p=0; m1p<m1pos; m1p++)
   {
     for (m2p=0; m2p<m2pos; m2p++)
     {
       if (rand < prbout[m1p][m2p])
       {
	 temp.push_back(m1Out[m1p]);
	 temp.push_back(m2Out[m2p]);
	 return temp;
       }   
       else
       {
	 rand -= prbout[m1p][m2p];
       }
     }     
   }   

  throw G4HadronicException(__FILE__, __LINE__,  "Should never get here ");
  return temp;
}


G4double G4Clebsch::NormalizedClebschGordan(G4int J, G4int M, 
					    G4int J1, G4int J2,
					    G4int m_1, G4int m_2) const
{
  // Calculate the normalized Clebsch-Gordan coefficient, that is the prob 
  // of isospin decomposition of (J,m) into J1, J2, m1, m2

  G4double cleb = 0.;

  if(J1 == 0 || J2 == 0) return cleb; 
  
  G4double sum = 0.0;

  // Loop over all J1,J2,Jtot,m1,m2 combinations

  for(G4int m1Current=-J1; m1Current<=J1;  m1Current+=2) 
    {
      G4int m2Current = M - m1Current;
      
      G4double prob = ClebschGordan(J1, m1Current, J2, m2Current, J);
      sum += prob;
      if (m2Current == m_2 && m1Current == m_1) cleb += prob;
    }

  // Normalize probs to 1 
  if (sum > 0.) cleb /= sum; 

  return cleb;
}

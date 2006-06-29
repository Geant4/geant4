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
// $Id: G4DNACrossSectionDataSetTest.cc,v 1.2 2006-06-29 19:43:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Riccardo Capra <capra@ge.infn.it>
//
// History:
// -----------
// 30 Jun 2005  RC         Created
                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                                                                                                          
#include "G4DNACrossSectionDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include <fstream>
#include <unistd.h>

#define EnergyUMIS 2 //eV
#define DataUMIS 3 //barn

int main(int argc, char **argv)
{
 const char *inFile;

 if (argc>2)
 {
  G4cout << "Syntax: " << argv[0] << " [filename]" << G4endl;
  G4cout << G4endl;
  G4cout << "(!) NOTE: some files will be written in this directory (pattern: G4DNACrossSectionDataSetTestFile*.dat)" << G4endl;
  return 0;
 }
 else if (argc==2)
  inFile=argv[1];
 else
 {
  inFile="G4DNACrossSectionDataSetTestFile";
  std::ofstream out("./G4DNACrossSectionDataSetTestFile.dat");
  out << "1.00000000e+05   8.14416702e+02   6.76315789e+02   4.88340749e+02   8.39880339e+01   7.89635115e-02   2.06314024e+03\n"
         "2.00000000e+05   5.28359350e+02   4.63221618e+02   3.44771389e+02   7.30146866e+01   2.52677590e-01   1.40961972e+03\n\n"
     "    3.00000000e+05   3.95840046e+02   3.53745700e+02   2.65748336e+02\t5.98188864e+01 4.09845238e-01   1.07556281e+03\n"
         "4.00000000e+05   3.19343880e+02   2.88487842e+02   2.17686875e+02   5.03590031e+01   5.32869987e-01   8.76410470e+02\r\n"
         "5.00000000e+05   2.69156735e+02   2.44940828e+02   1.85347963e+02   4.34715722e+01   6.24418579e-01   7.43541517e+02\n"
         "1.00000000e+06   1.55452794e+02   1.43970277e+02   1.09881243e+02   2.64352725e+01   8.10553505e-01   4.36550140e+02     \n"
         "1.50000000e+06   1.11668214e+02   1.04282494e+02   7.98560948e+01\t1.93490562e+01   8.14261386e-01   3.15970120e+02\n"
         "2.00000000e+06   8.80312507e+01   8.26279732e+01   6.33649352e+01   1.54081437e+01   7.71586134e-01   2.50203889e+02\n"
      "   2.50000000e+06   7.30696756e+01 6.88018974e+01   5.28437021e+01   1.28756775e+01   7.19674081e-01   2.08310627e+02\n\r"
         "3.00000000e+06   6.26612557e+01   5.91740007e+01   4.55044481e+01   1.11003552e+01   6.69053959e-01   1.79109114e+02\n"
      "   3.50000000e+06   5.50134904e+01   5.20584986e+01   4.00630585e+01   9.78183471e+00   6.23068350e-01   1.57539951e+02\n"
         "4.00000000e+06   4.91088724e+01   4.65679165e+01   3.58555442e+01   8.76079088e+00   5.82043180e-01   1.40875167e+02\n\n    # comment\n"
         "4.50000000e+06   4.44233138e+01   4.21789727e+01   3.25041581e+01   7.94478293e+00   5.45740594e-01   1.27596968e+02\n#comment\r\n"
         "5.00000000e+06   4.05969944e+01   3.85988273e+01   2.97569104e+01   7.27564958e+00   5.13600786e-01 1.16741982e+02\n"
         "6.00000000e+06   3.47166562e+01   3.30808108e+01   2.55256187e+01   6.24407052e+00   4.59695847e-01   1.00026852e+02\n"
         "7.00000000e+06   3.03909338e+01   2.90156224e+01   2.24051103e+01   5.48359988e+00   4.16287624e-01   8.77115540e+01\n"
         "8.00000000e+06   2.70812879e+01   2.58941898e+01   2.00065525e+01\t4.89722388e+00   3.80739730e-01   7.82599938e+01\n"
         "9.00000000e+06   2.44599065e+01   2.34131260e+01   1.80955062e+01   4.43069001e+00   3.51212438e-01   7.07504411e+01\n"
         "1.00000000e+07   2.23229847e+01   2.13902553e+01   1.65394526e+01   4.05000347e+00   3.26180582e-01   6.46288766e+01 #comment     \n"
         "1.10000000e+07   2.05485697e+01   1.97086000e+01   1.52450755e+01   3.73311927e+00 3.04729195e-01   5.95400936e+01#comment\n"
         "1.20000000e+07   1.90475083e+01   1.82869443e+01   1.41492895e+01   3.46491210e+00   2.86140746e-01   5.52347951e+01\n"
         "1.30000000e+07   1.77638183e+01   1.70693653e+01   1.32099385e+01   3.23474125e+00 \t  2.69857202e-01   5.15477205e+01       \n"
    "     1.40000000e+07   1.66514478e+01   1.60114678e+01   1.23958510e+01   3.03496833e+00   2.55467589e-01   4.83492025e+01\n"
         "1.50000000e+07   1.56760259e+01   1.50837486e+01   1.16792471e+01   2.85979032e+00   2.42653835e-01   4.55414657e+01\n"
         "1.60000000e+07   1.48138030e+01   1.42630410e+01   1.10462104e+01   2.70486989e+00   2.31166719e-01   4.30590910e+01    \n\n";
         
  if (!out.good())
  {
   G4cout << argv[0] << ": Failed to create \"./G4DNACrossSectionDataSetTestFile.dat\" file" << G4endl;
   return 0;
  }         

  out.close();
  setenv("G4LEDATA", ".", 1);                                                                                                   
 }                                                                                                                              

 G4LogLogInterpolation *algorithm(new G4LogLogInterpolation);
 
 G4DNACrossSectionDataSet dataset(algorithm, EnergyUMIS, DataUMIS);
 
 G4cout << argv[0] << ": G4DNACrossSectionDataSet::Load: " << G4endl;

 try
 {
  if (!dataset.LoadData(inFile))
  {
   G4cout << " - Failed with false" << G4endl;
   return 0;
  }
 }
 catch (...)
 {
  G4cout << " - Failed with exception" << G4endl;
  return 0;
 }
 
 if (argc==1)
  unlink("./G4DNACrossSectionDataSetTestFile.dat");
  
 G4cout << " - Passed" << G4endl;
 
 G4cout << argv[0] << ": G4DNACrossSectionDataSet::PrintData: " << G4endl;
 dataset.PrintData();
 
 setenv("G4LEDATA", ".", 1);
 G4cout << argv[0] << ": G4DNACrossSectionDataSet::Save: " << G4endl;
 try
 {
  if (!dataset.SaveData("G4DNACrossSectionDataSetTestFile-2"))
  {
   G4cout << " - Failed with false" << G4endl;
   return 0;
  }
 }
 catch (...)
 {
  G4cout << " - Failed with exception" << G4endl;
  return 0;
 }
 
 G4DNACrossSectionDataSet dataset2(algorithm->Clone(), EnergyUMIS, DataUMIS);
 try
 {
  if (!dataset2.LoadData("G4DNACrossSectionDataSetTestFile-2"))
  {
   G4cout << " - Failed with file missing" << G4endl;
   return 0;
  }
 }
 catch (...)
 {
  G4cout << " - Failed with file missing" << G4endl;
  return 0;
 }
 
 unlink("./G4DNACrossSectionDataSetTestFile-2.dat");

 G4cout << " - Passed" << G4endl;
 
 G4cout << argv[0] << ": Comparing data:" << G4endl;
 
 size_t n(dataset.NumberOfComponents());
 
 if (n!=dataset2.NumberOfComponents())
 {
  G4cout << " - Failed different number of components" << G4endl;
  return 0;
 }
 
 while (n>0)
 {
  n--;
  
  const G4DataVector &e(dataset.GetComponent(n)->GetEnergies(0));
  const G4DataVector &e2(dataset2.GetComponent(n)->GetEnergies(0));
  const G4DataVector &d(dataset.GetComponent(n)->GetData(0));
  const G4DataVector &d2(dataset2.GetComponent(n)->GetData(0));

  if (e.size()!=e2.size() || d.size()!=d2.size() || e.size()!=d.size())
  {
   G4cout << " - Failed different number of rows" << G4endl;
   return 0;
  }
  
  size_t m(e.size());
  
  while (m>0)
  {
   m--;
   
   if (e[m]!=e2[m])
   {
    G4cout << " - Failed different energies (line " << (m+1) << " - " << e[m] << ", " << e2[m] << ")" << G4endl;
    return 0;
   }

   if (d[m]!=d2[m])
   {
    G4cout << " - Failed different data (line " << (m+1) << ", column " << (n+2) << " - " << d[m] << ", " << d2[m] << ")" << G4endl;
    return 0;
   }
  }
 }
 
 G4cout << " - Passed" << G4endl;
 G4cout << G4endl;
 G4cout << "(!) NOTE: Files G4DNACrossSectionDataSetTestFile*.dat were written and removed from this directory" << G4endl;
 
 return 0;
}

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
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries, 
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class: 
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//

#include "G4DNACPA100LogLogInterpolation.hh"

// Constructor

G4DNACPA100LogLogInterpolation::G4DNACPA100LogLogInterpolation()
{ }

// Destructor

G4DNACPA100LogLogInterpolation::~G4DNACPA100LogLogInterpolation()
{ }

G4VDataSetAlgorithm* G4DNACPA100LogLogInterpolation::Clone() const 
{ return new G4DNACPA100LogLogInterpolation; }


G4double G4DNACPA100LogLogInterpolation::Calculate(G4double x, G4int bin, 
  const G4DataVector& points, 
  const G4DataVector& data) const
{
  //G4cout << "G4DNACPA100LogLogInterpolation is performed (2 arguments) " << G4endl;
  G4int nBins = data.size() - 1;
//G4double oldresult = 0.;
  G4double value = 0.;
  if (x < points[0])
    {
      value = 0.;
    }
  else if (bin < nBins)
    {
      G4double e1 = points[bin];
      G4double e2 = points[bin+1];
      G4double d1 = data[bin];
      G4double d2 = data[bin+1];
// Check of e1, e2, d1 and d2 values to avoid floating-point errors when estimating the interpolated value below -- S.I., Jun. 2008
      if ((d1 > 0.) && (d2 > 0.) && (e1 > 0.) && (e2 > 0.))
        {
// Streamline the Log-Log Interpolation formula in order to reduce the required number of log10() function calls
// Variable oldresult contains the result of old implementation of Log-Log interpolation -- M.G.P. Jun. 2001
//       oldresult = (std::log10(d1)*std::log10(e2/x) + std::log10(d2)*std::log10(x/e1)) / std::log10(e2/e1);
//       oldresult = std::pow(10.,oldresult);
// Variable value contains the result of new implementation, after streamlining the math operation -- N.A.K. Oct. 2008
         value = std::log10(d1)+(std::log10(d2/d1)/std::log10(e2/e1)*std::log10(x/e1));
         value = std::pow(10.,value);
// Test of the new implementation result (value variable) against the old one (oldresult) -- N.A.K. Dec. 2008
//       G4double diffResult = value - oldresult;
//       G4double relativeDiff = 1e-11;
// Comparison of the two values based on a max allowable relative difference
//       if ( std::fabs(diffResult) > relativeDiff*std::fabs(oldresult) )
//        {
// Abort comparison when at least one of two results is infinite
//           if ((!std::isinf(oldresult)) && (!std::isinf(value)))
//            {
//              G4cout << "G4DNACPA100LogLogInterpolation> Old Interpolated Value is:" << oldresult << G4endl;
//              G4cout << "G4DNACPA100LogLogInterpolation> New Interpolated Value is:" << value << G4endl << G4endl;
//              G4cerr << "G4DNACPA100LogLogInterpolation> Error in Interpolation:" << G4endl;
//              G4cerr << "The difference between new and old interpolated value is:" << diffResult << G4endl << G4endl;
//            }
//        }
        }
      else value = 0.;
    }
  else
    {
      value = data[nBins];
    }
  return value;
}


// Nicolas A. Karakatsanis: New implementation of log-log interpolation after directly loading 
//                          logarithmic values from G4EMLOW dataset

G4double G4DNACPA100LogLogInterpolation::Calculate(G4double x, G4int bin, 
  const G4DataVector& points,
  const G4DataVector& data,
  const G4DataVector& log_points, 
  const G4DataVector& log_data) const
{
  G4int nBins = data.size() - 1;
  G4double value = 0.;
  G4double log_x = std::log10(x);
  if (x < points[0])
    {
      value = 0.;
    }
  else if (bin < nBins)
    {
      G4double log_e1 = log_points[bin];
      //G4double log_e2 = log_points[bin+1];
      G4double log_d1 = log_data[bin];
      G4double log_d2 = log_data[bin+1];
      
      //G4cout << "x = " << x << " , logx = " << log_x  << " , bin = " << bin << G4endl; 
      //G4cout << "e1 = " << points[bin] << " d1 = " << data[bin] << G4endl;
      //G4cout << "e2 = " << points[bin+1] << " d2 = " << data[bin+1] << G4endl;
      //G4cout << "loge1 = " << log_e1 << " logd1 = " << log_d1 << G4endl;
      //G4cout << "loge2 = " << log_e2 << " logd2 = " << log_d2 << G4endl;
      //G4cout << "interpol " << log_d1 + (log_d2 - log_d1)*(log_x - log_e1)/(log_e2 - log_e1) << " " << G4endl;


      // Values e1, e2, d1 and d2 are the log values of the corresponding
      // original energy and data values. Simple linear interpolation performed
      // on loagarithmic data should be equivalent to log-log interpolation

      // CPA100 specific
       //value = log_d1 + (log_d2 - log_d1)*(log_x - log_e1)/(log_e2 - log_e1);
      // value = log_d1;

       value = log_d2; // UPPER VALUE INTERPOLATION
       if (log_x == log_e1) value = log_d1; // IN CASE OF EQUALITY

// Delogarithmize to obtain interpolated value
      value = std::pow(10.,value);
   }
 else
   {
     value = data[nBins];
   }
   
  return value;
}

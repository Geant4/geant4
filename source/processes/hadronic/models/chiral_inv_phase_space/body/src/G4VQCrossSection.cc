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
//
// $Id: G4VQCrossSection.cc,v 1.1 2004-12-14 16:01:20 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// CHIPS virtual class: G4VQCrossSection for the collision cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-04
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 27-Nov-04
// 
//===============================================================================================
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************

///#define debug
//#define edebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4VQCrossSection.hh"

// Initialization of the
G4int     G4VQCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4VQCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4VQCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4VQCrossSection::lastCS=0.;// Last value of the Cross Section
G4double  G4VQCrossSection::lastP=0.;  // Last used in cross section TheMomentum
//G4int     G4VQCrossSection::lastN1=0;  // Last used in cross section TheNumOfBin1
//G4int     G4VQCrossSection::lastF1=0;  // Last used in cross section TheFirstBin1
//G4int     G4VQCrossSection::lastL1=0;  // Last used in cross section TheLastBin1
//G4int     G4VQCrossSection::lastN2=0;  // Last used in cross section TheNumOfBin1
//G4int     G4VQCrossSection::lastF2=0;  // Last used in cross section TheFirstBin1
//G4int     G4VQCrossSection::lastL2=0;  // Last used in cross section TheLastBin1
//G4double  G4VQCrossSection::lastBP=0.; // Last value of the Boundary Momentum
//G4double  G4VQCrossSection::lastMP=0.; // Last value of the Maximum Momentum

G4int     G4VQCrossSection::lastI=0;          // The last position in the DAMDB
G4double  G4VQCrossSection::tolerance=.001;   // The last position in the DAMDB

// Set the new tolerance (abs(p_old/p_new-1)<tolerance)
void G4VQCrossSection::setTolerance(G4double tol)
//   ============================================
{
		tolerance=tol;
}

// Gives the threshold energy for different isotopes (can be improved in the derived class)
G4double G4VQCrossSection::ThresholdEnergy(G4int Z, G4int N) {return Z*0.*N;}

// The main member function giving the collision cross section (P is in MeV/c, CS is in mb)
G4double G4VQCrossSection::GetCrossSection(G4double Momentum, G4int targZ, G4int targN )
{
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <G4double> colTH; // Vector of energy thresholds for the reaction
  static std::vector <G4double> colCS; // Vector of last cross sections for the reaction
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
  G4bool in=false;                     // By default the isotope is found in the DAMDB
  if(targN!=lastN || targZ!=lastZ)     // This nucleus was not the last used isotop
  {
    in = false;                        // Now by default the isotope isn't found in DAMDB  
    lastP    = 0.;                     // New momentum history (nothing to compare with)
    lastN    = targN;                  // The last N of the calculated nucleus
    lastZ    = targZ;                  // The last Z of the calculated nucleus
    lastI    = colN.size();            // Size of the Associative Memory DB in the heap
    if(lastI) for(G4int i=0; i<lastI; i++) if(colN[i]==targN && colZ[i]==targZ)
	   { // The nucleus is found in DAMDB
      in = true;                       // This is the case when the isotop is found in DB
      lastP  =colP [i];                // Last Momentum  (A-dependent)
      lastTH =colTH[i];                // Last THreshold (A-dependent)
      lastCS =colCS[i];                // Last CrossSect (A-dependent)
      if(Momentum<=lastTH) return 0.;  // Momentum is below the Threshold value
      else if(fabs(lastP/Momentum-1.)<tolerance) return lastCS*millibarn; // Use last CS
      lastI  = i;                      // Make the found isotope to be current isotope
      lastCS=CalculateCrossSection(-1,lastI,lastN,lastZ,Momentum);//read&update DB, calc.CS
      break;                           // Go out of the LOOP
	   }
	   if(!in)                            // This nucleus has not been calculated previously
	   {
      lastTH = ThresholdEnergy(targZ, targN); // The last Threshold Energy
      // lastI==colN.size() frome above
      lastCS = CalculateCrossSection(0,lastI,lastN,lastZ,Momentum);// calcCS + createDAMDB
#ifdef pdebug
      G4cout<<"G4VQCS::GetCrossSection: P="<<Momentum<<", Z="<<targZ<<",N="<<targN<<G4endl;
#endif
      colN.push_back(targN);
      colZ.push_back(targZ);
      colP.push_back(Momentum);
      colTH.push_back(lastTH);
      colCS.push_back(lastCS);
      return lastCS*millibarn;
	   } // End of creation of the new set of parameters
  } // End of parameters udate
  else if(Momentum<=lastTH) return 0.; // Momentum is below the Threshold value
  else if(fabs(lastP/Momentum-1.)<tolerance) return lastCS*millibarn; // Use the last CS
  else lastCS=CalculateCrossSection(1,lastI,lastN,lastZ,Momentum); // Update DB, calc. CS
  colP[lastI]=Momentum;
  colCS[lastI]=lastCS;
  return lastCS*millibarn;
}

G4double G4VQCrossSection::GetExchangeEnergy() {return 0.;}

G4double G4VQCrossSection::GetExchangeQ2(G4double) {return 0.;}

G4int G4VQCrossSection::GetExchangePDGCode(G4double, G4double) {return 0;}

G4double G4VQCrossSection::GetVirtualFactor(G4double nu, G4double Q2) {return 0.*nu*Q2;}

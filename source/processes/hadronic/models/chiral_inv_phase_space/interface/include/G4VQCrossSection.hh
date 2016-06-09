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
// $Id: G4VQCrossSection.hh,v 1.3 2005/11/30 16:26:42 mkossov Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// GEANT4 virtual class: G4VQCrossSection -- header file
// M.V. Kossov, CERN-ITEP(Moscow), 4-FEB-2004
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 27-Nov-04
//
// Short description: this G4 virtual class is made for the cross section
// classes of the CHIPS model, which calculate the cross section for the
// particular Element (virtual GetCrossSection member function). Each of the
// CHIPS cross section classes creates its own Dynamic Associative Memory
// Data Base (DAMDB) for the already used isotopes. For all of them thay use the
// same algorithm. Common member functions of this algorithm can be in this
// basic virtual class. Any CHIPS cross section class MUST inherit from this virtual
// G4VQCrossSection class. In the G4QCollision class the general G4VQCrossSection*
// pointer is connected to this or that CHIPS cross section class (depending on the
// projectile particle), so each of the CHIPS cross section class must be
// an evolving singletone. The singletone nature can not be realized in the
// virtual class. So each derived CS class must have
//  static G4VQCrossSection* GetPointer(); // Gives a pointer to the singletone
// static function, which is defined in the *.cc file as
//     // Returns Pointer to the G4VQCrossSection class
//     G4VQCrossSection* G4VQCrossSection::GetPointer()
//     {
//       static  G4QXCrossSection theCrossSection; //***Static body of the Cross Section***
//       return &theCrossSection;
//     }
// the line
//   //virtual static G4VQCrossSection* GetPointer(); // Gives a pointer to the singletone
// Reminds about this necesity, but in C++ the virtual static function can not be
// realised, so the static function can not be realised in the interface. Developers
// must take care of this themselves because this member fuction is called to get a pointer
// to the singletone in the G4QCollision class. So there is an agreement to
// make a separate CS class for each projectile particle, e.g. while the (pi-)d
// and (pi+)d (as well as [n,z] and [z,n]) cross sections) are almost equal,
// they must be calculated in different classes: G4QPiMinusCrossSection and
// G4QPiPlusCrossSections. For the ion-nuclear cross sections there should exist only
// one G4QIonCrossSection class with a huge (#0f isotopes times #of already produced
// ions) DAMDB or a general analitic formula with parameters. --- December 2004 ---
// -----------------------------------------------------------------------
// At present (25.11.04) for the test purposes this virtual class is created
// for ohly G4QPhotonCrossSection, G4QElectronCrossSection, G4QMuonCrossSection,
// G4QTauCrossSection and G4QProtonCrossSection (only for pp collisions now).
// ****************************************************************************************
// ********* This HEADER is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************

#ifndef G4VQCrossSection_h
#define G4VQCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include <vector>
#include "Randomize.hh"

class G4VQCrossSection
{
protected:

  G4VQCrossSection() {;} // for each particle a separate instance of G4QCollision should be
                         // used (and inside use a separate instance of G4Q*CrossSection)

public:
  virtual ~G4VQCrossSection() {;}// for each particle a separate instance of G4QCollision
  //@@ can be improved in future)// should be used and inside a separate istance of CS's

  virtual G4double GetCrossSection(G4double Momentum, G4int Z, G4int N);

  //virtual static G4VQCrossSection* GetPointer()=0; // Gives a pointer to the singletone

  static void setTolerance(G4double tol); // Set NewTolerance for TheSameCroSec

  virtual G4double ThresholdEnergy(G4int Z, G4int N); // Gives 0 by default

  // Define in the derived class, F=0 - create DAMDB, F=-1 - read DAMDB, F=1 - update DAMDB
  virtual G4double CalculateCrossSection(G4int F, G4int I, G4int Z, G4int N, G4double P)=0;

  virtual G4double GetLastTOTCS(); // Get the last calculated total cross-section

  virtual G4double GetLastQELCS(); // Get the last calculated quasi-elastic cross-section

  virtual G4double GetDirectPart(G4double Q2); // Direct interaction with quark-partons

  virtual G4double GetNPartons(G4double Q2); // #of quark-partons in non-perturbative PhSp

  // Subroutines for the t-chanel processes with a leader (DIS, Elastic, Quasielastic etc.)

  virtual G4double GetExchangeEnergy(); // Returns energy of the t-chanel particle (gam,pi)

  virtual G4double GetExchangeQ2(G4double nu); // Returns mass (-t or Q2) of the exchange

  virtual G4double GetVirtualFactor(G4double nu, G4double Q2); // Returns a ReductionFactor

  virtual G4double GetQEL_ExchangeQ2(); // Get randomized Q2 for quasi-elastic scattering

  virtual G4double GetNQE_ExchangeQ2(); // Get randomized Q2 for non quasi-elastic scat.

  virtual G4int GetExchangePDGCode(); // PDGCode of the Exchange Particle

  // Body: Basic Parameters of DAMDB (each derived class can add it's own values)
  // -----------------------------------------------------------------------------
  // The basic scheme of the DAMDB coveres the cross section for isotopes with fixed
  // Z (lastZ - number of protons) and N (lastN - number of neutrons) from the
  // Threshold momentum (TH) up to infinity. The cross section is first (Tab.1)
  // tabulated from the threshold till the boundary momentum (BP). The Tab.1 is
  // the function of the momentum (p) with the N1 elements. The N1 elements can be
  // not all different from zero. The first non-zero element is F1, the last non-zero
  // element is L1. If TH#0 the Tab.1 can be skipped. It is defined by N1=F1=L1=0 and
  // BP=TH. The Tab.1 is the function of the ln(p) with N2 elements (F2 is the first
  // non-zero element, L2 is the last non-zero element) from BP up tp MP. Both Tab.1
  // and Tab.2 are calculated when the projectile of the class meet the corresponding
  // ion. After that the tables are stored in the DAMDB for the fast calculations. To
  // avoid a complete calculation of the tables in the low energy calculation case,
  // the lastP momentum is used. The tables are calculated only till the momentum,
  // which already appeared in the simulation for this projectile and this isotope.
  // If the momentum above MP appeared, then the extrapolation function is calculated.
  // So, if lastP>MP it means that the cross section is defined for all energies above
  // TH. All parameters and pointers to arrays MUST be stored (F=0), updated (F=1) and
  // retrieved (F=-1) by the derived class in the CalculateCrossSection(F,I,N,Z,P)
  // function. The parameters are used for the immediate result: if the cross section is
  // calculated for the same Z, N, and fabs(p-lastP)/lastP<.001 (? - a parameter), the same
  // cross section (lastCS) is returned, if p<lastTH, then the 0 cross section is returned.
  // It helps to avoid double counting. The derived class can have only the approximation
  // functions, but such class is too slow, as it calculates the arythmetic equations each
  // time, when it is necessary to get a new cross section. So it is reasonable to
  // precalculate the tables, store them in memory, remember the pointers to these
  // functions and just interpolate them in the range of the most frequent energies (use
  // a LinearFit inline function of this virtual class for that). Starting some high
  // momentum (PM) the functional calculations are unavoidable, but fortunately they are
  // not frequent. In case of the ion-nuclear cross section the functional approach can
  // be reasonable, because tabulated cross-sections demand too much memory.
  //
  // -----------------------------------------------------------------------------
protected:
  G4double LinearFit(G4double X, G4int N, G4double* XN, G4double* YN);

  G4double EquLinearFit(G4double X, G4int N, G4double X0, G4double DX, G4double* Y);
protected:
  static G4int     lastN;    // The last N of calculated nucleus
  static G4int     lastZ;    // The last Z of calculated nucleus
  static G4double  lastP;    // Last used in the cross section Momentum
  static G4double  lastTH;   // Last value of the Momentum Threshold
  static G4double  lastCS;   // Last value of the Cross Section
  //static G4int     lastF1;   // Last used in the cross section TheFirstBin in Tab.1
  //static G4int     lastL1;   // Last used in the cross section TheLastBin in Tab.1
  //static G4int     lastN1;   // Last used in the cross section TheLastBin in Tab.1
  //static G4int     lastF2;   // Last used in the cross section TheFirstBin in Tab.2
  //static G4int     lastL2;   // Last used in the cross section TheLastBin in Tab.2
  //static G4int     lastN2;   // Last used in the cross section TheLastBin in Tab.2
  //static G4double  lastBP;   // Last value of the Boundary Momentum
  //static G4double  lastMP;   // Last value of the Maximum Momentum

private:
  static G4int     lastI;      // The last position in the DAMDB
  static G4double  tolerance;  // relative tolerance in momentum to get old CroSec
};

#endif

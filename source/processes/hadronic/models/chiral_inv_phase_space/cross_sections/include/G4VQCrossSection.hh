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
// $Id: G4VQCrossSection.hh,v 1.1 2009-11-16 18:15:42 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// Short description: a basic class for all CHIPS reaction cross-sections.
// -----------------------------------------------------------------------

#ifndef G4VQCrossSection_h
#define G4VQCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include <vector>
#include "Randomize.hh"

class G4VQCrossSection
{
protected:

  G4VQCrossSection() {;} // for each particle a separate instance of G4QCollision should be
                         // used (and inside use a separate instance of G4Q*CrossSection)

public:
  virtual ~G4VQCrossSection() {;}// for each particle separate instance of G4QXCrossSection
  //@@ can be improved in future)// should be used and inside a separate istance of CS's
  // Set the new tolerance (abs(p_old/p_new-1)<tolerance)
  static void setTolerance(G4double tol){tolerance=tol;}// Set NewTolerance for SameCrosSec

  // At present momentum (pMom) must be in GeV (@@ Units)
  virtual G4double GetCrossSection(G4bool, G4double, G4int, G4int, G4int pPDG=0)
                                                                   {return G4double(pPDG);}

  virtual G4double ThresholdEnergy(G4int Z, G4int N, G4int PDG=0); // Gives 0 by default

  // Define in the derived class, F=0 - create AMDB, F=-1 - read AMDB, F=1 - update AMDB
  virtual G4double CalculateCrossSection(G4bool CS, G4int F, G4int I, G4int PDG, G4int tgZ,
                                         G4int tgN, G4double pMom)=0;//*** PURE VIRTUAL ***

  virtual G4double GetLastTOTCS(); // LastCalculated total cross-section (total elastic)

  virtual G4double GetLastQELCS(); // LastCalculated quasielastic cross-section (quasifree)

  virtual G4double GetDirectPart(G4double Q2); // DirectInteraction with QuarkPartons (nuA)

  virtual G4double GetNPartons(G4double Q2); // #ofQuarkPartons in nonPerturbatPhaseSp(nuA)

  // Subroutines for the t-chanel processes with a leader (DIS, Elastic, Quasielastic etc.)

  virtual G4double GetExchangeEnergy(); // Returns energy of the t-chanel particle (gam,pi)

  virtual G4double GetExchangeT(G4int tZ, G4int tN, G4int pPDG); // -t=Q2 for hadronic

  virtual G4double GetSlope(G4int tZ, G4int tN, G4int pPDG); // B-slope of the maim maximum

  virtual G4double GetHMaxT();          // max(-t=Q2)/2 for hadronic (MeV^2)

  virtual G4double GetExchangeQ2(G4double nu=0); // Q2 for lepto-nuclear reactions

  virtual G4double GetVirtualFactor(G4double nu, G4double Q2); // ReductionFactor (leptA)

  virtual G4double GetQEL_ExchangeQ2(); // Get randomized Q2 for quasi-elastic scattering

  virtual G4double GetNQE_ExchangeQ2(); // Get randomized Q2 for non quasi-elastic scat.

  virtual G4int GetExchangePDGCode(); // PDGCode of the Exchange Particle (Pi0 by default)

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

  static G4double  tolerance;// relative tolerance in momentum to get old CroSec
};

#endif

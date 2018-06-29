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
// $Id: Pythia6.hh 100687 2016-10-31 11:20:33Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/include/Pythia6.hh
/// \brief Definition of the Pythia6 class

// 
// ----------------------------------------------------------------------------

// ******************************************************************************
// ******************************************************************************
// **                                                                          **
// **                                                                          **
// **              *......*                  Welcome to the Lund Monte Carlo!  **
// **         *:::!!:::::::::::*                                               **
// **      *::::::!!::::::::::::::*          PPP  Y   Y TTTTT H   H III   A    **
// **    *::::::::!!::::::::::::::::*        P  P  Y Y    T   H   H  I   A A   **
// **   *:::::::::!!:::::::::::::::::*       PPP    Y     T   HHHHH  I  AAAAA  **
// **   *:::::::::!!:::::::::::::::::*       P      Y     T   H   H  I  A   A  **
// **    *::::::::!!::::::::::::::::*!       P      Y     T   H   H III A   A  **
// **      *::::::!!::::::::::::::* !!                                         **
// **      !! *:::!!:::::::::::*    !!       This is PYTHIA version 6.418      **
// **      !!     !* -><- *         !!       Last date of change:  9 Jun 2008  **
// **      !!     !!                !!                                         **
// **      !!     !!                !!       Now is  0 Jan 2000 at  0:00:00    **
// **      !!                       !!                                         **
// **      !!        lh             !!       Disclaimer: this program comes    **
// **      !!                       !!       without any guarantees. Beware    **
// **      !!                 hh    !!       of errors and use common sense    **
// **      !!    ll                 !!       when interpreting results.        **
// **      !!                       !!                                         **
// **      !!                                Copyright T. Sjostrand (2008)     **
// **                                                                          **
// ** An archive of program versions and documentation is found on the web:    **
// ** http://www.thep.lu.se/~torbjorn/Pythia.html                              **
// **                                                                          **
// ** When you cite this program, the official reference is to the 6.4 manual: **
// ** T. Sjostrand, S. Mrenna and P. Skands, JHEP05 (2006) 026                 **
// ** (LU TP 06-13, FERMILAB-PUB-06-052-CD-T) [hep-ph/0603175].                **
// **                                                                          **
// ** Also remember that the program, to a large extent, represents original   **
// ** physics research. Other publications of special relevance to your        **
// ** studies may therefore deserve separate mention.                          **
// **                                                                          **
// ** Main author: Torbjorn Sjostrand; Department of Theoretical Physics,      **
// **   Lund University, Solvegatan 14A, S-223 62 Lund, Sweden;                **
// **   phone: + 46 - 46 - 222 48 16; e-mail: torbjorn@thep.lu.se              **
// ** Author: Stephen Mrenna; Computing Division, GDS Group,                   **
// **   Fermi National Accelerator Laboratory, MS 234, Batavia, IL 60510, USA; **
// **   phone: + 1 - 630 - 840 - 2556; e-mail: mrenna@fnal.gov                 **
// ** Author: Peter Skands; Theoretical Physics Department,                    **
// **   Fermi National Accelerator Laboratory, MS 106, Batavia, IL 60510, USA; **
// **   and CERN/PH, CH-1211 Geneva, Switzerland;                              **
// **   phone: + 41 - 22 - 767 24 59; e-mail: skands@fnal.gov                  **
// **                                                                          **
// **                                                                          **
// ******************************************************************************

#ifndef PYTHIA_Pythia6
#define PYTHIA_Pythia6

#include <vector>

int const KNDCAY  =  8000; //should be 4000 for pythia61

/// PYJETS common-block 
struct Pyjets_t
{
  int    N;
  int    NPAD;
  int    K[5][4000];
  double P[5][4000];
  double V[5][4000];
};

/// PYDAT1 common-block 
struct Pydat1_t
{
  int    MSTU[200];
  double PARU[200];
  int    MSTJ[200];
  double PARJ[200];
};

/// PYDAT3 common-block 
struct Pydat3_t
{
  int    MDCY[3][500];
  int    MDME[2][KNDCAY];
  double BRAT[KNDCAY];
  int    KFDP[5][KNDCAY];
};

/// Structure for Pythia6 particle properties
struct Pythia6Particle  
{
   Pythia6Particle(
     int ks, int kf, int parent, int firstChild, int lastChild,
     float px, float py, float pz, float energy, float mass,
     float vx, float vy, float vz, float time, float lifetime)
     : fKS(ks), fKF(kf), 
       fParent(parent), fFirstChild(firstChild), fLastChild(lastChild),
       fPx(px), fPy(py), fPz(pz), fEnergy(energy), fMass(mass),
       fVx(vx), fVy(vy), fVz(vz), fTime(time), fLifetime(lifetime) {}

   int    fKS;            // status of particle       ( LUJETS K[1] )
   int    fKF;            // KF flavour code          ( LUJETS K[2] )
   int    fParent;        // parrent's id             ( LUJETS K[3] )
   int    fFirstChild;    // id of first child        ( LUJETS K[4] )
   int    fLastChild;     // id of last  child        ( LUJETS K[5] )

   float  fPx;            // X momenta [GeV/c]        ( LUJETS P[1] )
   float  fPy;            // Y momenta [GeV/c]        ( LUJETS P[2] )
   float  fPz;            // Z momenta [GeV/c]        ( LUJETS P[3] )
   float  fEnergy;        // Energy    [GeV]          ( LUJETS P[4] )
   float  fMass;          // Mass      [Gev/c^2]      ( LUJETS P[5] )

   float  fVx;            // X vertex  [mm]           ( LUJETS V[1] )
   float  fVy;            // Y vertex  [mm]           ( LUJETS V[2] )
   float  fVz;            // Z vertex  [mm]           ( LUJETS V[3] )
   float  fTime;          // time of procuction [mm/c]( LUJETS V[4] )
   float  fLifetime;      // proper lifetime [mm/c]   ( LUJETS V[5] )
};   

typedef std::vector<Pythia6Particle*> ParticleVector;

/// The C++ interface class to Pythia6 
///
/// According to TPythia6 class from Root:
/// (The TPythia6 class is an interface class to F77 routines in Pythia6                //
/// CERNLIB event generators, written by T.Sjostrand.)                         
/// http://root.cern.ch/
/// see http://root.cern.ch/root/License.html
///
/// The complete Pythia6 documentation can be found at:
/// http://home.thep.lu.se/~torbjorn/pythiaaux/recent.html
/// 

class Pythia6
{
  public:

   // ****** constructors and destructor
   Pythia6();
   virtual ~Pythia6();

   static Pythia6 *Instance();

   // ****** TPYTHIA routines
   //
   int   Pycomp(int kf);
   void  Py1ent(int line, int kf, double pe, double theta, double phi);
   ParticleVector*  ImportParticles();
   int   ImportParticles(ParticleVector* particles, const char* option="");

   // ****** /PYDAT1/
   //
   void  SetMSTJ(int i, int m   ) { fPydat1->MSTJ[i-1] = m; }

   // ****** /PYDAT3/
   //
   int   GetMDCY(int i, int j) { return fPydat3->MDCY[j-1][i-1]; }
   int   GetKFDP(int i, int j) { return fPydat3->KFDP[j-1][i-1]; }
   void  SetMDCY(int i, int j, int m) { fPydat3->MDCY[j-1][i-1] = m; }
   void  SetMDME(int i, int j, int m) { fPydat3->MDME[j-1][i-1] = m; }
   
  private:
   static  Pythia6* fgInstance;

   ParticleVector*  fParticles;
   Pyjets_t*        fPyjets;
   Pydat1_t*        fPydat1;
   Pydat3_t*        fPydat3;
};

#endif


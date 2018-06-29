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
// $Id: Pythia6.cc 100687 2016-10-31 11:20:33Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/src/Pythia6.cc
/// \brief Implementation of the Pythia6 class

// ----------------------------------------------------------------------------
// According to TPythia6 class from Root:
// (The TPythia6 class is an interface class to F77 routines in Pythia6                //
// CERNLIB event generators, written by T.Sjostrand.)                         
// http://root.cern.ch/
// see http://root.cern.ch/root/License.html
//
// The complete Pythia6 documentation can be found at:
// http://home.thep.lu.se/~torbjorn/pythiaaux/recent.html
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

#include "Pythia6.hh"

#include <iostream>
#include <cstdlib>
#include <cstring>

#ifndef WIN32
# define pycomp pycomp_
# define py1ent py1ent_
# define type_of_call
#else
# define pycomp PYCOMP
# define py1ent PY1ENT
# define type_of_call _stdcall
#endif

// pythia6 functions
extern "C" {
  int  type_of_call pycomp(int *kf);
  void type_of_call py1ent(int&, int&, double&, double&, double&);
  void*  pythia6_common_address(const char*);
}

// Direct declaration of pythia6 common blocks
// extern "C" {
//   extern Pyjets_t pyjets_;
//   extern Pydat1_t pydat1_;
//   extern Pydat3_t pydat3_;
// }

Pythia6*  Pythia6::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Pythia6* Pythia6::Instance() 
{
/// Static access method

   if ( ! fgInstance ) fgInstance = new Pythia6();

   return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Pythia6::Pythia6()  
  : fParticles(0),
    fPyjets(0),
    fPydat1(0),
    fPydat3(0)
{
/// Pythia6 constructor: creates a vector of Pythia6Particle in which it will 
/// store all particles. Note that there may be only one functional Pythia6 
/// object at a time, so it's not use to create more than one instance of it.
  
   // Protect against multiple objects.   All access should be via the
   // Instance member function. 
   if ( fgInstance ) {
      std::cerr << "There's already an instance of Pythia6" << std::endl;
      exit (1);
   }   
  
   fParticles = new ParticleVector();

   // Initialize common-blocks 
   fPyjets = (Pyjets_t*) pythia6_common_address("PYJETS");
   fPydat1 = (Pydat1_t*) pythia6_common_address("PYDAT1");
   fPydat3 = (Pydat3_t*) pythia6_common_address("PYDAT3");

   // Alternative way to initialize common-blocks
   // usind direct declaration of pythia6 common blocks
   // fPyjets = &pyjets_;
   // fPydat1 = &pydat1_;
   // fPydat3 = &pydat3_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Pythia6::~Pythia6()
{
/// Destroy the object, delete and dispose all Pythia6Particles currently on 
/// list.

   if ( fParticles ) {
      ParticleVector::const_iterator it;
      for ( it = fParticles->begin(); it != fParticles->end(); it++ )
        delete  *it;
      delete fParticles;
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Pythia6::Pycomp(int kf) 
{
/// Interface with fortran routine pycomp

   return pycomp(&kf);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Pythia6::Py1ent(int ip, int kf, double pe, double theta, double phi)
{
/// Add one entry to the event record, i.e. either a parton or a
/// particle. 
///
///  IP:   normally line number for the parton/particle. There are two
///        exceptions:
/// 
///        If IP = 0: line number 1 is used and PYEXEC is called. 
///        If IP < 0: line -IP is used, with status code K(-IP,2)=2
///                   rather than 1; thus a parton system may be built
///                   up by filling all but the last parton of the
///                   system with IP < 0.   
///  KF:   parton/particle flavour code (PDG code)
///  PE:   parton/particle energy. If PE is smaller than the mass,
///        the parton/particle is taken to be at rest.  
///  THETA:
///  PHI:  polar and azimuthal angle for the momentum vector of the
///        parton/particle. 

   py1ent(ip, kf, pe, theta, phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Pythia6::ImportParticles(ParticleVector* particles, const char* option)
{
///  Default primary creation method. It reads the /HEPEVT/ common block which
///  has been filled by the GenerateEvent method. If the event generator does
///  not use the HEPEVT common block, This routine has to be overloaded by
///  the subclasses.
///  The function loops on the generated particles and store them in
///  the TClonesArray pointed by the argument particles.
///  The default action is to store only the stable particles (ISTHEP = 1)
///  This can be demanded explicitly by setting the option = "Final"
///  If the option = "All", all the particles are stored.

   if ( particles == 0 ) return 0;
   
   ParticleVector::const_iterator it;
   for ( it = particles->begin(); it != particles->end(); it++ )
     delete  *it;
   particles->clear();
   
   int numpart = fPyjets->N;
   int nparts=0;
   if (!strcmp(option,"") || !strcmp(option,"Final")) {
      for (int i = 0; i<numpart; i++) {

        if (fPyjets->K[0][i] == 1) {
          //
          //  Use the common block values for the TParticle constructor
          //
          particles->push_back(
            new Pythia6Particle(
                            fPyjets->K[0][i] ,
                            fPyjets->K[1][i] ,
                            fPyjets->K[2][i] ,
                            fPyjets->K[3][i] ,
                            fPyjets->K[4][i] ,
                            fPyjets->P[0][i] ,
                            fPyjets->P[1][i] ,
                            fPyjets->P[2][i] ,
                            fPyjets->P[3][i] ,
                            fPyjets->P[4][i] ,
                            fPyjets->V[0][i] ,
                            fPyjets->V[1][i] ,
                            fPyjets->V[2][i] ,
                            fPyjets->V[3][i] ,
                            fPyjets->V[4][i]));

          //     if(gDebug) printf("%d %d %d! ",i,fPyjets->K[1][i],numpart);
          nparts++;
       }
     }
   } 
   else if (!strcmp(option,"All")) {
      for (int i = 0; i<numpart; i++) {
          particles->push_back(
            new Pythia6Particle(
                            fPyjets->K[0][i] ,
                            fPyjets->K[1][i] ,
                            fPyjets->K[2][i] ,
                            fPyjets->K[3][i] ,
                            fPyjets->K[4][i] ,
                            fPyjets->P[0][i] ,
                            fPyjets->P[1][i] ,
                            fPyjets->P[2][i] ,
                            fPyjets->P[3][i] ,
                            fPyjets->P[4][i] ,
                            fPyjets->V[0][i] ,
                            fPyjets->V[1][i] ,
                            fPyjets->V[2][i] ,
                            fPyjets->V[3][i] ,
                            fPyjets->V[4][i]));
      }
      nparts=numpart;
   }

   return nparts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

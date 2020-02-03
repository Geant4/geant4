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
//
//
// class description:
//
// This is a concrete class of G4VPrimaryGenerator.
// This class object reads an ASCII file which contains particles generated
// by a physics generator which supports /HEPEVT/ common block.
// 
// The format of ASCII file must be equivalent to the following sample
// Fortran code:
//
//***********************************************************
//      SUBROUTINE HEP2G4
//*
//* Output /HEPEVT/ event structure to G4HEPEvtInterface
//*
//***********************************************************
//      PARAMETER (NMXHEP=2000)
//      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
//     >JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
//      DOUBLE PRECISION PHEP,VHEP
//*
//      WRITE(6,*) NHEP
//      DO IHEP=1,NHEP
//       WRITE(6,10) 
//     >  ISTHEP(IHEP),IDHEP(IHEP),JDAHEP(1,IHEP),JDAHEP(2,IHEP),
//     >  PHEP(1,IHEP),PHEP(2,IHEP),PHEP(3,IHEP),PHEP(5,IHEP)
//10    FORMAT(I4,I10,I5,I5,4(1X,D15.8))
//      ENDDO
//*
//      RETURN
//      END
//
// The position and time of the primary interaction must be set by the
// corresponding set methods of G4VPrimaryGenerator base class, otherwise
// zero will be set.

#ifndef G4HEPEvtInterface_h
#define G4HEPEvtInterface_h 1

#include <fstream>
#include <vector>

#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4HEPEvtParticle.hh"

class G4PrimaryVertex;
class G4Event;

class G4HEPEvtInterface:public G4VPrimaryGenerator
{
  public: // with description
    G4HEPEvtInterface(const char* evfile, G4int vl=0);
//    G4HEPEvtInterface(G4String evfile);
    // Constructors, "evfile" is the file name (with directory path).
  
  public:
    ~G4HEPEvtInterface();

    void GeneratePrimaryVertex(G4Event* evt);

  private:
    G4int vLevel;
    G4String fileName;
    std::ifstream inputFile;
    std::vector<G4HEPEvtParticle*> HPlist;
};

#endif

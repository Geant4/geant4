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
// $Id: G4QChipolino.hh,v 1.24 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QChipolino ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Chipolino (Double Hadron) in CHIPS Model
// ----------------------------------------------------------
// Short description: In the CHIPS model not only hadrons are considered,
// but the di-hadrons, which can not be convereged to the quark content
// of only one hadron (e.g. pi+pi+, K+p, Delta++p etc). This kind of
// hadronic states, which can be easily decayed in two hadrons, is called
// Chipolino-particle in the model.
// ---------------------------------------------------------------------- 

#ifndef G4QChipolino_h
#define G4QChipolino_h 1

#include "globals.hh"
#include "G4QPDGCode.hh"

class G4QChipolino
{
public:
  // Constructors
  G4QChipolino();                                      // Default Constructor
  G4QChipolino(G4QContent& QContent);                  // Construction by Quark Content
  G4QChipolino(const G4QChipolino& right);             // Copy constructor by value
  G4QChipolino(G4QChipolino* right);                   // Copy constructor by pointer

  ~G4QChipolino();                                     // Public Destructor

  // Operators
  const G4QChipolino& operator=(const G4QChipolino& right);
  G4bool              operator==(const G4QChipolino& right) const;
  G4bool              operator!=(const G4QChipolino& right) const;

  // Selectors
  G4double              GetMass();            // Get ChipolinoMass (MinDoubleHadronMass)
  G4double              GetMass2();           // Get mass^2 of the Chipolino
  G4QPDGCode            GetQPDG1();           // Get 1-st QPDG of the Chipolino
  G4QPDGCode            GetQPDG2();           // Get 2-nd QPDG of the Chipolino
  G4QContent            GetQContent() const;  // Get private quark content of the Chipolino
  G4QContent            GetQContent1() const; // Get 1-st quark content of the Chipolino
  G4QContent            GetQContent2() const; // Get 2-st quark content of the Chipolino

  // Modifiers
  void SetHadronQPDG(const G4QPDGCode& QPDG); // Set QPDG of 1-st Hadron of the Chipolino
  void SetHadronPDGCode(const G4int& PDGCode);// Set PDGCode of 1st Hadron of the Chipolino
  void SetHadronQCont(const G4QContent& QC);  // Set QContent of 1st Hadron of theChipolino

private:  
  G4QPDGCode            theQPDG1;             // QPDG of the 1-st Hadron of the Chipolino
  G4QPDGCode            theQPDG2;             // QPDG of the 2-nd Hadron of the Chipolino
  G4QContent            theQCont;             // QuarkContent of the whole Chipolino
  G4QContent            theQCont1;            // QuarkCont. of the 1st Hadron of Chipolino
  G4double              minM;                 // Minimal Mass of Chipolino
};

std::ostream& operator<<(std::ostream& lhs, G4QChipolino& rhs);
//std::ostream& operator<<(std::ostream& lhs, const G4QChipolino& rhs);
inline G4bool G4QChipolino::operator==(const G4QChipolino& rhs) const {return this==&rhs;}
inline G4bool G4QChipolino::operator!=(const G4QChipolino& rhs) const {return this!=&rhs;}
 
inline G4double   G4QChipolino::GetMass()      {return minM;}
inline G4double   G4QChipolino::GetMass2()     {return minM*minM;}
inline G4QPDGCode G4QChipolino::GetQPDG1()     {return theQPDG1;}
inline G4QPDGCode G4QChipolino::GetQPDG2()     {return theQPDG2;}
inline G4QContent G4QChipolino::GetQContent1() const {return theQCont1;}
inline G4QContent G4QChipolino::GetQContent2() const {return theQCont - theQCont1;}
inline G4QContent G4QChipolino::GetQContent()  const {return theQCont;}


inline void G4QChipolino::SetHadronQPDG(const G4QPDGCode& newQPDG)
{
  theQPDG1  = newQPDG;
  theQCont1 = theQPDG1.GetQuarkContent();
  G4QContent theQCont2 = theQCont - theQCont1;
  theQPDG2  = G4QPDGCode(theQCont2.GetSPDGCode());
}

inline void G4QChipolino::SetHadronPDGCode(const G4int& PDGCode)
{
  theQPDG1.SetPDGCode(PDGCode);
  SetHadronQPDG(theQPDG1);
}
inline void G4QChipolino::SetHadronQCont(const G4QContent& QCont)
                          {SetHadronPDGCode(QCont.GetSPDGCode());}

#endif






// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QChipolino.hh,v 1.2 2000-09-10 13:58:55 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QChipolino_h
#define G4QChipolino_h 1

// ------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QChipolino ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Chipolino (Double Hadron) in CHIPS Model
// ----------------------------------------------------------

#include "globals.hh"
#include "G4QContent.hh"
#include "G4QPDGCode.hh"

class G4QChipolino
{
public:
  // Constructors
  G4QChipolino();                                      // Default Constructor
  G4QChipolino(G4QContent& QContent);                  // Construction by Quark Content
  G4QChipolino(const G4QChipolino& right);             // Copy constructor by value
  G4QChipolino(G4QChipolino* right);                   // Copy constructor by pointer

  ~G4QChipolino();                                     // Destructor

  // Operators
  const G4QChipolino& operator=(const G4QChipolino& right);
  G4int              operator==(const G4QChipolino& right) const;
  G4int              operator!=(const G4QChipolino& right) const;

  // Selectors
  G4double              GetMass();            // Get mass of the Chipolino (MinDoubleHadronMass)
  G4double              GetMass2();           // Get mass^2 of the Chipolino
  G4QPDGCode            GetQPDG1();           // Get 1-st QPDG of the Chipolino
  G4QPDGCode            GetQPDG2();           // Get 2-nd QPDG of the Chipolino
  G4QContent            GetQContent() const;  // Get private quark content of the Chipolino
  G4QContent            GetQContent1() const; // Get 1-st quark content of the Chipolino
  G4QContent            GetQContent2() const; // Get 2-st quark content of the Chipolino

  // Modifiers
  void SetHadronQPDG(const G4QPDGCode& QPDG); // Set QPDG of 1-st Hadron of the Chipolino
  void SetHadronPDGCode(const G4int& PDGCode);// Set PDGCode of 1-st Hadron of the Chipolino
  void SetHadronQCont(const G4QContent& QC);  // Set QContent of 1-st Hadron of the Chipolino

private:  
  G4QPDGCode            theQPDG1;             // QPDG of the 1-st Hadron of the Chipolino
  G4QPDGCode            theQPDG2;             // QPDG of the 2-nd Hadron of the Chipolino
  G4QContent            theQCont;             // Quark Content of the whole Chipolino
  G4QContent            theQCont1;            // Quark Content of the 1-st Hadron of Chipolino
  G4double              minM;                 // Minimal Mass of Chipolino
};

ostream& operator<<(ostream& lhs, G4QChipolino& rhs);
//ostream& operator<<(ostream& lhs, const G4QChipolino& rhs);
inline G4int G4QChipolino::operator==(const G4QChipolino& right) const {return this==&right;}
inline G4int G4QChipolino::operator!=(const G4QChipolino& right) const {return this!=&right;}
 
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






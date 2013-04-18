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
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN02SteppingAction_h
#define ExN02SteppingAction_h 1

#include "G4UserSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4String.hh"
#include "G4Types.hh"
#include <list>
#include <sstream>
class G4VLoggingFilter {
public:
    virtual ~G4VLoggingFilter() {};
    virtual G4String Name() const = 0;
    virtual G4bool AcceptStep( const G4Step* aStep ) = 0;
    virtual G4VLoggingFilter* Clone() const = 0;
    virtual void  Dump(std::ostream& out) const;
};
std::ostream& operator<<(std::ostream& out, const G4VLoggingFilter& logfileter);

class G4StepInfo {
private:
    G4String name;
    std::list<G4VLoggingFilter*> filters;
    G4int counter;
    G4double sum;
    G4double sum2;
    G4double sum3;
    G4double sum4;
    std::list<float> vals;
    static G4bool debug;
    void Clear();
public:
    G4StepInfo(const G4StepInfo& rhs);
    G4StepInfo(const char* aName);
    ~G4StepInfo();
    void AddFilter(G4VLoggingFilter* f);
    G4bool LogStep(const G4Step* aStep , G4double weigth = 1 );
    void SetName( const char* newName ) { name = newName; }
    void Dump(std::ostream& out) const;
    G4double Mean() const ;
    G4double Variance() const ;
    G4double ExcKurtosis() const ;
    G4double Skewness() const ;
    G4int GetCounter() const { return counter; }
};
std::ostream& operator<<(std::ostream& out, const G4StepInfo& stepinfo);

class G4DistributionInfo {
private:
    G4String name;
    long counter;
    std::list<G4StepInfo*> infos;
public:
    void AddInfoLogger( G4StepInfo* logger ) { infos.push_back(logger); }
    void Dump( std::ostream& out ) const;
    G4DistributionInfo( const char* aName ) : name(aName) { }
    virtual ~G4DistributionInfo();
    G4bool LogStep( const G4Step* aStep , G4double quantity );
};
std::ostream& operator<<(std::ostream& out,const G4DistributionInfo& distinfo);

class ExN02SteppingAction : public G4UserSteppingAction
{
  public:
    ExN02SteppingAction();
   ~ExN02SteppingAction();
    void Dump() const;
    void UserSteppingAction(const G4Step*);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

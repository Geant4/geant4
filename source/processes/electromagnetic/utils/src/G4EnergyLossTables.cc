// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//  
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnergyLossTables.cc,v 1.9 1999-11-12 14:10:19 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------------- 
// first version created by P.Urban , 06/04/1998
// modifications + "precise" functions added by L.Urban , 27/05/98
// modifications , TOF functions , 26/10/98, L.Urban
// cache mechanism in order to gain time, 11/02/99, L.Urban
// bug fixed , 12/04/99 , L.Urban
// 10/11/99: moved from RWT hash dictionary to STL map, G.Barrand, M.Maire
// -------------------------------------------------------------------

#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EnergyLossTablesHelper G4EnergyLossTables::t  ;
const G4ParticleDefinition* G4EnergyLossTables::lastParticle = NULL ; 
G4double G4EnergyLossTables::QQPositron = eplus*eplus ;
G4double G4EnergyLossTables::Chargesquare ;
G4int    G4EnergyLossTables::oldIndex = -1 ;
G4double G4EnergyLossTables::rmin = 0. ;
G4double G4EnergyLossTables::rmax = 0. ;
G4double G4EnergyLossTables::Thigh = 0. ;

G4EnergyLossTables::helper_map G4EnergyLossTables::dict;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EnergyLossTables::Register(
  const G4ParticleDefinition* p,
  const G4PhysicsTable* tDEDX,
  const G4PhysicsTable* tRange,
  const G4PhysicsTable* tInverseRange,
  const G4PhysicsTable* tLabTime,
  const G4PhysicsTable* tProperTime,
  G4double lowestKineticEnergy,
  G4double highestKineticEnergy,
  G4double massRatio,
  G4int NumberOfBins)
{
  dict[p]= G4EnergyLossTablesHelper(tDEDX, tRange,tInverseRange,
                    tLabTime,tProperTime,lowestKineticEnergy,
		    highestKineticEnergy, massRatio,NumberOfBins);
           
  t = GetTables(p) ;    // important for cache !!!!!
  lastParticle = p ;
  Chargesquare = (p->GetPDGCharge())*(p->GetPDGCharge())/
                  QQPositron ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 G4double G4EnergyLossTables::GetPreciseDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial)
{
  if( aParticle != lastParticle)
  {
    t= GetTables(aParticle);
    lastParticle = aParticle;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
  }
  const G4PhysicsTable*  dEdxTable= t.theDEDXTable;

  G4int materialIndex = aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t.theMassRatio;
  G4double dEdx;
  G4bool isOut;

  if (scaledKineticEnergy<t.theLowestKineticEnergy) {

     dEdx = (*dEdxTable)(materialIndex)->GetValue(
              t.theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>t.theHighestKineticEnergy) {

     dEdx = (*dEdxTable)(materialIndex)->GetValue(
	      t.theHighestKineticEnergy,isOut);

  } else {
    
      dEdx = (*dEdxTable)(materialIndex)->GetValue(
                          scaledKineticEnergy,isOut) ;

  }

  return dEdx*Chargesquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 G4double G4EnergyLossTables::GetPreciseRangeFromEnergy(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial)
{
  if( aParticle != lastParticle)
  {
    t= GetTables(aParticle);
    lastParticle = aParticle;
    Chargesquare = (aParticle->GetPDGCharge())*
                   (aParticle->GetPDGCharge())/
                    QQPositron ;
  }
  const G4PhysicsTable* rangeTable= t.theRangeTable;
  const G4PhysicsTable*  dEdxTable= t.theDEDXTable;

  G4int materialIndex = aMaterial->GetIndex();

  G4double Thighr = t.theHighestKineticEnergy*t.theLowestKineticEnergy/
                   (*rangeTable)(materialIndex)->
                   GetLowEdgeEnergy(1) ;

  G4double scaledKineticEnergy = KineticEnergy*t.theMassRatio;
  G4double Range;
  G4bool isOut;

  if (scaledKineticEnergy<t.theLowestKineticEnergy) {

    Range = scaledKineticEnergy/t.theLowestKineticEnergy*
            (*rangeTable)(materialIndex)->GetValue(
              t.theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>Thighr) {

    Range = (*rangeTable)(materialIndex)->GetValue(
	      Thighr,isOut)+
            (scaledKineticEnergy-Thighr)/
            (*dEdxTable)(materialIndex)->GetValue(
              Thighr,isOut);

  } else {
    
     Range = (*rangeTable)(materialIndex)->GetValue(
                       scaledKineticEnergy,isOut) ;

  }

  return Range/(Chargesquare*t.theMassRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......





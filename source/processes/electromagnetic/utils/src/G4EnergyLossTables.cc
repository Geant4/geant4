// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//  
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnergyLossTables.cc,v 1.3 1999-02-17 08:51:43 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// first version created by P.Urban , 06/04/1998
// modifications + "precise" functions added by L.Urban , 27/05/98
// modifications , TOF functions , 26/10/98, L.Urban
// cache mechanism in order to gain time, 11/02/99, L.Urban

#include "G4EnergyLossTables.hh"

G4EnergyLossTablesHelper G4EnergyLossTables::t  ;
const G4ParticleDefinition* G4EnergyLossTables::lastParticle = NULL ; 
G4double G4EnergyLossTables::Chargesquare ;

RWTValHashDictionary<const G4ParticleDefinition*, G4EnergyLossTablesHelper>
G4EnergyLossTables::dict(G4EnergyLossTables::HashFun);

#if defined(GNU_GCC) && ((__GNUC__ * 10 + __GNUC_MINOR__) < 28)
template RWTValHashDictionary
  <G4ParticleDefinition const *, G4EnergyLossTablesHelper>
  ::RWTValHashDictionary(unsigned int (*)(G4ParticleDefinition
  const *const &), unsigned int);
#endif

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
  Chargesquare = (p->GetPDGCharge())*(p->GetPDGCharge());
}

///////////////////////////////////////////////////////////////////////

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
                   (aParticle->GetPDGCharge()) ;
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
                   (aParticle->GetPDGCharge()) ;
  }
  const G4PhysicsTable* rangeTable= t.theRangeTable;
  const G4PhysicsTable*  dEdxTable= t.theDEDXTable;

  G4int materialIndex = aMaterial->GetIndex();

  G4double Thigh = t.theHighestKineticEnergy*t.theLowestKineticEnergy/
                   (*rangeTable)(materialIndex)->
                   GetLowEdgeEnergy(1) ;

  G4double scaledKineticEnergy = KineticEnergy*t.theMassRatio;
  G4double Range;
  G4bool isOut;

  if (scaledKineticEnergy<t.theLowestKineticEnergy) {

    Range = scaledKineticEnergy/t.theLowestKineticEnergy*
            (*rangeTable)(materialIndex)->GetValue(
              t.theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>Thigh) {

    Range = (*rangeTable)(materialIndex)->GetValue(
	      Thigh,isOut)+
            (scaledKineticEnergy-Thigh)/
            (*dEdxTable)(materialIndex)->GetValue(
              Thigh,isOut);

  } else {
    
     Range = (*rangeTable)(materialIndex)->GetValue(
                       scaledKineticEnergy,isOut) ;

  }

  return Range/(Chargesquare*t.theMassRatio);
}





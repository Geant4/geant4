// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedBaryonConstructor.cc,v 1.1 1999-01-07 16:10:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
//      History: first implementation, based on object model of
//      10 oct 1998  H.Kurashige
// ---------------------------------------------------------------


#include "G4ExcitedBaryonConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"


G4ExcitedBaryonConstructor::G4ExcitedBaryonConstructor(G4int nStates,
						       G4int isoSpin):
         NumberOfStates(nStates),iIsoSpin(isoSpin),
         iConjugation(0), iGParity(0), leptonNumber(0), baryonNumber(1),
	 type("baryon")
{
 
}

G4ExcitedBaryonConstructor::~G4ExcitedBaryonConstructor()
{
}

void G4ExcitedBaryonConstructor::Construct(G4int idx)
{
  if (idx < 0 ) {
    for (G4int state=0; state< NumberOfStates; state +=1) {
       ConstructParticle(state);
       ConstructAntiParticle(state);
     }
  } else if (idx < NumberOfStates) {
    ConstructParticle(idx);
    ConstructAntiParticle(idx);
  } else {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>1) {
      G4cerr << "G4ExcitedBaryonConstructor::Construct()";
      G4cerr << "   illegal index os state = " << idx << endl;
    }
#endif
  }
}


#include "G4ExcitedBaryons.hh"

void G4ExcitedBaryonConstructor::ConstructParticle(G4int idx)
{
  if (!Exist(idx) ) return;

  //    Construct Resonace particles as dynamic object
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table 
  
  
  G4String name;
  G4ParticleDefinition* particle;
  
  for (G4int iIso3 = -1*iIsoSpin; iIso3 <= iIsoSpin; iIso3 +=2) {
    name= GetName(iIso3, idx);

    particle = new G4ExcitedBaryons(            
		 name,    GetMass(idx), GetWidth(idx),    GetCharge(iIso3), 
	GetiSpin(idx), GetiParity(idx),  iConjugation,       
	     iIsoSpin,           iIso3,      iGParity,
                 type,    leptonNumber,  baryonNumber, GetEncoding( iIso3,idx),
                false,             0.0,   NULL
				    );
    particle->SetDecayTable(CreateDecayTable( name, iIso3, idx, false));
  }
}

void G4ExcitedBaryonConstructor::ConstructAntiParticle(G4int idx)
{
  if (!Exist(idx) ) return;

  //    Construct Resonace particles as dynamic object
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table 
  
  
  G4String name;
  G4ParticleDefinition* particle;
  
  for (G4int iIso3 = -1*iIsoSpin; iIso3 <= iIsoSpin; iIso3 +=2) {
    name = GetName(iIso3, idx);
    name = "anti_" + name;

    particle = new G4ExcitedBaryons(            
		 name,    GetMass(idx), GetWidth(idx),  -1.0*GetCharge(iIso3), 
	GetiSpin(idx), GetiParity(idx),  iConjugation,       
	     iIsoSpin,        -1*iIso3,      iGParity,
                 type,    leptonNumber, 
				       -1*baryonNumber, 
				                   -1*GetEncoding( iIso3,idx),
                false,         0.0,   NULL
				    );
    particle->SetDecayTable(CreateDecayTable( name, iIso3, idx, true));
  }
   
}

G4double  G4ExcitedBaryonConstructor::GetCharge(G4int iIsoSpin3)
{
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4double charge = 0.0;
  static G4double quark_charge[7] = 
  {
    0., -1./3., +2./3., -1./3., +2./3., -1./3., +2./3.
  };
  
  for (G4int idx=0; idx<3; idx+=1){
    charge += quark_charge[GetQuarkContents(idx, iIsoSpin3)]*eplus;
  }
  return charge;
}

G4int     G4ExcitedBaryonConstructor::GetEncoding(G4int iIsoSpin3, G4int idxState)
{
  G4int encoding = GetEncodingOffset(idxState);
  encoding += 1000*GetQuarkContents(0, iIsoSpin3);
  encoding +=  100*GetQuarkContents(1, iIsoSpin3);
  encoding +=   10*GetQuarkContents(2, iIsoSpin3);
  if (GetiSpin(idxState) <9) {
    encoding += GetiSpin(idxState) +1;
  } else {
    encoding += (GetiSpin(idxState) +1)*10000000;
  }
  return encoding;
}




















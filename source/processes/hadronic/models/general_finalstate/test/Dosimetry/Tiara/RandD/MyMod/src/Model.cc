#include "INCModel.hh"
#include "NuclModel.h"
#include "Model.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "PreCompoundTransitions.hh"
#include "PreCompoundEmission.hh"
#include "G4LorentzVector.hh"
#include "G4PreCompoundEmission.hh"
#include "G4PreCompoundTransitions.hh"
#include "G4PreCompoundParameters.hh"


#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "PreCompoundParameters.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"

  Model::~Model() {}
  
struct PROJ{
  unsigned char type;
  double Energy;
  VECTOR Mom;
  VECTOR pos;
};

static PROJ Projectiles[300];
static unsigned nProj;
static unsigned n=20;
static unsigned nNucSl = 8;
extern const double g_dProtonMass = 988.;
extern const double g_dNeutronMass = 989.;

double Model::GetTau()
{
  if(nProj==0) return 0;
  double res=1e+12,res1;
  for(unsigned i=0;i<nProj;i++){
    res1 = m_INCModel.GetTau(&Projectiles[i],n);
    if(res1 < res) res = res1;
  }
  return res;
}

void Model::Emit(unsigned nPart)
{
  //Ako izleze izvyn jadroto - tja se slaga kato secondary
  if(nPart <nProj){
    float fTmp;
    if(Projectiles[nPart].type==PROTON) G4cout<<"Emiting proton";
    else G4cout<<"Emitin neutron";
    //    G4cout<<" energy: "<<Projectiles[nPart].Energy<<G4endl;
    G4ReactionProduct* pNewProduct = new G4ReactionProduct(((Projectiles[nPart].type==PROTON) ? ((G4VBaryon*)G4Proton::ProtonDefinition()) : ((G4VBaryon*)G4Neutron::NeutronDefinition())));
    //    fTmp = sqrt(Projectiles[nPart].Mom.x*Projectiles[nPart].Mom.x +
    //Projectiles[nPart].Mom.y*Projectiles[nPart].Mom.y+
    //	Projectiles[nPart].Mom.z*Projectiles[nPart].Mom.z);
    pNewProduct->SetMomentum(G4ThreeVector(Projectiles[nPart].Mom.x,Projectiles[nPart].Mom.y,Projectiles[nPart].Mom.z));
    fTmp = Projectiles[nPart].Energy;
    pNewProduct->SetKineticEnergy(fTmp*MeV);
    m_INCModel.AddExcitationEnergy(-pNewProduct->GetKineticEnergy());
    G4cout<<"Left excitation energy: "<<m_INCModel.GetExcitationEnergy()<<G4endl;
    if(Projectiles[nPart].type==PROTON){
      fTmp -= 1.44*m_INCModel.GetCurrZ()*eplus/MeV*eplus/MeV/coulomb/coulomb/(1.07*pow(m_INCModel.GetCurrA(),1./3.)+1.88);
    }
    fTmp -= G4NucleiProperties::GetBindingEnergy(m_INCModel.GetCurrA(),m_INCModel.GetCurrZ());
    G4cout<<" energy: "<<fTmp<<G4endl;
    pNewProduct->SetKineticEnergy(fTmp*MeV);
    m_pResultVector->push_back(pNewProduct);
    m_INCModel.AddRecoilMomentum(&Projectiles[nPart]);
    m_INCModel.SetCurrA(m_INCModel.GetCurrA()-1);
    if(Projectiles[nPart].type==PROTON){
      m_INCModel.SetCurrZ(m_INCModel.GetCurrZ()-1);
    }
    for(unsigned i=nPart;i<nProj-1;i++){
      memcpy(&Projectiles[i],&Projectiles[i+1],sizeof(PROJ));
    }
  }
  if(--nProj<0) nProj=0;
}

void Model::ScaleRadial(unsigned nPart,double RefrCoef)
{
  VECTOR vTmp,vTmp1;
  if(RefrCoef>10000) G4cout<<"RefrCoef: "<<RefrCoef<<G4endl;
  vTmp.x = Projectiles[nPart].pos.x;
  vTmp.y = Projectiles[nPart].pos.y;
  vTmp.z = Projectiles[nPart].pos.z;
  double fTmp = sqrt(vTmp.x*vTmp.x + vTmp.y*vTmp.y + vTmp.z*vTmp.z)/*,fTmp1*/;
  vTmp.x /= fTmp;
  vTmp.y /= fTmp;
  vTmp.z /= fTmp;
  fTmp = Projectiles[nPart].Mom.x*vTmp.x + Projectiles[nPart].Mom.y*vTmp.y+
	  Projectiles[nPart].Mom.z*vTmp.z;

  //vTmp1 - tangencialnata chast
  vTmp1.x = Projectiles[nPart].Mom.x - fTmp*vTmp.x;
  vTmp1.y = Projectiles[nPart].Mom.y - fTmp*vTmp.y;
  vTmp1.z = Projectiles[nPart].Mom.z - fTmp*vTmp.z;
  vTmp.x = fTmp*vTmp.x;
  vTmp.y = fTmp*vTmp.y;
  vTmp.z = fTmp*vTmp.z;
  vTmp.x *= RefrCoef;
  vTmp.y *= RefrCoef;
  vTmp.z *= RefrCoef;
  Projectiles[nPart].Mom.x = vTmp.x + vTmp1.x;
  Projectiles[nPart].Mom.y = vTmp.y + vTmp1.y;
  Projectiles[nPart].Mom.z = vTmp.z + vTmp1.z;
  fTmp = Projectiles[nPart].Mom.x*Projectiles[nPart].Mom.x+
    Projectiles[nPart].Mom.y*Projectiles[nPart].Mom.y+
    Projectiles[nPart].Mom.z*Projectiles[nPart].Mom.z;
  Projectiles[nPart].Energy = fTmp/2/((Projectiles[nPart].type==PROTON) ? g_dProtonMass : g_dNeutronMass);
}

double Model::GetPotentialDiff(unsigned nPart)
{
  unsigned nPos1,nPos2;
  double m_NuclRad = 1.07*pow(m_INCModel.GetCurrA(),1./3.)+1;
  double m_NuclThick = 0.54;
  double Rad = sqrt(Projectiles[nPart].pos.x*Projectiles[nPart].pos.x+
		    Projectiles[nPart].pos.y*Projectiles[nPart].pos.y+
		    Projectiles[nPart].pos.z*Projectiles[nPart].pos.z);
  nPos1 = (unsigned)(Rad*(nNucSl-0.5)/m_NuclRad);
  nPos2 = (unsigned)(Rad*(nNucSl+0.5)/m_NuclRad);
  double Pot1,Pot2;
  if(Projectiles[nPart].type == PROTON){
    Pot1 = -(52.+33.*(m_INCModel.GetCurrA()-2*m_INCModel.GetCurrZ())/m_INCModel.GetCurrA())/
      (1+exp((((double)nPos1-(double)nNucSl)/(double)nNucSl*m_NuclRad)/m_NuclThick));
    Pot1 += (m_INCModel.GetCurrZ()*(eplus/MeV*eplus/MeV)/coulomb/coulomb/2/m_NuclRad)*
      (3 - (double)nPos1*(double)nPos1/(double)nNucSl/(double)nNucSl);
    Pot2 = -(52.+33.*(m_INCModel.GetCurrA()-2*m_INCModel.GetCurrZ())/m_INCModel.GetCurrA())/
      (1+exp((((double)nPos2-(double)nNucSl)/(double)nNucSl*m_NuclRad)/m_NuclThick));
    Pot2 += (m_INCModel.GetCurrZ()*(eplus/MeV*eplus/MeV)/coulomb/coulomb/2/m_NuclRad)*
      (3-(double)nPos2*(double)nPos2/(double)nNucSl/(double)nNucSl);
  }
  else{
    Pot1 = -(52.-33.*(m_INCModel.GetCurrA()-2*m_INCModel.GetCurrZ())/m_INCModel.GetCurrA());
    Pot2 = Pot1/(1+exp((((double)nPos2-(double)nNucSl)/(double)nNucSl*m_NuclRad)/m_NuclThick));
    Pot1 /= 1+exp((((double)nPos1-(double)n)*m_NuclRad)/m_NuclThick);
  }
  //sega: Pot1 - Pot2 ili Pot2-Pot1;
  float fTmp;
  if(nPos1<nPos2){
    fTmp = Projectiles[nPart].pos.x*Projectiles[nPart].Mom.x+
      Projectiles[nPart].pos.y*Projectiles[nPart].Mom.y+
      Projectiles[nPart].pos.z*Projectiles[nPart].Mom.z;
    if(fTmp<0) return Pot1-Pot2;
    return Pot2-Pot1;
  }
  fTmp = Projectiles[nPart].pos.x*Projectiles[nPart].Mom.x+
    Projectiles[nPart].pos.y*Projectiles[nPart].Mom.y+
    Projectiles[nPart].pos.z*Projectiles[nPart].Mom.z;
  if(fTmp<0) return Pot2-Pot1;
  return Pot1-Pot2; 
}

double Model::GetRefrCoef(unsigned nPart)
{
  double p,En,fTmp;
  p = Projectiles[nPart].Mom.x*Projectiles[nPart].pos.x+
    Projectiles[nPart].Mom.y*Projectiles[nPart].pos.y+
    Projectiles[nPart].Mom.z*Projectiles[nPart].pos.z;
  fTmp = sqrt(Projectiles[nPart].pos.x*Projectiles[nPart].pos.x+Projectiles[nPart].pos.y*Projectiles[nPart].pos.y+
	      Projectiles[nPart].pos.z*Projectiles[nPart].pos.z);
  p /= fTmp;
  //p *= p;
  En = 2*((Projectiles[nPart].type==PROTON) ? g_dProtonMass : g_dNeutronMass)*(fTmp=GetPotentialDiff(nPart));
  if(fTmp>10000) G4cout<<"PotDiff: "<<fTmp<<G4endl;
  En = 1 + En/p;
  if(En < 0) return 1;
  else return sqrt(1-En/p);
}

bool Model::Propagate(unsigned nPart,double tau)
{
  double RealPath=tau;
  if(m_INCModel.CheckCrossing(&Projectiles[nPart],nNucSl,RealPath)){
    if(m_INCModel.Interact(&Projectiles[nPart],&Projectiles[nProj],RealPath)) nProj++;
    //ScaleRadial(nPart,GetRefrCoef(nPart));
    /*    Projectiles[nPart].Energy += GetPotentialDiff(nPart);
    RealPath = sqrt(2*Projectiles[nPart].Energy*
    (Projectiles[nPart].type==PROTON) ? g_dProtonMass : g_dNeutronMass);*/
  }
  else{
    if(m_INCModel.Interact(&Projectiles[nPart],&Projectiles[nProj],RealPath)) nProj++;
  }
  if(m_INCModel.IsOutside(&Projectiles[nPart])) return true;
  return false;
}

void Model::MainCycle()
{
  double currTau;
  unsigned i,j;
  do{
    currTau = GetTau();
    for(i=0;i<nProj;i++){
      if(Propagate(i,currTau)){
	Emit(i);
      }
    }
    for(i=0;i<nProj;i++){
      if(Projectiles[i].Energy < m_INCModel.GetMinimumProjectileEnergy((Projectiles[i].type==PROTON)?1:0,1,Projectiles[i].pos)){
	m_INCModel.AddParticle();
	m_INCModel.AddExcitationEnergy(Projectiles[i].Energy);
	m_INCModel.AddRecoilMomentum(&Projectiles[i]);
	if(Projectiles[i].type==PROTON) m_INCModel.SetCurrZ(m_INCModel.GetCurrZ()+1);
	m_INCModel.SetCurrA(m_INCModel.GetCurrA()+1);
	for(j=i;j<nProj-1;j++){
	  memcpy(&Projectiles[j],&Projectiles[j+1],sizeof(PROJ));
	}
	nProj--;
      }
    }
  }
  while(nProj != 0);
}

void Model::PerformPreequilibrium()
{
  G4Fragment m_CurrFragm;
  m_CurrFragm.SetA(m_INCModel.GetCurrA());
  m_CurrFragm.SetZ(m_INCModel.GetCurrZ());
  m_CurrFragm.SetNumberOfParticles(m_INCModel.GetParticles());
  m_CurrFragm.SetNumberOfCharged(m_INCModel.GetCharged());
  m_CurrFragm.SetNumberOfHoles(m_INCModel.GetHoles());
  G4cout<<"Preequilibrium: excitationEnergy: "<<m_INCModel.GetExcitationEnergy()<<G4endl;
  //energija i impuls na jadroto
  //  VECTOR vMom1;
  //  double En = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(m_INCModel.GetCurrZ(),m_INCModel.GetCurrA())+ m_INCModel.GetExcitationEnergy()*MeV;
  // En *= En;
  //En += vMom.x*vMom.x+vMom.y*vMom.y+vMom.z*vMom.z;
  G4double En = m_INCModel.GetExcitationEnergy() +
    G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(m_INCModel.GetCurrZ(),m_INCModel.GetCurrA());
  En *= En;
  VECTOR vMom;
  m_INCModel.GetRestMomentum(vMom);
  //  m_CurrFragm.SetExcitationEnergy(m_INCModel.GetExcitationEnergy()*MeV);
  //  m_CurrFragm.SetMomentum(G4LorentzVector(G4ThreeVector(vMom.x,vMom.y,vMom.z),sqrt(En+G4ThreeVector(vMom.x,vMom.y,vMom.z).mag2())));
  m_CurrFragm.SetExcitationEnergy(m_INCModel.GetExcitationEnergy()*MeV);
  G4ReactionProductVector * result = DeExciteEq(m_CurrFragm);
  //  m_INCModel.GetRestMomentum(vMom1);
  for(G4ReactionProductVector::iterator i= result->begin(); i != result->end(); ++i){
    if((*i)==NULL) continue;
    /*    vMom.x = (*i)->GetMomentum().x()/MeV;
	  vMom.y = (*i)->GetMomentum().y()/MeV;
	  vMom.z = (*i)->GetMomentum().z()/MeV;
	  vMom.x += vMom1.x;
	  vMom.y += vMom1.y;
	  vMom.z += vMom1.z;
	  En = (*i)->GetDefinition()->GetPDGMass()/MeV;
	  En*=En;
	  En += vMom.x*vMom.x+vMom.y*vMom.y+vMom.z*vMom.z;
	  (*i)->SetMomentum(G4LorentzVector(G4ThreeVector(vMom.x,vMom.y,vMom.z),sqrt(En)));
	  m_INCModel.AddRecoilMomentum(vMom);*/
    m_pResultVector->push_back(*i);
  }
  delete result;
}

G4ReactionProductVector* Model::DeExciteEq(G4Fragment& fNuc)
{
  G4ReactionProductVector * Result = new G4ReactionProductVector;
  if (fNuc.GetA() < 5) {
    G4ReactionProduct * theRP = new G4ReactionProduct(G4ParticleTable::GetParticleTable()->
                                                      GetIon(fNuc.GetZ(),fNuc.GetA(),
                                                             fNuc.GetExcitationEnergy()));
    theRP->SetMomentum(fNuc.GetMomentum().vect());
    theRP->SetTotalEnergy(fNuc.GetMomentum().e());   
    Result->push_back(theRP);
    return Result;
  }
  if(fNuc.GetExcitationEnergy()<100/*m_INCModel.GetExcitationEnergy()<100*/){
    PreCompoundEmission aEmission;
    PreCompoundTransitions aTransition;
    G4double a;
    G4int i=0;
    for (;;i=0) {
      aEmission.Initialize(fNuc);
      //   G4double g = 0.595*fNuc.GetA()*(a=PreCompoundParameters::GetAddress()->GetLevelDensity(fNuc));
      G4double g = 1.2158542*(a=PreCompoundParameters::GetAddress()->GetLevelDensity(fNuc));
      G4int EquilibriumExcitonNumber = G4int(sqrt(2.0*g*fNuc.GetExcitationEnergy())+0.5);
      G4bool ThereIsTransition = false;
      do{
	i++;
	if(fNuc.GetNumberOfExcitons() < EquilibriumExcitonNumber && fNuc.GetA() > 4 && i<100000){
	  G4double TotalEmissionProbability = aEmission.GetTotalProbability(fNuc,a);
	  if (fNuc.GetNumberOfExcitons() <= 0){
	    return Result;
	  }
	  aTransition.Init(fNuc,a);
	  G4double TotalTransitionProbability = aTransition.GetTotalProbability();
	  if (G4UniformRand() > TotalEmissionProbability/(TotalTransitionProbability+TotalEmissionProbability)){
	    ThereIsTransition = true;
	    aTransition.PerformTransition(fNuc,(TotalEmissionProbability==0));
	  }
	  else{
	    ThereIsTransition = false;
	    G4ReactionProduct* vProd = aEmission.PerformEmission(fNuc,a);
	    /*	    G4ThreeVector vMom(vProd->GetMomentum().x(),vProd->GetMomentum().y(),vProd->GetMomentum().z());
		    VECTOR vMom1;
		    m_INCModel.GetRestMomentum(vMom1);
		    vMom = G4ThreeVector(vMom.x()+vMom1.x,vMom.y()+vMom1.y,
		    vMom.z()+vMom1.z);
		    G4double En = vProd->GetDefinition()->GetPDGMass();
		    En *= En;
		    En += vMom.mag2();
		    vProd->SetMomentum(vMom);
		    vProd->SetKineticEnergy(sqrt(En)-vProd->GetDefinition()->GetPDGMass());*/
	    Result->push_back(vProd);
	  }
	}
	else{
	  DeExcite(fNuc);
	  return Result;
	}
      }
      while(ThereIsTransition);
    }
  }
  else{
    G4double a;
    G4int i=0;
    for (;;i=0) {
      G4PreCompoundEmission aEmission(fNuc);
      aEmission.Initialize(fNuc);
      //      G4double g = 0.595*fNuc.GetA()*G4PreCompoundParameters::GetAddress()->GetLevelDensity();
      G4double g = 1.2158542*(a=PreCompoundParameters::GetAddress()->GetLevelDensity(fNuc));
      G4int EquilibriumExcitonNumber = G4int(sqrt(2.0*g*fNuc.GetExcitationEnergy())+0.5);
      G4bool ThereIsTransition = false;
      do{
	i++;
	if(fNuc.GetNumberOfExcitons() < EquilibriumExcitonNumber && fNuc.GetA() > 4 && i<100000){
	  G4double TotalEmissionProbability = aEmission.GetTotalProbability(fNuc);
	  if (fNuc.GetNumberOfExcitons() <= 0){
	    return Result;
	  }
	  G4PreCompoundTransitions aTransition(fNuc);
	  G4double TotalTransitionProbability = aTransition.GetTotalProbability();
	  if (G4UniformRand() > TotalEmissionProbability/(TotalTransitionProbability+TotalEmissionProbability)){
	    ThereIsTransition = true;
	    fNuc = aTransition.PerformTransition(fNuc);
	  }
	  else{
	    ThereIsTransition = false;
	    G4ReactionProduct* vProd = aEmission.PerformEmission(fNuc);
	    /*	    G4ThreeVector vMom(vProd->GetMomentum().x(),vProd->GetMomentum().y(),
			       vProd->GetMomentum().z());
			       VECTOR vMom1;
			       m_INCModel.GetRestMomentum(vMom1);
			       vMom = G4ThreeVector(vMom.x()+vMom1.x,vMom.y()+vMom1.y,
			       vMom.z()+vMom1.z);
			       G4double En = vProd->GetDefinition()->GetPDGMass();
			       En*= En;
			       En += vMom.mag2();
			       vProd->SetMomentum(vMom);
			       vProd->SetKineticEnergy(sqrt(En)-vProd->GetDefinition()->GetPDGMass());*/
	    Result->push_back(vProd);
	  }
	}
	else{
	  DeExcite(fNuc);
	  return Result;
	}
      }
      while(ThereIsTransition);
    }
  }
}

void Model::DeExcite(G4Fragment& fNuc)
{
  VECTOR vMom1,vMom;
  double En;
  if(fNuc.GetA()<5){
    G4ReactionProductVector* pResult = new G4ReactionProductVector();
    G4ReactionProduct * theRP = new G4ReactionProduct(G4ParticleTable::GetParticleTable()->
                                                      GetIon(fNuc.GetZ(),fNuc.GetA(),
                                                             fNuc.GetExcitationEnergy()));
   theRP->SetMomentum(fNuc.GetMomentum().vect());
   theRP->SetTotalEnergy(fNuc.GetMomentum().e());   
   m_pResultVector->push_back(theRP);
   //   return pResult;
  }
  G4LorentzVector pRestVector(fNuc.GetMomentum());
  G4ThreeVector pBoostVector = pRestVector.boostVector();
  pRestVector = pRestVector.rest4Vector();
  fNuc.SetMomentum(pRestVector);
  G4ReactionProductVector* pVect = m_pExcitation->BreakItUp(fNuc);
  //m_INCModel.GetRestMomentum(vMom1);
  G4ReactionProductVector::iterator i;
  for(i=pVect->begin();i!=pVect->end();i++){
    /*    vMom.x = (*i)->GetMomentum().x()/MeV;
	  vMom.y = (*i)->GetMomentum().y()/MeV;
	  vMom.z = (*i)->GetMomentum().z()/MeV;
	  vMom.x += vMom1.x;
	  vMom.y += vMom1.y;
	  vMom.z += vMom1.z;
	  En = (*i)->GetDefinition()->GetPDGMass()/MeV;
	  En *= En;
	  En += vMom.x*vMom.x+vMom.y*vMom.y+vMom.z*vMom.z;
	  (*i)->SetMomentum(G4ThreeVector(vMom.x,vMom.y,vMom.z));
	  (*i)->SetKineticEnergy(sqrt(En)-(*i)->GetDefinition()->GetPDGMass());*/
    pRestVector = G4LorentzVector((*i)->GetMomentum(),(*i)->GetTotalEnergy());
    pRestVector.boost(pBoostVector);
    (*i)->SetMomentum(pRestVector.vect());
    (*i)->SetTotalEnergy(pRestVector.mag());
    m_pResultVector->push_back(*i);
  }
  delete pVect;
}

G4VParticleChange* Model::ApplyYourself(const G4Track& pTrack,G4Nucleus& pNucl)
{
  /*
    m_pResultVector = new G4ReactionProductVector();
    G4ParticleChange* theResult = new G4ParticleChange;
    theResult->Initialize(pTrack);
    nProj=1;
    G4double fTmp;
    m_INCModel.SetNucleus(&pNucl);
    if(pTrack.GetDefinition() == G4Proton::ProtonDefinition()){
    Projectiles[0].type = PROTON;
    Projectiles[0].Energy = pTrack.GetDynamicParticle()->GetKineticEnergy()/MeV;
    Projectiles[0].Mom.y = Projectiles[0].Mom.z = Projectiles[0].Mom.x = sqrt(2*Projectiles[0].Energy*g_dProtonMass);//sqrt((fTmp = pTrack.GetDynamicParticle()->GetTotalEnergy()/MeV)*fTmp - g_dProtonMass*g_dProtonMass);
    m_INCModel.SetCurrZ(m_INCModel.GetCurrZ()+1);
    }
    else{
    Projectiles[0].type = NEUTRON;
    Projectiles[0].Energy = pTrack.GetDynamicParticle()->GetKineticEnergy()/MeV;
    Projectiles[0].Mom.y = Projectiles[0].Mom.z = Projectiles[0].Mom.x = sqrt(2*Projectiles[0].Energy*g_dNeutronMass);//sqrt((fTmp=pTrack.GetDynamicParticle()->GetTotalEnergy()/MeV)*fTmp - g_dNeutronMass*g_dNeutronMass);
    }
    G4cout<<"Projectiles energy: "<<Projectiles[0].Energy<<G4endl;
    m_INCModel.SetCurrA(m_INCModel.GetCurrA()+1);
    Projectiles[0].Mom.x *= pTrack.GetDynamicParticle()->GetMomentumDirection().x();
    Projectiles[0].Mom.y *= pTrack.GetDynamicParticle()->GetMomentumDirection().y();
    Projectiles[0].Mom.z *= pTrack.GetDynamicParticle()->GetMomentumDirection().z();
    fTmp = sqrt(Projectiles[0].Mom.x*Projectiles[0].Mom.x+
    Projectiles[0].Mom.y*Projectiles[0].Mom.y+
    Projectiles[0].Mom.z*Projectiles[0].Mom.z);
    fTmp = (1.07*pow(pNucl.GetN(),1./3.))/fTmp;
    Projectiles[0].pos.x = -Projectiles[0].Mom.x*fTmp;
    Projectiles[0].pos.y = -Projectiles[0].Mom.y*fTmp;
    Projectiles[0].pos.z = -Projectiles[0].Mom.z*fTmp;
    MainCycle();
    PerformPreequilibrium();*/
  //Second Try. Samo precompound modela

  m_pResultVector = new G4ReactionProductVector();
  static G4ParticleChange theResult;
  theResult.Initialize(pTrack);

  G4Fragment fNuc;
  fNuc.SetA(pNucl.GetN()/*+pTrack.GetDynamicParticle()->GetDefinition()->
			  GetBaryonNumber()*/);

  fNuc.SetZ(pNucl.GetZ()/*+pTrack.GetDynamicParticle()->GetDefinition()->
			  GetPDGCharge()*/);

  fNuc.SetNumberOfParticles(1+pTrack.GetDynamicParticle()->GetDefinition()->
			    GetBaryonNumber());
  fNuc.SetNumberOfCharged(pTrack.GetDynamicParticle()->GetDefinition()->
			  GetPDGCharge()+(int)(G4UniformRand()+0.5));
  fNuc.SetNumberOfHoles(1);
  G4double En = fNuc.GetGroundStateMass()/*G4NucleiProperties::GetNuclearMass(pNucl.GetN(),pNucl.GetZ())*/+
    pTrack.GetDynamicParticle()->GetTotalEnergy();
  fNuc.SetA(fNuc.GetA()+pTrack.GetDynamicParticle()->GetDefinition()->
	    GetBaryonNumber());
  fNuc.SetZ(fNuc.GetZ()+pTrack.GetDynamicParticle()->GetDefinition()->
	    GetPDGCharge());
  G4ThreeVector p = pTrack.GetDynamicParticle()->Get4Momentum().vect();
  G4LorentzVector pRestVector = G4LorentzVector(p,En);
  G4ThreeVector pBoostVector = pRestVector.boostVector();
  pRestVector = pRestVector.rest4Vector();
  fNuc.SetMomentum(pRestVector);
  G4ReactionProductVector* pVec = DeExciteEq(fNuc);
  G4ReactionProductVector::iterator i;
  for(i=pVec->begin();i!=pVec->end();++i)
    m_pResultVector->push_back(*i);
  delete pVec;

  theResult.SetStatusChange(fStopAndKill);
  theResult.SetNumberOfSecondaries(m_pResultVector->size());
  for(i= m_pResultVector->begin(); i != m_pResultVector->end(); ++i){
    pRestVector = G4LorentzVector((*i)->GetMomentum(),(*i)->GetTotalEnergy());
    pRestVector = pRestVector.boost(pBoostVector);
    G4DynamicParticle * aNew = new G4DynamicParticle((*i)->GetDefinition(),pRestVector);
    theResult.AddSecondary(aNew);
    if(aNew->GetDefinition()->GetBaryonNumber() != 0)
      G4cout<<" emited "<< aNew->GetDefinition()->GetPDGCharge()<<" "
	    <<aNew->GetDefinition()->GetBaryonNumber()<<" "
	    <<aNew->GetKineticEnergy()<<G4endl;
    delete (*i);
  }
  delete m_pResultVector;
  return &theResult;
}

#include "NuclModel.h"
#include "INCModel.hh"
#include "G4Fragment.hh"
#include "G4Nucleus.hh"
#include "Randomize.hh"
#include "G4NucleiProperties.hh"
#include "G4FermiMomentum.hh"

struct PROJ{
  unsigned char type;
  double Energy;
  VECTOR Mom;
  VECTOR pos;
};

const double g_dProtonMass = 988;
const double g_dNeutronMass = 989;

INCModel::INCModel(G4Nucleus* pNucl)
{
  m_cParticles = m_cCharged = m_cHoles = 0;
  if(pNucl!=NULL)
    m_pNucleus = new G4Nucleus(*pNucl);
  else
    m_pNucleus = NULL;
}

void INCModel::SetNucleus(G4Nucleus* pNuc){
  m_cParticles = m_cCharged = m_cHoles = 0;
  if(m_pNucleus!=NULL) delete m_pNucleus;
  m_pNucleus = new G4Nucleus(*pNuc);
}
bool INCModel::Interact(PROJ* pProj,PROJ* pNewProj,double tau)
{
  //Opredelja dali shte vzaimodejstva za vreme tau i vzaimodejstva ako she go pravi
  //Pyrvo : izbor na vzaimodejstvashta chastica
  PROJ SecProj;
  float fTmp1,fTmp2,fTmp;
  float fMass1,fMass2;
 loop:
  if(G4UniformRand() < GetCurrZ()/GetCurrA()) SecProj.type = PROTON;
  else SecProj.type = NEUTRON;
  fTmp1 = pi*G4UniformRand();
  fTmp2 = 2*pi*G4UniformRand();
  SecProj.Mom.x = sin(fTmp1)*sin(fTmp2);
  SecProj.Mom.y = sin(fTmp1)*cos(fTmp2);
  SecProj.Mom.z = cos(fTmp1);
  SecProj.pos.x = pProj->pos.x;
  SecProj.pos.y = pProj->pos.y;
  SecProj.pos.z = pProj->pos.z;
  //sega skorostta i energijata;
  fTmp = sqrt(SecProj.pos.x*SecProj.pos.x + SecProj.pos.y*SecProj.pos.y + SecProj.pos.z*SecProj.pos.z);
  SecProj.Energy = sqrt(GetPotential(fTmp,GetCurrA(),GetCurrZ())/(2*(fMass2 = (SecProj.type==PROTON)? g_dProtonMass : g_dNeutronMass)));
  fTmp1 = sqrt(SecProj.Energy*2*fMass2);
  SecProj.Mom.x *= fTmp1;
  SecProj.Mom.y *= fTmp1;
  SecProj.Mom.z *= fTmp1;
  fMass1 = (pProj->type==PROTON) ? g_dProtonMass : g_dNeutronMass;
  //Sega opredeljame energijata v centyra na masite
  double Ecm = pProj->Mom.x*pProj->Mom.x+pProj->Mom.y*pProj->Mom.y+pProj->Mom.z*pProj->Mom.z;
  Ecm += SecProj.Mom.x*SecProj.Mom.x + SecProj.Mom.y*SecProj.Mom.y + SecProj.Mom.z*SecProj.Mom.z;
  Ecm += 2*(pProj->Mom.x*SecProj.Mom.x + pProj->Mom.y*SecProj.Mom.y + pProj->Mom.z*SecProj.Mom.z);
  Ecm /= 2*(fMass1+fMass2);

  //Trjabva ni plytnostta i sigma;
  double rhosigma = 0.19487677;
  double Elab = (pProj->Mom.x - fMass1/fMass2*SecProj.Mom.x)*(pProj->Mom.x - fMass1/fMass2*SecProj.Mom.x);
  Elab += (pProj->Mom.y - fMass1/fMass2*SecProj.Mom.y)*(pProj->Mom.y - fMass1/fMass2*SecProj.Mom.y);
  Elab += (pProj->Mom.z - fMass1/fMass2*SecProj.Mom.z)*(pProj->Mom.z - fMass1/fMass2*SecProj.Mom.z);
  Elab /= 2*fMass1;
  fTmp = (fMass1/fMass2*SecProj.Mom.x-pProj->Mom.x)*(fMass1/fMass2*SecProj.Mom.x-pProj->Mom.x);
  fTmp += (fMass1/fMass2*SecProj.Mom.y-pProj->Mom.y)*(fMass1/fMass2*SecProj.Mom.y-pProj->Mom.y);
  fTmp += (fMass1/fMass2*SecProj.Mom.z-pProj->Mom.z)*(fMass1/fMass2*SecProj.Mom.z-pProj->Mom.z);
  fTmp = sqrt(fTmp);
  fTmp /= sqrt(pProj->Mom.x*pProj->Mom.x + pProj->Mom.y*pProj->Mom.y + pProj->Mom.z*pProj->Mom.z);
  rhosigma *= fTmp*GetXSection(Elab,(pProj->type==SecProj.type)?EQUAL_TYPE : DIFF_TYPE)*fermi*fermi/millibarn;

  //Sega ostava da smetnem N(delta(a));
  fTmp1 = sqrt(pProj->Mom.x*pProj->Mom.x + pProj->Mom.y*pProj->Mom.y + pProj->Mom.z*pProj->Mom.z);
  fTmp = fTmp1/fMass1;
  double a = fTmp*tau;
  double path = G4UniformRand();
  if(path <= rhosigma*a){
    //interaction
    path = path/rhosigma;
    //move particle path units forward;
    SecProj.pos.x = pProj->pos.x += pProj->Mom.x/fTmp1*path;
    SecProj.pos.y = pProj->pos.y += pProj->Mom.y/fTmp1*path;
    SecProj.pos.z = pProj->pos.z += pProj->Mom.z/fTmp1*path;
    //Izberi impuls i smetni energijata
    VECTOR vSumOfMom;
    vSumOfMom.x = pProj->Mom.x + SecProj.Mom.x;
    vSumOfMom.y = pProj->Mom.y + SecProj.Mom.y;
    vSumOfMom.z = pProj->Mom.z + SecProj.Mom.z;
    double costheta,phi = 2*pi*G4UniformRand();
    costheta = GetAngle(pProj->Energy,fMass1,SecProj.Energy,fMass2);
    G4cout<<flush;
    VECTOR vNewMom;
    vNewMom.x = pProj->Mom.x - SecProj.Mom.x*fMass1/fMass2;
    vNewMom.y = pProj->Mom.y - SecProj.Mom.y*fMass1/fMass2;
    vNewMom.z = pProj->Mom.z - SecProj.Mom.z*fMass1/fMass2;
    fTmp1 = sqrt(vNewMom.x*vNewMom.x + vNewMom.y*vNewMom.y + vNewMom.z*vNewMom.z);
    vNewMom.x /= fTmp1;
    vNewMom.y /= fTmp1;
    vNewMom.z /= fTmp1;
    ModifyMomentum(costheta,phi,vNewMom);
    double dNewEn;
    dNewEn = 2*fTmp1*fTmp1*costheta*costheta/fMass1/((1+fMass2/fMass1)*(1+fMass2/fMass1));
    if(dNewEn>100000){
      G4cout<<"a: "<<a<<"  path: "<<path<<" Momentum "<<vNewMom.x<<","<<vNewMom.y<<","<<vNewMom.z<<G4endl;
      G4cout<<"costheta: "<<costheta<<" rhosigma "<<rhosigma<<G4endl;
      G4cout<<"Energy1: "<<pProj->Energy<<" energy2: "<<SecProj.Energy<<G4endl;
    }
    if(dNewEn<GetFermiEnergy(pProj->pos) || (pProj->Energy + SecProj.Energy - dNewEn < GetFermiEnergy(SecProj.pos))){
      //zabraneno samo momenta se promenja. No ne za dylgo
      // sled malko she byde probvaj pak
      fTmp = sqrt(pProj->Mom.x*pProj->Mom.x+pProj->Mom.y*pProj->Mom.y+pProj->Mom.z*pProj->Mom.z);
      /*      pProj->pos.x = pProj->Mom.x/fTmp*path;
	      pProj->pos.y = pProj->Mom.y/fTmp*path;
	      pProj->pos.z = pProj->Mom.z/fTmp*path;*/
      fTmp /= fMass1;
      tau -= (a-path)/fTmp;
      if(tau>0)
	goto loop;
      return false;
    }
    else{
      //Promenja se i energijata
      double dNewEn1 = pProj->Energy + SecProj.Energy - dNewEn;
      fTmp1 = sqrt(2*fMass1*dNewEn);
      vNewMom.x *= fTmp1;
      vNewMom.y *= fTmp1;
      vNewMom.z *= fTmp1;
      //Mestim go kakto treve
      vNewMom.x += fMass1/fMass2*SecProj.Mom.x;
      vNewMom.y += fMass1/fMass2*SecProj.Mom.y;
      vNewMom.z += fMass1/fMass2*SecProj.Mom.z;
      if(dNewEn1 < GetMinimumProjectileEnergy((SecProj.type==PROTON)?1:0,1,SecProj.pos)){
	//obrazuva se dvojka chastica dupka. Energijata na vyzbujdane se uvelichava
	pProj->Mom.x = vNewMom.x;
	pProj->Mom.y = vNewMom.y;
	pProj->Mom.z = vNewMom.z;
	pProj->Energy = dNewEn;
	m_cParticles++;
	if(SecProj.type==PROTON) m_cCharged++;
	m_cHoles++;
	m_pNucleus->AddExcitationEnergy(dNewEn1);
	G4cout<<"Particle/hole pair produced "<<m_cParticles<<" particles and "<<m_cHoles<<" holes"<<" excitation energy: "<<m_pNucleus->GetEnergyDeposit()<<G4endl;
	return false;
      }
      //obrazuva se nova chastica
      m_cHoles++;
      pNewProj->Mom.x = vSumOfMom.x - vNewMom.x;
      pNewProj->Mom.y = vSumOfMom.y - vNewMom.y;
      pNewProj->Mom.z = vSumOfMom.z - vNewMom.z;
      pNewProj->Energy = dNewEn1;
      pNewProj->type = SecProj.type;
      pNewProj->pos.x = SecProj.pos.x;
      pNewProj->pos.y = SecProj.pos.y;
      pNewProj->pos.z = SecProj.pos.z;
      pProj->Mom.x = vNewMom.x;
      pProj->Mom.y = vNewMom.y;
      pProj->Mom.z = vNewMom.z;
      pProj->Energy = dNewEn;
      //SetCurrA(GetCurrA()-1);
      //if(pNewProj->type==PROTON) SetCurrZ(GetCurrZ()-1);
      G4cout<<"Current number of particles: "<<m_cParticles<<" and "<<m_cCharged<<" are charged"<<G4endl;
      G4cout<<"Current number of holes: "<<m_cHoles<<" and the excitation energy "<<m_pNucleus->GetEnergyDeposit()<<G4endl;
      return true;
    } 
  }
  //njama interaction;
  pProj->pos.x += pProj->Mom.x/fTmp1*a;
  pProj->pos.y += pProj->Mom.y/fTmp1*a;
  pProj->pos.z += pProj->Mom.z/fTmp1*a;
  return false;
}

double INCModel::GetFermiEnergy(VECTOR& v)
{
  /*  double fTmp = sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
      fTmp = GetDensity(NormalizeDensity(m_dCurrA),fTmp,m_dCurrA);
      G4FermiMomentum fMom;
      fMom.Init(m_dCurrA,m_dCurrZ);
      return sqrt(fMom.GetFermiMomentum(fTmp)/MeV*fMom.GetFermiMomentum(fTmp)/MeV + (fTmp=G4NucleiProperties::GetNuclearMass(m_dCurrA,m_dCurrZ))*fTmp);
  */
  double fTmp = 35; //sqrt(m_pNucleus->GetFermiMomentum()/MeV*m_pNucleus->GetFermiMomentum()/MeV + 988*988);
  return fTmp;
}

double INCModel::GetMinimumProjectileEnergy(unsigned nZ,unsigned nA,VECTOR& vPos)
{
  double ret;
  /*  ret = G4NucleiProperties::GetNuclearMass(m_pNucleus->GetA(),m_pNucleus->GetZ());
      ret -= G4NucleiProperties::GetNuclearMass(m_pNucleus->GetA()-nA,m_pNucleus->GetZ()-nZ);*/
  ret = G4NucleiProperties::GetBindingEnergy(GetCurrA(),GetCurrZ());
  if(ret<0) ret = 0;
  ret += GetPotential(sqrt(vPos.x*vPos.x+vPos.y*vPos.y+vPos.z*vPos.z),nA,nZ);
  if(nZ > 0)
    ret += 1.44*GetCurrZ()*nZ*eplus/MeV*eplus/MeV/coulomb/coulomb/(1.07*pow(GetCurrA(),1./3.)+1.88);
  // G4cout<<" Minimum energy: "<<ret<<G4endl;
  return ret;
}

bool INCModel::CheckCrossing(PROJ* pProj,unsigned nSlices,double& tau)
{
  unsigned nSlice1,nSlice2;
  double Rad = sqrt(pProj->pos.x*pProj->pos.x + pProj->pos.y*pProj->pos.y+pProj->pos.z*pProj->pos.z);
  double ARad = 1.07*pow(GetCurrA(),1./3.)+1;
  double RadS,Discr,T1,T2,A,AB,B,fMass;
  fMass = (pProj->type==PROTON) ? g_dProtonMass : g_dNeutronMass;
  nSlice1 = (unsigned)(Rad/ARad*(double)nSlices);
  A = pProj->pos.x*pProj->pos.x + pProj->pos.y*pProj->pos.y + pProj->pos.z*pProj->pos.z; 
  B = pProj->Mom.x*pProj->Mom.x + pProj->Mom.y*pProj->Mom.y + pProj->Mom.z*pProj->Mom.z;
  B/= fMass*fMass;
  AB = pProj->pos.x*pProj->Mom.x + pProj->pos.y*pProj->Mom.y + pProj->pos.z*pProj->Mom.z; 
  AB /= fMass;
  if(nSlice1 !=0){
    nSlice2 = nSlice1-1;
    RadS = ARad*(double)nSlice2/(double)nSlices;
    Discr = 4*(AB*AB - (A-RadS*RadS)*B);
    if(Discr>=0){
      if(Discr==0){
	T1 = -AB/B;
	if(T1>=0 && T1<tau){
	  tau = T1;
	  return true;
	}
      }
      else{
	T1 = (-2*AB + sqrt(Discr))/2*B;
	T2 = (-2*AB - sqrt(Discr))/2*B;
	if(T1 > T2 && T2 > 0) T1 = T2;
	if(T1 > 0 && T1<tau){
	  tau = T1;
	  return true;
	}
      }
    }
  }
  if(nSlice1< nSlices){
    nSlice2 = nSlice1+1;
    RadS = ARad*(double)nSlice2/(double)nSlices;
    if(fabs(RadS-Rad)<1e-10 && nSlice2<nSlices){
      nSlice2++;
      RadS += ARad/(double)nSlices;
    }
    Discr = 4*(AB*AB - (A-RadS*RadS)*B);
    if(Discr>=0){
      if(Discr==0){
	T1 = -AB/B;
	if(T1>=0 && T1<tau){
	  tau = T1;
	  return true;
	}
      }
      else{
	T1 = (-2*AB + sqrt(Discr))/2*B;
	T2 = (-2*AB - sqrt(Discr))/2*B;
	if(T1 > T2 && T2 > 0) T1 = T2;
	if(T1 > 0 && T1<tau){
	  tau = T1;
	  return true;
	}
      }
    }
  }
  return false;
}

bool INCModel::IsOutside(PROJ* pProj)
{
  double Rad = sqrt(pProj->pos.x*pProj->pos.x + pProj->pos.y*pProj->pos.y+pProj->pos.z*pProj->pos.z);
  if(fabs(Rad-1.07*pow(GetCurrA(),1./3.)-1)<1e-10){
    Rad = pProj->pos.x*pProj->Mom.x + pProj->pos.y*pProj->Mom.y + pProj->pos.z*pProj->Mom.z;
    if(Rad<0) return false;
    return true;
  }
  else
    if(Rad < 1.07*pow(GetCurrA(),1./3.)+1) return false;
  return true;
}

double INCModel::GetTau(PROJ* pProj,double nSlices)
{
  double sigmap,sigman,lambda,res,fTmp;
  if(pProj->type==PROTON){
    sigmap = GetXSection(pProj->Energy,EQUAL_TYPE)*fermi*fermi/millibarn;
    sigman = GetXSection(pProj->Energy,DIFF_TYPE)*fermi*fermi/millibarn;
  }
  else{
    sigmap = GetXSection(pProj->Energy,DIFF_TYPE)*fermi*fermi/millibarn;
    sigman = GetXSection(pProj->Energy,EQUAL_TYPE)*fermi*fermi/millibarn;
  }
  fTmp = sqrt(pProj->pos.x*pProj->pos.x+pProj->pos.y*pProj->pos.y+pProj->pos.z*pProj->pos.z);
  lambda = GetCurrA()/5.1314481/(GetCurrZ()*sigmap + (GetCurrA()-GetCurrZ())*sigman);
  fTmp = sqrt(pProj->Mom.x*pProj->Mom.x+pProj->Mom.y*pProj->Mom.y+pProj->Mom.z*pProj->Mom.z);
  res = lambda/(nSlices*fTmp/((pProj->type == PROTON) ? g_dProtonMass : g_dNeutronMass));
  return res;
}

void INCModel::AddRecoilMomentum(PROJ* pProj)
{
  m_pNucleus->AddMomentum(G4ThreeVector(pProj->Mom.x,pProj->Mom.y,pProj->Mom.z));
}
void INCModel::AddRecoilMomentum(VECTOR& vMom)
{
  m_pNucleus->AddMomentum(G4ThreeVector(vMom.x,vMom.y,vMom.z));
}


void INCModel::GetRestMomentum(VECTOR& vMom)
{
  G4ThreeVector* vThreeVector; //HPW to get to compile  = m_pNucleus->GetMomentum();
  vMom.x = vThreeVector->x();
  vMom.y = vThreeVector->y();
  vMom.z = vThreeVector->z();
}

#ifndef __INC_MODEL_DEFINED__
#define __INC_MODEL_DEFINED__

class G4Nucleus;
struct PROJ;
#include "G4Nucleus.hh"
#include "globals.hh"

struct VECTOR{
  double x;
  double y;
  double z;
};

class INCModel
{
  G4Nucleus* m_pNucleus;
  unsigned   m_cParticles;
  unsigned   m_cHoles;
  unsigned   m_cCharged;
public:
  INCModel(G4Nucleus* pNuc=NULL);
  ~INCModel(){if (m_pNucleus) delete m_pNucleus;};
  void SetNucleus(G4Nucleus* pNuc);
  void AddExcitationEnergy(double En){m_pNucleus->AddExcitationEnergy(En);};
  void AddParticle(){m_cParticles++;};
  void DecParticle(){m_cParticles--;};
  unsigned GetParticles(){return m_cParticles;};
  void AddCharged(){m_cCharged++;};
  void DecCharged(){m_cCharged--;};
  unsigned GetCharged(){return m_cCharged;};
  void SetCurrZ(double Z){m_pNucleus->SetParameters(m_pNucleus->GetN(),Z);};
  double GetCurrZ(){return m_pNucleus->GetZ();};
  void SetCurrA(double A){m_pNucleus->SetParameters(A,m_pNucleus->GetZ());};
  double GetCurrA(){return m_pNucleus->GetN();};
  void AddHole(){m_cHoles++;};
  unsigned GetHoles(){return m_cHoles;};
  bool Interact(PROJ* pProj,PROJ* pNewProj,double tau); //true ako e syzdadena nova chastica za sledene
  bool CheckCrossing(PROJ* pProj,unsigned nSlices,double& tau);
  bool IsOutside(PROJ* pProj);
  double GetTau(PROJ* pProj,double nSlices);
  void AddRecoilMomentum(PROJ* pProj);
  void AddRecoilMomentum(VECTOR& vMom);
  void GetRestMomentum(VECTOR& vMom);
  double GetExcitationEnergy(){return m_pNucleus->GetEnergyDeposit();};
public:
  double GetFermiEnergy(VECTOR& v);
  double GetMinimumProjectileEnergy(unsigned nZ,unsigned nA,VECTOR& vPos);
};

#endif

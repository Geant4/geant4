
#include "MyNewClass.hh"
#include "G4NucleiProperties.hh"

void HadronElastic::CalcComb()
{
  unsigned i,j,k,l;
  m_Comb[0] = m_Comb[1] = m_Comb[2] = 1;
  l = 1;
  for(i=3;i<200;i++){
    k = (i*(i-1))>>1;
    m_Comb[k] = 1;
    for(j=1;j<i-1;j++) m_Comb[k+j] = m_Comb[l+j]+m_Comb[l+j-1];
    m_Comb[k+j] = 1;
    l = k;
  }
}

G4double HadronElastic::GetComb(unsigned A,unsigned n)
{
  if(n>A) return 0;
  unsigned k = (A*(A-1))>>1;
  return m_Comb[k+n-1];
}

G4double HadronElastic::GetMomentumHint(G4double num,G4double maxP,G4Nucleus& Nucl,int Z)
{
  G4double q = num;
  G4int A = Nucl.GetN();
  G4double rad = 1.25*fermi*pow(A,1./3.),rad2a;
  rad = 2./3.*(rad*rad-1*fermi*fermi)/(1-1./A);
  rad2a = rad + 2*0.24142*fermi*fermi;
  int i;
  G4double currRes=0,Tmp,Tmp1=42.7*millibarn*sqrt(1+0.28*0.28),Tmp2=1,Tmp3=0;
  for(i=1;i<=A;i++){
    Tmp = GetComb((unsigned)A,(unsigned)i);
    if(i%2 ==1){
      currRes += /*cos(15.64*i)*/1/i*exp(-(rad2a/(4*i)-rad/(4*A))*q)*Tmp1*Tmp;
      Tmp3 += /*sin(15.64*i)*/1/i*exp(-(rad2a/(4*i)-rad/(4*A))*q)*Tmp1*Tmp;
    }
    else{
      currRes -= /*cos(15.64*i)*/1/i*exp(-(rad2a/(4*i)-rad/(4*A))*q)*Tmp1*Tmp;
      Tmp3 -= /*sin(15.64*i)*/1/i*exp(-(rad2a/(4*i)-rad/(4*A))*q)*Tmp1*Tmp;
    }
    Tmp1 *= 42.7*millibarn*sqrt(1+0.28*0.28)/(2*pi*rad2a);
  }
  return currRes*currRes+Tmp3*Tmp3;
}

void HadronElastic::PrintComb()
{
  unsigned i,j,k;
  for(i=1;i<10;i++){
    k = (i*(i-1))>>1;
    for(j=0;j<i;j++){
      G4cout<<m_Comb[k+j]<<" ";
    }
    G4cout<<G4endl;
  }
}

G4double HadronElastic::GetMaximalValue(G4double maxP,G4Nucleus& Nucl)
{
  G4int A = Nucl.GetN();
  G4double rad = 1.25*fermi*pow(A,1./3.),rad2a;
  rad = 2./3.*(rad*rad-1*fermi*fermi)/(1-1./A);
  rad2a = rad + 2*0.24142*fermi*fermi;
  int i;
  G4double currRes1=0,Tmp,Tmp1=42.7*millibarn*sqrt(1+0.28*0.28),Tmp2=1,Tmp3=0;
  for(i=1;i<=A;i++){
    Tmp = GetComb((unsigned)A,(unsigned)i);
    if(i%2 ==1){
      currRes1 += Tmp*Tmp1/i/*cos(15.64*i)*/;
      Tmp3 += Tmp1/i*Tmp/*sin(15.64*i)*/;
    }
    else{
      currRes1 -= Tmp1/i*Tmp/*cos(15.64*i)*/;
      Tmp3 -= Tmp1/i*Tmp/*sin(15.64*i)*/;
    }
    Tmp1 *= 42.7*millibarn*sqrt(1+0.28*0.28)/(2*pi*rad2a);
  }
  return (currRes1*currRes1 + Tmp3*Tmp3);
}

G4VParticleChange* HadronElastic::ApplyYourself(const G4Track& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Initialize(aTrack);
  G4DynamicParticle* aSecondary=NULL,*aPart = aTrack.GetDynamicParticle();
  G4double maxP = aPart->GetTotalMomentum();
  maxP *= maxP;
  //Promjana na energijata
  if(targetNucleus.GetN() < 1.2){
    G4ParticleDefinition* aParticleType = aPart->GetDefinition();
    if (aParticleType == G4PionPlus::PionPlus())
      aSecondary = LightMedia.PionPlusExchange(aPart, targetNucleus);
    else if (aParticleType == G4PionMinus::PionMinus())
      aSecondary = LightMedia.PionMinusExchange(aPart, targetNucleus);
    else if (aParticleType == G4KaonPlus::KaonPlus())
      aSecondary = LightMedia.KaonPlusExchange(aPart, targetNucleus);
    else if (aParticleType == G4KaonZeroShort::KaonZeroShort())
      aSecondary = LightMedia.KaonZeroShortExchange(aPart,targetNucleus);
    else if (aParticleType == G4KaonZeroLong::KaonZeroLong())
      aSecondary = LightMedia.KaonZeroLongExchange(aPart, targetNucleus);
    else if (aParticleType == G4KaonMinus::KaonMinus())
      aSecondary = LightMedia.KaonMinusExchange(aPart, targetNucleus);
    else if (aParticleType == G4Proton::Proton())
      aSecondary = LightMedia.ProtonExchange(aPart, targetNucleus);
    else if (aParticleType == G4AntiProton::AntiProton())
      aSecondary = LightMedia.AntiProtonExchange(aPart, targetNucleus);
    else if (aParticleType == G4Neutron::Neutron())
      aSecondary = LightMedia.NeutronExchange(aPart, targetNucleus);
    else if (aParticleType == G4AntiNeutron::AntiNeutron())
      aSecondary = LightMedia.AntiNeutronExchange(aPart, targetNucleus);
    else if (aParticleType == G4Lambda::Lambda())
      aSecondary = LightMedia.LambdaExchange(aPart, targetNucleus);
    else if (aParticleType == G4AntiLambda::AntiLambda())
      aSecondary = LightMedia.AntiLambdaExchange(aPart, targetNucleus);
    else if (aParticleType == G4SigmaPlus::SigmaPlus())
      aSecondary = LightMedia.SigmaPlusExchange(aPart, targetNucleus);
    else if (aParticleType == G4SigmaMinus::SigmaMinus())
      aSecondary = LightMedia.SigmaMinusExchange(aPart, targetNucleus);
    else if (aParticleType == G4AntiSigmaPlus::AntiSigmaPlus())
      aSecondary = LightMedia.AntiSigmaPlusExchange(aPart,targetNucleus);
    else if (aParticleType == G4AntiSigmaMinus::AntiSigmaMinus())
      aSecondary= LightMedia.AntiSigmaMinusExchange(aPart,targetNucleus);
    else if (aParticleType == G4XiZero::XiZero())
      aSecondary = LightMedia.XiZeroExchange(aPart, targetNucleus);
    else if (aParticleType == G4XiMinus::XiMinus())
      aSecondary = LightMedia.XiMinusExchange(aPart, targetNucleus);
    else if (aParticleType == G4AntiXiZero::AntiXiZero())
      aSecondary = LightMedia.AntiXiZeroExchange(aPart, targetNucleus);
    else if (aParticleType == G4AntiXiMinus::AntiXiMinus())
      aSecondary = LightMedia.AntiXiMinusExchange(aPart, targetNucleus);
    else if (aParticleType == G4OmegaMinus::OmegaMinus())
      aSecondary = LightMedia.OmegaMinusExchange(aPart, targetNucleus);
    else if (aParticleType == G4AntiOmegaMinus::AntiOmegaMinus())
      aSecondary= LightMedia.AntiOmegaMinusExchange(aPart,targetNucleus);
    else if (aParticleType == G4KaonPlus::KaonPlus())
      aSecondary = LightMedia.KaonPlusExchange(aPart, targetNucleus);
  }
  if (aSecondary) {
    aSecondary->SetMomentum(aPart->GetMomentum());
    theParticleChange.SetStatusChange(fStopAndKill);
    theParticleChange.AddSecondary(aSecondary);
  }

  //  G4double rand1,rand2,tranverse,Rad,tranverse1,qbeta2,exp1,exp2;
  int i;
  static double beta = 0.63*fermi;
  double Rad = 1.5*fermi*pow(targetNucleus.GetN(),1./3.);
  double rand1,rand2,qbeta2,tranverse,tranverse1;
  // - ot statijata za eikonalnoto priblijenie
  do{
    rand1 = -log(1-G4UniformRand()*(1-exp(-2*pi*beta*maxP)))/(2*pi*beta);
    rand2 = G4UniformRand()*exp(-2*pi*beta*rand1)/(40*40*millibarn*millibarn+0.126*0.126);
    qbeta2 = rand1*rand1*beta*beta;
    tranverse = exp(-Rad/beta)/((qbeta2+1)*(qbeta2+1));
    i=1;
    do{
      i++;
      if(i%2==0){
	tranverse1 = -exp(-i*Rad/beta)*i/((qbeta2+i*i)*(qbeta2+i*i));
      }
      else{
	tranverse1 = exp(-i*Rad/beta)*i/((qbeta2+i*i)*(qbeta2+i*i));
      }
      tranverse += tranverse1;
    }
    while(fabs(tranverse1) > 1e-15);

    /* -- neshto sym zamisljal ot eikonalnoto priblijenie - izglejda za visoki energii
      i=0;
      exp1=exp2=0;
      do{
      tranverse1 = exp(-2*i*beta*rand1);
      exp1 += tranverse1;
      i++;
      }
      while(tranverse1 > 1e-10);
      i=1;
      do{
      tranverse1 = exp(-2*i*beta*rand1)*i;
      exp2 += tranverse1;
      i++;
      }
      while(tranverse1 > 1e-10);
      G4cout<<"tranv: "<<tranverse1<<" exp1: "<<exp1<<" exp2: "<<exp2<<G4endl;
    */
    /*  
	tranverse *= 2*beta*beta*beta*rand1;
	tranverse += (sin(rand1*Rad)*Rad - pi*beta*cos(rand1*Rad)*cosh(pi*beta*rand1)/sinh(pi*beta*rand1))/sinh(pi*beta*rand1);
	tranverse *= tranverse;
	tranverse /= (40*40*millibarn*millibarn + 0.126*0.126);
    */
    //Rad*sin(rand1*Rad)/sinh(pi*beta*rand1)-pi*beta*cosh(pi*beta*rand1)*cos(rand1*Rad)/(sinh(pi*beta*rand1)*sinh(pi*beta*rand1));

    // -- Glauber not working formula - rand1 e vinagi nula
    G4double Maximum=GetMaximalValue(maxP,targetNucleus);
    G4double rand1,rand2,tranverse,help;
    help = 0.24142*fermi*fermi/2/targetNucleus.GetN();
    //    G4cout<<"help = "<<help<<G4endl;
    do{
      rand1 = -log(1-G4UniformRand()*(1-exp(help*maxP)))/help;
      rand2 = G4UniformRand()*exp(-help*rand1)*Maximum;
      tranverse = GetMomentumHint(rand1,maxP,targetNucleus,
				  aPart->GetDefinition()->GetPDGCharge());
      if(tranverse - exp(-help*rand1)*Maximum > 1e-5){
	G4cout<<"Error: "<<tranverse<<" > "<<exp(-help*rand1)*Maximum<<G4endl;
      }
      if(rand1!=0) G4cout<<"ypaaaaa "<<sqrt(4*maxP)<<" "<<sqrt(rand1)<<" "<<sqrt(rand1/4/maxP)<<G4endl;
    }
    while(tranverse < rand2);
    tranverse = rand1;

    /*
      if(tranverse - exp(-2*pi*beta*rand1)/(2*pi*beta)*1e+6 > 1e-15){
      G4cout<<"tranverse: "<<tranverse<<" is greater than: "<<exp(-2*pi*beta*rand1)/(2*pi*beta)*1e+6<<G4endl;*/
      //      assert(tranverse <= exp(-2*pi*beta*rand1)*Common*(pi*beta+Rad)*(pi*beta+Rad))
    //}
  }
  while(tranverse < rand2/(2*pi*beta)*1e+6);
  tranverse = rand1;
  G4double gamma = aPart->GetMass()/(G4NucleiProperties::GetNuclearMass(targetNucleus.GetN(),targetNucleus.GetZ()));
  maxP = sqrt(maxP);
  G4double theta,phi;
  theta = 2*tranverse/(maxP*maxP);
  if(theta + gamma < 1e-5) theta = 1;
  else theta = cos(atan(sin(acos(theta)/(theta + gamma))));
//theta = 2*cos(asin(sqrt(tranverse)/2/maxP));
  //  theta = 1-0.5*tranverse/(maxP*maxP);
  //  theta = (gamma+theta)/sqrt(1+gamma*gamma+2*gamma*theta);
  if(fabs(theta) > 1) theta = 1;
  G4cout<<"Scattering at: "<<acos(theta)*180/pi<<" deg"<<G4endl;
  phi = twopi*G4UniformRand();
  if(theta < 1e-10){
    return &theParticleChange;
  }
  G4double momx,momy,momz,pxinc,pyinc,pzinc,pxnew,pynew,pznew;
  
  pxinc = maxP*(aPart->GetMomentumDirection().x());
  pyinc = maxP*(aPart->GetMomentumDirection().y());
  pzinc = maxP*(aPart->GetMomentumDirection().z());

  momx = maxP*sqrt(1-theta*theta)*sin(phi);
  momy = maxP*sqrt(1-theta*theta)*cos(phi);
  momz = maxP*theta;

  G4double pt2 = pxinc*pxinc + pyinc*pyinc;
  if (pt2 > 0.) {
    G4double cost = pzinc/maxP;
    G4double sint1 = sqrt(abs((1. - cost )*(1.+cost)));
    G4double sint2 = sqrt(pt2)/maxP;
    G4double sint = 0.5*(sint1 + sint2);
    G4double ph = pi*0.5;
    if (pyinc < 0.) ph = pi*1.5;
    if (abs(pxinc) > 1.e-6) ph = atan2(pyinc, pxinc);
    G4double cosp = cos(ph);
    G4double sinp = sin(ph);
    if (verboseLevel > 1) {
      G4cout << "cost sint " << cost << " " << sint << G4endl;
      G4cout << "cosp sinp " << cosp << " " << sinp << G4endl;
    }
    pxnew = cost*cosp*momx - sinp*momy + sint*cosp*momz;
    pynew = cost*sinp*momx + cosp*momy + sint*sinp*momz;
    pznew =     -sint*momx                 +cost*momz;
  }
  else {
    pxnew = momx;
    pynew = momy;
    pznew = momz;
  }
  momx /= maxP;
  momy /= maxP;
  momz /= maxP;
  if(aSecondary){
    aSecondary->SetMomentumDirection(G4ThreeVector(momx,momy,momz));
    //    aSecondary->SetKineticEnergy(Lost);
  }
  else{
    theParticleChange.SetMomentumChange(momx,momy,momz);
    //    theParticleChange.SetEnergyChange(Lost);
  }
  return &theParticleChange;
}

#define RUN
//#define DEBUG

#include "EquilibriumEvaporator.h"
#include "InuclNuclei.h"
#include "LorentzConvertor.h"

CollisionOutput EquilibriumEvaporator::collide(InuclParticle* bullet,
                     InuclParticle* target) {

// simple implementation of the equilibium evaporation a la Dostrowski

const double huge = 50.;
const double small = -50.;
const double one_third = 1./3.;
const double two_thirds = 2./3.;
const double prob_cut_off = 1.e-15;
const double Q1[6] = { 0., 0., 2.23, 8.49, 7.72, 28.3 };
const double AN[6] = { 1., 1., 2., 3., 3., 4. };
const double Q[6] = { 0., 1., 1., 1., 2., 2. };
const double G[6] = { 2., 2., 6., 6., 6., 4. };
const double BE = 0.0063;
const double fisssion_cut = 1000.;
#ifdef RUN
const double cut_off_energy = 0.1;
#else
const double cut_off_energy = 100000.;
#endif
const double BF = 0.0242;
const int itry_max = 1000;
const int itry_global_max = 1000;
const double small_ekin = 1.e-6;
const int itry_gam_max = 100;

vector<double> W(8);
vector<double> A1(6);
vector<double> Z1(6);
vector<double> u(6);
vector<double> V(6);
vector<double> TM(6);

double coul_coeff;

CollisionOutput output;

if(InuclNuclei* nuclei_target = dynamic_cast<InuclNuclei*>(target)) {

  double A = nuclei_target->getA();
  double Z = nuclei_target->getZ();
  vector<double> PEX = nuclei_target->getMomentum();
  double EEXS = nuclei_target->getExitationEnergy();
  if(EEXS < 0.) cout << " after noeq: eexs " << EEXS << endl;
  InuclElementaryParticle dummy(small_ekin,1);
  LorentzConvertor toTheNucleiSystemRestFrame;
  toTheNucleiSystemRestFrame.setBullet(dummy.getMomentum(),dummy.getMass());
  vector<double> ppout(4,0.);
  
  if(timeToBigBang(A,Z,EEXS)) {
    cout << " big bang in eql start " << endl;
    return theBigBanger->collide(0,target);
  }
   else {  
   
    if(A >= 100.) {
      coul_coeff = 1.4;
    }
     else {
      coul_coeff = 1.2;
    };   
   
    InuclNuclei dummy_nuc;
    double EEXS_new;
    vector<double> pin = PEX;
    pin[0] += 0.001*EEXS;
    bool try_again = true;  
    bool fission_open = true;
    double nuc_mass;  
    int itry_global = 0;
    while(try_again && itry_global < itry_global_max) {
      itry_global++;
#ifdef DEBUG
      cout << " A " << A << " Z " << Z << " EEXS " << EEXS << endl;
#endif
      nuc_mass = dummy_nuc.getNucleiMass(A,Z); 
      PEX[0] = sqrt(PEX[1]*PEX[1] + PEX[2]*PEX[2] + PEX[3]*PEX[3] +
                                             nuc_mass*nuc_mass);  	
      toTheNucleiSystemRestFrame.setTarget(PEX,nuc_mass);
      toTheNucleiSystemRestFrame.toTheTargetRestFrame();

      if(timeToBigBang(A,Z,EEXS)) { // big bang
      
#ifdef DEBUG
        cout << " big bang in eql step " << endl;
#endif
        InuclNuclei nuclei(PEX,A,Z);        
        nuclei.setExitationEnergy(EEXS);      
        CollisionOutput explosion = theBigBanger->collide(0,&nuclei);
	output.addOutgoingParticles(explosion.getOutgoingParticles());
        return output;
	
      }
       else { // normal chain
        
	if(EEXS > cut_off_energy) { 

          double E0 = getE0(A);
          double parlev = getPARLEVDEN(A,Z);
          double u1 = parlev*A;

          pair<vector<double>,vector<double> > parms = paraMaker(Z);
	  vector<double> AK = parms.first;
	  vector<double> CPA = parms.second;
        
	  double DM0 = bindingEnergy(A,Z);   

	  for(int i = 0; i < 6; i++) {
	    A1[i] = A - AN[i];
	    Z1[i] = Z - Q[i];
	    u[i] = parlev*A1[i];

            TM[i] = -0.1;
	    if(goodRemnant(A1[i],Z1[i])) {
	      double QB = DM0 - bindingEnergy(A1[i],Z1[i]) - Q1[i];
	      V[i] = coul_coeff* Z*Q[i] * AK[i] / (1. + EEXS/E0) /
	         (pow(A1[i],one_third) + pow(AN[i],one_third));
	      TM[i] = EEXS - QB - V[i]*A/A1[i];  
	    };
	  }; 
      
          double ue = 2.*sqrt(u1*EEXS);
	  double prob_sum = 0.;
	
	  if(TM[0] > cut_off_energy) {
	    double AL = getAL(A);
	    W[0] = BE*pow(A1[0],two_thirds)*G[0]*AL;
	    double TM1 = 2.*sqrt(u[0]*TM[0]) - ue;
	    if(TM1 > huge) {
	      TM1 = huge;
	    }
	     else if(TM1 < small) {
	      TM1 = small;
	    };
	    W[0] = W[0]*exp(TM1); 	     
	    prob_sum += W[0];
	  }
	   else {
	    W[0] = 0.;
	  }; 
      
          for(int i = 1; i < 6; i++) {
	    if(TM[i] > cut_off_energy) {
	      W[i] = BE*pow(A1[i],two_thirds)*G[i]*(1. + CPA[i]);
	      double TM1 = 2.*sqrt(u[i]*TM[i]) - ue;
	      if(TM1 > huge) {
	        TM1 = huge;
	      }
	       else if(TM1 < small) {
	        TM1 = small;
	      };
	      W[i] = W[i]*exp(TM1); 	     
	      prob_sum += W[i];
	    }
	     else {
	      W[i] = 0.;
	    }; 
	  };

//         fisson part
	  W[6] = 0.;
          if(A >= 100. && fission_open) {
	    double X2 = Z*Z/A;
	    double X1 = 1. - 2.*Z/A; 
	    double X = 0.019316 * X2 / (1. - 1.79*X1*X1);
	    double EF = EEXS - getQF(X,X2,A,Z,EEXS);
	  
            if(EF > 0.) {
	      double AF = u1*getAF(X,A,Z,EEXS);
	      double TM1 = 2.*sqrt(AF*EF) - ue;
	      if(TM1 > huge) {
	        TM1 = huge;
	      }
	       else if(TM1 < small) {
	        TM1 = small;
	      };
	      W[6] = BF*exp(TM1);
	      if(W[6] > fisssion_cut*W[0]) W[6] = fisssion_cut*W[0]; 	     
	      prob_sum += W[6];
	    };
	  };  
//   again time to decide what next
#ifdef DEBUG
          cout << " wn " << W[0] << " wp " << W[1] << " wd " << W[2] << endl
	       << " wh3 " << W[3] << " wt " << W[4] << " whe4 " << W[5] << endl
	       << " wfi " << W[6] << endl;
#endif
	  int icase = -1;
	  if(prob_sum < prob_cut_off) { // photon emission chain
            double UCR0 = 2.5 + 150./A;
            double T00=1./(sqrt(u1/UCR0)-1.25/UCR0);
	    int itry_gam = 0;
	    while(EEXS > cut_off_energy && try_again) {
	      itry_gam++;
	      int itry = 0;
	      double T04 = 4.*T00;
	      double FMAX;
	      if(T04 < EEXS) {
                FMAX = pow(T04,4)*exp((EEXS-T04)/T00);	        
	      }
	       else {
	        FMAX = pow(EEXS,4);
	      }; 
	      double S;
	      while(itry < itry_max) {
	        itry++;
		S = EEXS * inuclRndm();
		double X1 = pow(S,4)*exp((EEXS - S)/T00);
		if(X1 > FMAX*inuclRndm()) break;
	      };
	      if(itry < itry_max) {
//      new photon escape
		InuclElementaryParticle particle(10);
                double pmod = 0.001*S;
		vector<double> mom(4);
		pair<double,double> COS_SIN = randomCOS_SIN();
		double FI = randomPHI();
		double P1 = pmod*COS_SIN.second;
		mom[1] = P1*cos(FI);
		mom[2] = P1*sin(FI);
		mom[3] = pmod*COS_SIN.first;
		mom[0] = pmod;
		vector<double> mom_at_rest(4);
		for(int i = 1; i < 4; i++) mom_at_rest[i] = -mom[i];
                mom_at_rest[0] = sqrt(mom_at_rest[1]*mom_at_rest[1] +
		       mom_at_rest[2]*mom_at_rest[2] + 
		       mom_at_rest[3]*mom_at_rest[3] +
		     nuc_mass*nuc_mass); 
		vector<double> part_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom);
		part_mom[0] = sqrt(part_mom[1]*part_mom[1] +
		       part_mom[2]*part_mom[2] + part_mom[3]*part_mom[3]); 
		vector<double> ex_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom_at_rest);
		ex_mom[0] = sqrt(ex_mom[1]*ex_mom[1] + ex_mom[2]*ex_mom[2]
		       + ex_mom[3]*ex_mom[3] + nuc_mass*nuc_mass);			
		EEXS_new = 1000.*(PEX[0] + 0.001*EEXS - 
		                   part_mom[0] - ex_mom[0]);
		if(EEXS_new > 0.) { // everything ok
		  PEX = ex_mom;
		  EEXS = EEXS_new;
		  particle.setMomentum(part_mom);
		  output.addOutgoingParticle(particle);
	          for(int i = 0; i < 4; i++) ppout[i] += part_mom[i];
	  	}
		 else {
		  if(itry_gam == itry_gam_max) try_again = false;
		};
              }
	       else {
	        try_again = false;
	      }; 
	    };
	    try_again = false;
  	  }
	   else {
            double SL = prob_sum * inuclRndm();
	    double S1 = 0.;
            for(int i = 0; i < 7; i++) {
	      S1 += W[i]; 	
	      if(SL <= S1) {
	        icase = i;
	        break;
	      };
	    };
	    if(icase < 6) { // particle or light nuclei escape

              double uc = 2.*sqrt(u[icase]*TM[icase]);
	      double ur = (uc > huge ? exp(huge) : exp(uc));
	      double d1 = 1./ur;
	      double d2 = 1./(ur - 1.);
	    
              int itry1 = 0;
	      bool bad = true;
	      while(itry1 < itry_max && bad) {
	        itry1++; 
	      
	        int itry = 0;
	        double EPR = -1.;
	        double S;
	        while(itry < itry_max && EPR < 0.) {
	          itry++;
	          double uu = uc + log((1. - d1)*inuclRndm() + d2);
	          S = 0.5*(uc*uc - uu*uu)/u[icase];
	          EPR = TM[icase] - S*A/(A - 1.) + V[icase];
	        }; 
	    
	        if(EPR > 0. && S > 0.) { // real escape
		  S = 0.001*S;
                  if(icase < 2) { // particle escape
	            int ptype = 2 - icase;
		    InuclElementaryParticle particle(ptype);
		    double mass = particle.getMass();
//                       generate particle momentum
                    double pmod = sqrt((2.*mass + S)*S);
		    vector<double> mom(4);
		    pair<double,double> COS_SIN = randomCOS_SIN();
		    double FI = randomPHI();
		    double P1 = pmod*COS_SIN.second;
		    mom[1] = P1*cos(FI);
		    mom[2] = P1*sin(FI);
		    mom[3] = pmod*COS_SIN.first;
		    vector<double> mom_at_rest(4);
		    for(int i = 1; i < 4; i++) mom_at_rest[i] = -mom[i];
		    double new_nuc_mass = dummy_nuc.getNucleiMass(A1[icase],
		               Z1[icase]);
                    mom_at_rest[0] = sqrt(mom_at_rest[1]*mom_at_rest[1] +
		       mom_at_rest[2]*mom_at_rest[2] + 
		       mom_at_rest[3]*mom_at_rest[3] +
		                   new_nuc_mass*new_nuc_mass); 
		    mom[0] = sqrt(mom[1]*mom[1] + mom[2]*mom[2] +
		                   mom[3]*mom[3] + mass*mass);
		    vector<double> part_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom);
		    part_mom[0] = sqrt(part_mom[1]*part_mom[1] +
		       part_mom[2]*part_mom[2] + part_mom[3]*part_mom[3] +
		          	mass*mass);
		    vector<double> ex_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom_at_rest);
		    ex_mom[0] = sqrt(ex_mom[1]*ex_mom[1] + ex_mom[2]*ex_mom[2]
		          + ex_mom[3]*ex_mom[3] + new_nuc_mass*new_nuc_mass);			
		    EEXS_new = 1000.*(PEX[0] + 0.001*EEXS - 
		                         part_mom[0] - ex_mom[0]);
		    if(EEXS_new > 0.) { // everything ok
		      PEX = ex_mom;
		      EEXS = EEXS_new;
	              A = A1[icase];
	              Z = Z1[icase]; 	      
		      particle.setMomentum(part_mom);
		      output.addOutgoingParticle(particle);
	              for(int i = 0; i < 4; i++) ppout[i] += part_mom[i];
	  	      bad = false;
		    };
	          }
	           else {
		    InuclNuclei nuclei(AN[icase],Q[icase]);
		    double mass = nuclei.getMass();
//                       generate particle momentum
                    double pmod = sqrt((2.*mass + S)*S);
		    vector<double> mom(4);
		    pair<double,double> COS_SIN = randomCOS_SIN();
		    double FI = randomPHI();
		    double P1 = pmod*COS_SIN.second;
		    mom[1] = P1*cos(FI);
		    mom[2] = P1*sin(FI);
		    mom[3] = pmod*COS_SIN.first;
		    vector<double> mom_at_rest(4);
		    for(int i = 1; i < 4; i++) mom_at_rest[i] = -mom[i];
		    double new_nuc_mass = dummy_nuc.getNucleiMass(A1[icase],
		                Z1[icase]);
                    mom_at_rest[0] = sqrt(mom_at_rest[1]*mom_at_rest[1] +
		                mom_at_rest[2]*mom_at_rest[2] + 
		                mom_at_rest[3]*mom_at_rest[3] +
		         new_nuc_mass*new_nuc_mass); 
		    mom[0] = sqrt(mom[1]*mom[1] + mom[2]*mom[2] +
		           mom[3]*mom[3] + mass*mass);
		    vector<double> part_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom);
		    part_mom[0] = sqrt(part_mom[1]*part_mom[1] +
		       part_mom[2]*part_mom[2] + part_mom[3]*part_mom[3] +
		       	mass*mass);
		    vector<double> ex_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom_at_rest);
		    ex_mom[0] = sqrt(ex_mom[1]*ex_mom[1] + ex_mom[2]*ex_mom[2]
		       + ex_mom[3]*ex_mom[3] + new_nuc_mass*new_nuc_mass);			
		    EEXS_new = 1000.*(PEX[0] + 0.001*EEXS - 
		                   part_mom[0] - ex_mom[0]);
		    if(EEXS_new > 0.) { // everything ok
		      PEX = ex_mom;
		      EEXS = EEXS_new;
	              A = A1[icase];
	              Z = Z1[icase]; 	      
	              for(int i = 0; i < 4; i++) ppout[i] += part_mom[i];
		      nuclei.setMomentum(part_mom);
                      nuclei.setExitationEnergy(0.);
		      nuclei.setEnergy();
		      output.addTargetFragment(nuclei);
	  	      bad = false;
	  	    };
	          };
		};
	      };
              if(itry1 == itry_max || bad)  try_again = false;
	    }
	     else { // fission

              InuclNuclei nuclei(A,Z);        
              nuclei.setExitationEnergy(EEXS);
#ifdef DEBUG
	      cout << " fission: A " << A << " Z " << Z << " eexs " << EEXS <<
	       " Wn " << W[0] << " Wf " << W[6] << endl;
#endif
	      CollisionOutput foutput = theFissioner->collide(0,&nuclei);
	      vector<InuclNuclei> nuclea = foutput.getNucleiFragments();
	      if(nuclea.size() == 2) { // fission o'k
//                convert back to the lab
                for(int i = 0; i < 2; i++) {
		  vector<double> mom = nuclea[i].getMomentum();
		  mom = toTheNucleiSystemRestFrame.backToTheLab(mom);
		  nuclea[i].setMomentum(mom);
		  nuclea[i].setEnergy();
		};	      
	        CollisionOutput output1 = this->collide(0,&nuclea[0]);
		output.addOutgoingParticles(output1.getOutgoingParticles());
		output.addTargetFragments(output1.getNucleiFragments());
	        output1 = this->collide(0,&nuclea[1]);
		output.addOutgoingParticles(output1.getOutgoingParticles());
		output.addTargetFragments(output1.getNucleiFragments());
                return output;
	      }
	       else { // fission forbidden now
	        fission_open = false;
	      }; 
	    }; 
	  }; 	     
        }
         else {
          try_again = false;
        }; 
      }; 
    };
//   this time it's final nuclei
    if(itry_global == itry_global_max) cout << " ! itry_global " <<
        itry_global_max << endl;
    vector<double> pnuc(4);
    for(int i = 1; i < 4; i++) pnuc[i] = pin[i] - ppout[i];
    InuclNuclei nuclei(pnuc,A,Z);
    nuclei.setEnergy();
    pnuc = nuclei.getMomentum(); 
    double eout = pnuc[0] + ppout[0];  
    double eex_real = 1000.*(pin[0] - eout);        
    nuclei.setExitationEnergy(eex_real);
    output.addTargetFragment(nuclei);
  };
} 
 else {
    cout << " EquilibriumEvaporator -> target is not nuclei " << endl;    
}; 

  return output;
}		     

bool EquilibriumEvaporator::timeToBigBang(double a, double z, double e) const {

const double be_cut = 3.;

  bool bigb = true;
  if(a >= 12 && z >= 6. && z < 3.*(a-z)) {
    bigb = false;  
  }
   else {
    if(e < be_cut*bindingEnergy(a,z)) bigb = false;
  };
  return bigb;
}

bool EquilibriumEvaporator::goodRemnant(double a, double z) const {
  return a > 1. && z > 0. && a > z;
}

double EquilibriumEvaporator::getQF(double x, double x2, double a,
        double z, double e) const {

const double QFREP[72] = {  
//     TL201 *     *   *    *
//      1    2     3   4    5
     22.5, 22.0,21.0,21.0,20.,
//     BI209 BI207 PO210 AT213 *    TH234
//      6     7    8     9     10   11
     20.6, 20.6, 18.6, 15.8, 13.5,6.5,
//     TH233 TH232 TH231 TH230 TX229 PA233 PA232 PA231 PA230 U240
//     12    13    14    15    16    17    18    19    20    21
     6.65, 6.22, 6.27, 6.5,  6.7,  6.2,  6.25, 5.9,  6.1,  5.75,
//     U239 U238 U237  U236 U235 U234 U233 U232 U231
//     22   23   24    25   26   27   28   29   30
     6.46,5.7, 6.28, 5.8, 6.15,5.6, 5.8, 5.2, 5.8,
//     NP238 NP237 NP236 NP235 PU245 NP234  PU244 NP233
//     31    32    33    34    35    36     37    38
     6.2 , 5.9 , 5.9,  6.0,  5.8,  5.7,   5.4,  5.4,
//     PU242 PU241 PU240 PU239 PU238 AM247 PU236 AM245 AM244 AM243
//     39    40    41    42    43    44    45    46    47    48
     5.6,  6.1,  5.57, 6.3 , 5.5,  5.8,  4.7,  6.2,  6.4,  6.2,
//     AM242 AM241 AM240 CM250 AM239 CM249 CM248 CM247 CM246
//     49    50    51    52    53    54    55    56    57
     6.5,  6.2,  6.5,  5.3,  6.4,  5.7,  5.7,  6.2,  5.7,
//     CM245 CM244 CM243 CM242 CM241 BK250 CM240
//     58    59    60    61    62    63    64
     6.3,  5.8,  6.7,  5.8,  6.6,  6.1,  4.3,
//     BK249 CF252 CF250 CF248 CF246 ES254 ES253 FM254
//     65    66    67    68    69    70    71    72
     6.2,  3.8,  5.6,  4.0,  4.0,  4.2,  4.2,  3.5 };
     
const double XREP[72] = {
//      1      2     3      4      5
     0.6761,0.677,0.6788,0.6803,0.685,
//      6     7     8     9     10     11
     0.6889,.6914,.6991,.7068,0.725,.7391,
//     12  13    14   15   16    17  18    19    20    21
     0.74,.741,.742,.743,.744,.7509,.752,.7531,.7543,.7548,
//     22    23    24
     0.7557,.7566,.7576,
//      25     26   27    28    29   30   31    32     33    34
     0.7587,.7597,.7608,.762,.7632,.7644,.7675,.7686,.7697,.7709,
//      35    36    37    38    39   40    41
     0.7714,.7721,.7723,.7733,.7743,.7753,.7764,
//      42    43    44    45    46    47    48   49
     0.7775,.7786,.7801,.781,.7821,.7831,.7842,.7852,
//     50     51    52    53    54    55    56    57    58
     0.7864,.7875,.7880,.7887,.7889,.7899,.7909,.7919,.7930,
//      59    60    61    62    63    64
     0.7941,.7953,.7965,.7977,.7987,.7989,
//      65    66    67    68    69    70    71    72
     0.7997,.8075,.8097,.8119,.8143,.8164,.8174,.8274 };

const double G0 = 20.4;
const double XMIN = 0.6761;
const double XMAX = 0.8274;

double QFF;

if(x < XMIN || x > XMAX) {
  double X1 = 1. - 0.02*x2;
  double FX=(.73 + (3.33 * X1 - 0.66) * X1) * pow(X1,3);
  QFF = G0 * FX * pow(a,2./3.);
}
 else {
  for(int i = 1; i < 72; i++) {
    if(x <= XREP[i]) {
      QFF = QFREP[i-1] + (QFREP[i] - QFREP[i-1])*(x-XREP[i-1])/
              (XREP[i] - XREP[i-1]);
      break;
    };
  };
};
if(QFF < 0.) QFF = 0.;

return QFF; 
}

double EquilibriumEvaporator::getAF(double x, double a, double z, double e) const {
// ugly parameterisation to fit the experimental fission cs for Hg - Bi nuclei
   double AF = 1.285*(1. - e/1100.);
   if(AF < 1.06) AF = 1.06;
//   if(AF < 1.08) AF = 1.08;
   return AF;
}	

double EquilibriumEvaporator::getPARLEVDEN(double A, double Z) const {
const double par = 0.125;
return par;
}

double EquilibriumEvaporator::getE0(double A) const {
const double e0 = 200.;   
return e0;   
}

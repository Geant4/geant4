#define RUN
//#define DEBUG

#include <math.h>

#include "NonEquilibriumEvaporator.h"
#include "InuclElementaryParticle.h"
#include "InuclNuclei.h"
#include "LorentzConvertor.h"

CollisionOutput NonEquilibriumEvaporator::collide(InuclParticle* bullet,
                     InuclParticle* target) {

const double one_third = 1./3.;
const double a_cut = 5.;
const double z_cut = 3.;
#ifdef RUN
const double eexs_cut = 0.1;
#else
const double eexs_cut = 100000.;
#endif
const double coul_coeff = 1.4;
const int itry_max = 1000;
const double small_ekin = 1.e-6;
const double width_cut = 0.005;

CollisionOutput output;

if(InuclNuclei* nuclei_target = dynamic_cast<InuclNuclei*>(target)) {

//  initialization
  
  double A = nuclei_target->getA();
  double Z = nuclei_target->getZ();
  vector<double> PEX = nuclei_target->getMomentum();
  vector<double> pin = PEX;
  double EEXS = nuclei_target->getExitationEnergy();
  pin[0] += 0.001*EEXS;
  InuclNuclei dummy_nuc;
  ExitonConfiguration config = nuclei_target->getExitonConfiguration();  
  double QPP = config.protonQuasiParticles;
  double QNP = config.neutronQuasiParticles; 
  double QPH = config.protonHoles;
  double QNH = config.neutronHoles; 
  double QP = QPP + QNP;
  double QH = QPH + QNH;
  double QEX = QP + QH;
  InuclElementaryParticle dummy(small_ekin,1);
  LorentzConvertor toTheExitonSystemRestFrame;
  toTheExitonSystemRestFrame.setBullet(dummy.getMomentum(),dummy.getMass());
  double EFN = FermiEnergy(A,Z,0);
  double EFP = FermiEnergy(A,Z,1);
  double SNN1 = csNN(EFP);
  double SNN2 = csNN(EFN);
  double SPN1 = csPN(EFP);
  double SPN2 = csPN(EFN);
  double AR = A - QP;
  double ZR = Z - QPP;  
  int NEX = int(QEX + 0.5);
  vector<double> ppout(4,0.);
  
  bool try_again = NEX > 0 ? true : false;
  
  while(try_again) {

    if(A >= a_cut && Z >= z_cut && EEXS > eexs_cut) { // ok

#ifdef DEBUG
      cout << " A " << A << " Z " << Z << " EEXS " << EEXS << endl; 
#endif
//        update exiton system
      double nuc_mass = dummy_nuc.getNucleiMass(A,Z); 
      PEX[0] = sqrt(PEX[1]*PEX[1] + PEX[2]*PEX[2] + PEX[3]*PEX[3] +
                                             nuc_mass*nuc_mass);  	
      toTheExitonSystemRestFrame.setTarget(PEX,nuc_mass);
      toTheExitonSystemRestFrame.toTheTargetRestFrame();
      
      double MEL = getMatrixElement(A);
      double E0 = getE0(A);
      double PL = getParLev(A,Z);
      double parlev = PL/A;
      double EG = PL*EEXS;

      if(QEX < sqrt(2.*EG)) { // ok

        pair<double,double> parms = paraMakerTruncated(Z);
	double AK1 = parms.first;
	double CPA1 = parms.second;

	double VP = coul_coeff* Z* AK1 / (pow(A-1.,one_third) + 1.) /
	               (1. + EEXS/E0);
        double DM1 = bindingEnergy(A,Z);
	double BN = DM1 - bindingEnergy(A-1.,Z);
	double BP = DM1 - bindingEnergy(A-1.,Z-1.);
	
	double EMN = EEXS - BN;
	double EMP = EEXS - BP - VP*A/(A - 1.);
	double ESP;
	
	if(EMN > eexs_cut) { // ok
          
	  int icase = 0;
	  
	  if(NEX > 1) {
	    double APH = .25*(QP*QP+QH*QH+QP-3.*QH);
	    double APH1 = APH + 0.5*(QP+QH);
	    ESP = EEXS/QEX;
	    double MELE = MEL/ESP/pow(A,3);

	    if(ESP > 15.) {
	      MELE *= sqrt(15./ESP);
	    }
	     else if(ESP < 7.) {
	      MELE *= sqrt(ESP/7.);
	      if(ESP < 2.) MELE *= sqrt(ESP/2.);
	    };    
            
	    double F1 = EG - APH;
	    double F2 = EG - APH1;
	    
	    if(F1 > 0. && F2 > 0.) {
	      double F = F2/F1;
	      double M1 = 2.77*MELE*PL;

	      vector<double> D(3,0.);
	      D[0] = M1 * F2 * F2 * pow(F,NEX-1) / (QEX+1.);

	      if(D[0] > 0.) {
	        if(NEX >= 2) {
		  D[1] = 0.0462 / parlev / pow(A,one_third) * QP * EEXS / QEX;
		  if(EMP > eexs_cut) 
		    D[2] = D[1] * pow(EMP/EEXS,NEX) * (1.+CPA1);
		  D[1] *= pow(EMN/EEXS,NEX)*getAL(A);   
		  if(QNP < 1.) D[1] = 0.;
		  if(QPP < 1.) D[2] = 0.;

		  try_again = NEX > 1 && (D[1] > width_cut*D[0] || 
		           D[2] > width_cut*D[0]);
		  if(try_again) {
		    double D5 = D[0] + D[1] + D[2];
		    double SL = D5 * inuclRndm();
		    double S1 = 0.;
		    for(int i = 0; i < 3; i++) {
		      S1 += D[i]; 	
		      if(SL <= S1) {
		        icase = i;
			break;
		      };
		    };
		  }; 
		};
	      }
	       else {
                try_again = false;
	      }; 
	    }
	     else {
              try_again = false;
	    }; 
	  }; 
 	  
	  if(try_again) {
            if(icase > 0) { // N -> N - 1 with particle escape
	      double V;
	      int ptype;
	      double B;
	      if(A < 3.) try_again = false;
	      if(try_again) { 
	        if(icase == 1) { // neutron escape
	          if(QNP < 1.) { 
		    icase = 0;
		  }
		   else {
		    B = BN;
		    V = 0.;
		    ptype = 2;		  
		  };    
	        }
	         else { // proton esape
	          if(QPP < 1.) { 
		    icase = 0;
		  }
		   else {
		    B = BP;
		    V = VP;
		    ptype = 1;
		    if(Z - 1. < 1.) try_again = false;
		  };   
	        };
	        
		if(try_again && icase != 0) {
		  double EB = EEXS - B;
		  double E = EB - V*A / (A-1.);
		  if(E < 0.) {
		    icase = 0;
		  }
		   else {
		    double E1 = EB - V;
		    double EEXS_new = -1.;
		    double EPART;

		    int itry1 = 0;
		    bool bad = true;
		    while (itry1 < itry_max && icase > 0 && bad) {
		      itry1++;
		      
		      int itry = 0;		    
		      while (EEXS_new < 0. && itry < itry_max) {
		        itry++;
		        double R = inuclRndm();
		        double X;
		        if(NEX == 2) {
		          X = 1. - sqrt(R);
		        }
		         else {
		          double QEX2 = 1./QEX;
			  double QEX1 = 1./(QEX - 1.);
			  X = pow(0.5*R,QEX2);
			  for(int i = 0; i < 1000; i++) {
			    double DX = X * QEX1 * 
			      (1. + QEX2*X*(1. - R/pow(X,NEX))/(1. - X));
			    X -= DX;
			    if(fabs(DX/X) < 0.01) break;  
			  };
		        }; 
		        EPART = EB - X*E1;
		        EEXS_new = EB - EPART*A/(A - 1.);
		      };
                    
		      if(itry == itry_max || EEXS_new < 0.) {
		        icase = 0;
		      }
		       else { // real escape
		        
			InuclElementaryParticle particle(ptype);
		        double mass = particle.getMass();
		        EPART *= 0.001; // to the GeV
//                       generate particle momentum
                        double pmod = sqrt(EPART*(2.*mass + EPART));
		        vector<double> mom(4);
		        pair<double,double> COS_SIN = randomCOS_SIN();
		        double FI = randomPHI();
		        double P1 = pmod*COS_SIN.second;
		        mom[1] = P1*cos(FI);
		        mom[2] = P1*sin(FI);
		        mom[3] = pmod*COS_SIN.first;
		        vector<double> mom_at_rest(4);
		        for(int i = 1; i < 4; i++) mom_at_rest[i] = -mom[i];
                        double QPP_new = QPP;
		        double Z_new = Z;
		        
			if(ptype == 1) {
		          QPP_new -= 1.;
                          Z_new -= 1.;
		        };
		        double QNP_new = QNP;
		        if(ptype == 2) QNP_new -= 1.;
		        double A_new = A - 1.;
		     
		        double new_exiton_mass =
		                   dummy_nuc.getNucleiMass(A_new,Z_new);
	                mom_at_rest[0] = sqrt(mom_at_rest[1]*mom_at_rest[1] +
		            mom_at_rest[2]*mom_at_rest[2] + 
		            mom_at_rest[3]*mom_at_rest[3] +
		               new_exiton_mass*new_exiton_mass); 
		        mom[0] = sqrt(mom[1]*mom[1] + mom[2]*mom[2] +
		               mom[3]*mom[3] + mass*mass);
		        vector<double> part_mom = 
		            toTheExitonSystemRestFrame.backToTheLab(mom);
		        part_mom[0] = sqrt(part_mom[1]*part_mom[1] +
		          part_mom[2]*part_mom[2] + part_mom[3]*part_mom[3] +
		       	                mass*mass);
		        vector<double> ex_mom = 
		          toTheExitonSystemRestFrame.backToTheLab(mom_at_rest);
		        ex_mom[0] = sqrt(ex_mom[1]*ex_mom[1] + ex_mom[2]*ex_mom[2]
		       + ex_mom[3]*ex_mom[3] + new_exiton_mass*new_exiton_mass);			

//             check energy conservation and set new exitation energy
		        EEXS_new = 1000.*(PEX[0] + 0.001*EEXS - 
		                         part_mom[0] - ex_mom[0]);
		        if(EEXS_new > 0.) { // everything ok
		          particle.setMomentum(part_mom);
		          output.addOutgoingParticle(particle);
			  for(int i = 0; i < 4; i++) ppout[i] += part_mom[i];
		          A = A_new;
		          Z = Z_new;
			  PEX = ex_mom;
			  EEXS = EEXS_new;
		          NEX -= 1;
		          QEX -= 1;
		          QP -= 1.;
			  QPP = QPP_new;
			  QNP = QNP_new;
		          bad = false;
			};
		      };  	
		    };   
		    if(itry1 == itry_max) icase = 0;
		  };   
		};
	      }; 
	    };
	    if(icase == 0 && try_again) { // N -> N + 2 
	      double TNN = 1.6*EFN + ESP;
              double TNP = 1.6*EFP + ESP;
	      double XNUN = 1./(1.6 + ESP/EFN);
	      double XNUP = 1./(1.6 + ESP/EFP);
	      double SNN1 = csNN(TNP)*XNUP;
	      double SNN2 = csNN(TNN)*XNUN;
	      double SPN1 = csPN(TNP)*XNUP;
	      double SPN2 = csPN(TNN)*XNUN;
	      double PP = (QPP*SNN1 + QNP*SPN1)*ZR;
	      double PN = (QPP*SPN2 + QNP*SNN2)*(AR - ZR);
	      double PW = PP + PN;

	      NEX += 2;
	      QEX += 2.; 
	      QP += 1.;
	      QH += 1.;
	      AR -= 1.;
	      if(AR > 1.) {
		double SL = PW * inuclRndm();
                if(SL > PP) {
		  QNP += 1.;
		  QNH += 1.;
		}
		 else {
		  QPP += 1.;
		  QPH += 1.;
		  ZR -= 1.;
		  if(ZR < 2.) try_again = false;
		};    
	      }
	       else {
		try_again = false;
	      };
	    };
	  };
	}
	 else {
          try_again = false;
	};
      }
       else {
        try_again = false;
      }; 
    }
     else {
      try_again = false;
    };  
  };

//           everything finished, set output nuclei
//   the exitation energy has to be re-set properly for the energy
//      conservation
  vector<double> pnuc(4);
  for(int i = 1; i < 4; i++) pnuc[i] = pin[i] - ppout[i];
  InuclNuclei nuclei(pnuc,A,Z);
  nuclei.setEnergy();
  pnuc = nuclei.getMomentum(); 
  double eout = pnuc[0] + ppout[0];  
  double eex_real = 1000.*(pin[0] - eout);        
  nuclei.setExitationEnergy(eex_real);
  output.addTargetFragment(nuclei);
} 
 else {
  cout << " NonEquilibriumEvaporator -> target is not nuclei " << endl;    
}; 

return output;
}

double NonEquilibriumEvaporator::getMatrixElement(double A) const {
  double me;
  if(A > 150.) {
    me = 100.;
  }
   else if(A > 20.) {
    me = 140.;
  }
   else {
    me = 70.;
  };  
  return me;
}

double NonEquilibriumEvaporator::getE0(double A) const {
const double e0 = 200.;
return e0;   
}

double NonEquilibriumEvaporator::getParLev(double A, double Z) const {
const double par = 0.125;
double pl = 0.125 * A;
return pl; 

}


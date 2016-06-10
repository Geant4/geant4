#include "CLHEP/Random/RandGaussZiggurat.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include <iostream>
#include <cmath>	// for std::log()

namespace CLHEP {

CLHEP_THREAD_LOCAL unsigned long RandGaussZiggurat::kn[128], RandGaussZiggurat::ke[256];
CLHEP_THREAD_LOCAL float RandGaussZiggurat::wn[128],RandGaussZiggurat::fn[128],RandGaussZiggurat::we[256],RandGaussZiggurat::fe[256];
CLHEP_THREAD_LOCAL bool RandGaussZiggurat::ziggurat_is_init = false;

HepRandomEngine & RandGaussZiggurat::engine() {return RandGauss::engine();}

RandGaussZiggurat::~RandGaussZiggurat() {
}

std::string RandGaussZiggurat::name() const 
{
  return "RandGaussZiggurat";
}

bool RandGaussZiggurat::ziggurat_init()
{  
  const double rzm1 = 2147483648.0, rzm2 = 4294967296.;
  double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
  double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
  int i;

/* Set up tables for RNOR */
  q=vn/std::exp(-.5*dn*dn);
  kn[0]=(unsigned long)((dn/q)*rzm1);
  kn[1]=0;

  wn[0]=q/rzm1;
  wn[127]=dn/rzm1;

  fn[0]=1.;
  fn[127]=std::exp(-.5*dn*dn);

  for(i=126;i>=1;i--) {
    dn=std::sqrt(-2.*std::log(vn/dn+std::exp(-.5*dn*dn)));
    kn[i+1]=(unsigned long)((dn/tn)*rzm1);
    tn=dn;
    fn[i]=std::exp(-.5*dn*dn);
    wn[i]=dn/rzm1;
  }

/* Set up tables for REXP */
  q = ve/std::exp(-de);
  ke[0]=(unsigned long)((de/q)*rzm2);
  ke[1]=0;

  we[0]=q/rzm2;
  we[255]=de/rzm2;

  fe[0]=1.;
  fe[255]=std::exp(-de);

  for(i=254;i>=1;i--) {
    de=-std::log(ve/de+std::exp(-de));
    ke[i+1]= (unsigned long)((de/te)*rzm2);
    te=de;
    fe[i]=std::exp(-de);
    we[i]=de/rzm2;
  }
  ziggurat_is_init=true;
  
  //std::cout<<"Done RandGaussZiggurat::ziggurat_init()"<<std::endl;
  
  return true;
}

float RandGaussZiggurat::ziggurat_nfix(long hz,HepRandomEngine* anEngine)
{
  if(!ziggurat_is_init) ziggurat_init();
  const float r = 3.442620f;     /* The start of the right tail */
  float x, y;
  unsigned long iz=hz&127;
  for(;;)
  {  
    x=hz*wn[iz];      /* iz==0, handles the base strip */
    if(iz==0) {
      do { 
        /* change to (1.0 - UNI) as argument to std::log(), because CLHEP generates [0,1), 
           while the original UNI generates (0,1] */
        x=-std::log(1.0 - ziggurat_UNI(anEngine))*0.2904764; /* .2904764 is 1/r */
        y=-std::log(1.0 - ziggurat_UNI(anEngine));
      }	while(y+y<x*x);
      return (hz>0)? r+x : -r-x;
    }
    /* iz>0, handle the wedges of other strips */
    if( fn[iz]+(1.0 - ziggurat_UNI(anEngine))*(fn[iz-1]-fn[iz]) < std::exp(-.5*x*x) ) return x;

    /* initiate, try to exit for(;;) for loop*/
    hz=(signed)ziggurat_SHR3(anEngine);
    iz=hz&127;
    if((unsigned long)std::abs(hz)<kn[iz]) return (hz*wn[iz]);
  }
}

double RandGaussZiggurat::operator()() {
  return ziggurat_RNOR(localEngine.get()) * defaultStdDev + defaultMean;
}

double RandGaussZiggurat::operator()( double mean, double stdDev ) {
  return ziggurat_RNOR(localEngine.get()) * stdDev + mean;
}

void RandGaussZiggurat::shootArray( const int size, float* vect, float mean, float stdDev )
{
   for (int i=0; i<size; ++i) {
     vect[i] = shoot(mean,stdDev);
   }
}

void RandGaussZiggurat::shootArray( const int size, double* vect, double mean, double stdDev )
{
   for (int i=0; i<size; ++i) {
     vect[i] = shoot(mean,stdDev);
   }
}

void RandGaussZiggurat::shootArray( HepRandomEngine* anEngine, const int size, float* vect, float mean, float stdDev )
{
   for (int i=0; i<size; ++i) {
     vect[i] = shoot(anEngine,mean,stdDev);
   }
}

void RandGaussZiggurat::shootArray( HepRandomEngine* anEngine, const int size, double* vect, double mean, double stdDev )
{
   for (int i=0; i<size; ++i) {
     vect[i] = shoot(anEngine,mean,stdDev);
   }
}

void RandGaussZiggurat::fireArray( const int size, float* vect)
{
   for (int i=0; i<size; ++i) {
     vect[i] = fire( defaultMean, defaultStdDev );
   }
}

void RandGaussZiggurat::fireArray( const int size, double* vect)
{
   for (int i=0; i<size; ++i) {
     vect[i] = fire( defaultMean, defaultStdDev );
   }
}

void RandGaussZiggurat::fireArray( const int size, float* vect, float mean, float stdDev )
{
   for (int i=0; i<size; ++i) {
     vect[i] = fire( mean, stdDev );
   }
}

void RandGaussZiggurat::fireArray( const int size, double* vect, double mean, double stdDev )
{
   for (int i=0; i<size; ++i) {
     vect[i] = fire( mean, stdDev );
   }
}

std::ostream & RandGaussZiggurat::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  os << " " << name() << "\n";
  RandGauss::put(os);
  os.precision(pr);
  return os;
}

std::istream & RandGaussZiggurat::get ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != name()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read state of a "
    	      << name() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  RandGauss::get(is);
  return is;
}

}  // namespace CLHEP

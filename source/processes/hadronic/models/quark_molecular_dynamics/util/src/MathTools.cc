#include "MathTools.hh"
#include "Error.hh"

REAL GaussIntegration::Integral(REAL x1)
{
	int j;
	REAL xr,xm,dx,s;
	static REAL x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.97390652};
	static REAL w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.06667134};

	xm=0.5*(x1+x0);
	xr=0.5*(x1-x0);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*(ToBeIntegrated(xm+dx)+ToBeIntegrated(xm-dx));
	}
	return s *= xr;
}

REAL GaussLaguerre::Integral()
{
	int j;
	REAL s;
	static REAL x[]={0.0,0.222846604179,1.188932101673,2.992736326059,
		5.775143569105,9.837467418383,15.982873980602};
	static REAL w[]={0.0,0.573535507423,1.36925259071,2.26068459338,
		3.35052458236,4.88682680021,7.84901594560};

	s=0;
	for (j=1;j<=6;j++) {
		s += w[j]*ToBeIntegrated(x[j]);
	}
	return s;
}

bool InverseFunction::reducePrecision = true;

REAL InverseFunction::Inverse(REAL y) const 
{
        REAL z;
        try { 
	  REAL xx;
	  int n = DivideInterval(x,y,xx,0,1e-2);
	  if ( n==0 )
	    z = FindRoot(x.a,x.b,y,Epsilon); 
	  else 
	    if ( n>0 ) {
	      REAL z_array[2];
	      z_array[0] = FindRoot(x.a,xx,y,Epsilon);
	      z_array[1] = FindRoot(xx,x.b,y,Epsilon);
	      throw AmbiguousResult<REAL>(2,z_array);
	    }
	  else
	    throw "Zero no found";
	    //	    throw Zero_Not_Found(x.a,x.b,ToBeInverted(x.a)-y,ToBeInverted(x.b)-y);
	}
	catch (Zero_Not_Found& m) {
	  if ( !reducePrecision ) 
	    throw;
	  m.writeMessage(cerr);
	  cerr << "Checking for accuracy 1e-4...";
	  if (fabs(m.fa)<1e-4)
	    z = m.a;
	  else if (fabs(m.fb)<1e-4)
	    z = m.b;
	  else {
	    cerr << "FAILED!\n";
	    throw;
	  }
	  cerr << "OK\n" << endl;
	}
	catch (...) {
	  throw;
	}
        return z;
}

REAL InverseFunction::FindRoot(REAL x0_,REAL x1_,REAL y,double Accuracy) const
{
  const int N = 100;

//	cerr << "beginning FINDZERO: " << endl;
//	cerr << "FINDZERO: x1=" << x0 << ", x2=" << x1 << endl;
  REAL x,f2; 
  int n=0;
  REAL f0 = ToBeInverted(x0_)-y;
  if (f0 == 0.0) 
    return x0_;
  REAL f1 = ToBeInverted(x1_)-y;
  if (f1 == 0.0)
    return x1_;
  //  cerr << x0_ << "  " << x1_ << "  " << f0 << "  " << f1 << endl;
  if ( double(f0*f1) < 0.0 ) {
    do
      {
	f2 = ToBeInverted(x = (x0_+x1_)/2)-y;
	if (f2 == 0) {
	  x1_ = x;
	  break;
	}
	if (double(f0*f2) < 0.0)
	  {
	    x1_ = x;
	    f1 = f2;
	  }
	else
	  {
	    x0_ = x;
	    f0 = f2;
	  }
	//	cout << x0_ << "," << x1_ << ":" << fabs(f2) << endl;
      }
    while ((fabs(x1_-x0_) >= Accuracy ) && ++n<N);
    if (n>=N)
      {
	//	cerr << "FINDZERO: too much iterations" << endl;
	throw TooMuchIter(x0_,x1_,f0,f1);
      }
  }
  else
    {
      //      cerr << "no zero detected: y = " << y << ", " << f0 << "," << f1 << endl;
      //      x1_ = (fabs(f0) < fabs(f1)) ? x0_ : x1_	   ;
      if (f0 == 0)
	return x0_;
      else if (f1 == 0)
	return x1_;
      else
	throw Zero_Not_Found(x0_,x1_,f0,f1);
    }
  //	cerr << "ending FINDZERO: x = " << x1_ << endl;
  return x1_;

}

int InverseFunction::DivideInterval(IntervalRange xx,REAL y,REAL& Result,int n,double eps) const
{
  REAL z = 0.5*(xx.a+xx.b);
  REAL fa = ToBeInverted(xx.a)-y;
  REAL fb = ToBeInverted(xx.b)-y;
  //  cerr << n << "  " << xx.a << "  " << xx.b << "  " << fa+y << "  "  << fb+y << endl;
  if (fabs(xx.b-xx.a)>=eps && double(fa*fb)>0.0 ) {
    REAL y2;
    IntervalRange x1(xx.a,z);
    int n1 = DivideInterval(x1,y,Result,n+1,eps);
    int n2 = -1;
    if (n1 == -1) {
      IntervalRange x2(z,xx.b);
      n2 = DivideInterval(x2,y,y2,n+1,eps);
    }
    n = -1;
    if (n1>-1 && n2>-1)
      if (n1!=n2)
	throw Error();
      else
	n = n1;
    else {
      if (n1>-1)
	n = n1;
      if (n2>-1)
	n = n2;
    }
  }
  else {
    Result = xx.b;
  }
  return n;
};



REAL InverseFunctionWithNewton::FindRoot(REAL x0_,REAL,REAL y,double Accuracy) const
{
  REAL x = x0_, x1 = x;
  
  const int N = 100;

  //	cerr << "beginning FINDZERO: " << endl;
  //	cerr << "FINDZERO: x1=" << x << endl;
  int n=0;
  double f0 = 0.0;
  double f10 = 0.0;
  do {
    x1 = x;
    f0 = ToBeInverted(x)-y;
    if (f0 == 0.0) 
      return x;
    f10 = Derivative(x);
    if (f10==0.0) 
      throw Zero_Not_Found(x,x,f0,f10);
    x -= f0/f10;
    //    cerr << x << "  " << x1 << "  " << f0 << "  " << f10<< endl;
  }
  while ( fabs(f0) > Accuracy  && fabs(x-x1)>Accuracy && ++n<N);
  //  cout << x << "  " << f0 << endl;
  if ( fabs(f0) > Accuracy )
    {
      //	cerr << "FINDZERO: too much iterations" << endl;
      throw TooMuchIter(x0_,x0_,f0,f10);
    }
  return x;
}



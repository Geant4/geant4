#ifndef MYGAMMA_H
#define MYGAMMA_H



class MyGamma {
	public:
	MyGamma ();
	~MyGamma();
	
	double Gamma(double z);
	double  Gamma(double a,double x);
	private:
	double GamCf(double a,double x);
	double GamSer(double a,double x);
	
	// Abs
	static short  Abs(short d)  { return (d > 0) ? d : -d; }
	static int    Abs(int d)   { return (d > 0) ? d : -d; }
	static long   Abs(long d)  { return (d > 0) ? d : -d; }
	static float  Abs(float d){ return (d > 0) ? d : -d; }
	static double Abs(double d)   { return (d > 0) ? d : -d; }
	static double LnGamma(double z); 
	static double Log(double x){ return log(x); }
	static double Exp(double x){ return exp(x); }
	
	
};



#endif

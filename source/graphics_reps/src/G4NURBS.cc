// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBS.cc,v 1.1 1999-01-07 16:09:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// G4NURBS.cc
// Implementation of class G4NURBS
// OC 100796


#include "G4NURBS.hh"

// G4NURBS.hh includes globals.hh which includes a lot of others
// so no more includes required here

// stdlib required for the exit function
//#include <stdlib.h>

// memcpy
//#include <string.h>


////////////////////////////////////////////////////////////////////////
//    Here start the real world. Please, check your armored jacket.   //
////////////////////////////////////////////////////////////////////////



ostream & operator << (ostream & inout_outStream, const G4NURBS & in_kNurb)
		{
		inout_outStream
			// the magic could be changed for good reasons only
			<< "##ojc{NURBS}def[1.01.96.7]   Just a magic. Could be added to /etc/magic"
		  	<< "\n# NURBS Definition File (human and computer readable format)"
			<< "\n# :" << in_kNurb.Whoami()
			<< "\n# U order\tV order : " 
			<< '\n' << in_kNurb.GetUorder() << "\t\t" << in_kNurb.GetVorder();
		// number of knots and knots themselves for U and V
		for (G4NURBS::t_direction dir = G4NURBS::U; dir < G4NURBS::NofD; /*(*(int *)(&dir))++*/ dir=(G4NURBS::t_direction)(((int)(dir))+1) )
		 {
		 inout_outStream
			<< "\n# Number of knots along " << G4NURBS::Tochar(dir)
			<< '\n' << in_kNurb.GetnbrKnots(dir)
			<< "\n# " << G4NURBS::Tochar(dir) << " knots vector (as a column)";
		  {  // begin knots iteration
 		  G4double	oneKnot;
		  G4NURBS::KnotsIterator knotI(in_kNurb,dir);
		  G4bool otherKnots;
		  do 
		    {
		     otherKnots = knotI.pick(&oneKnot);
		     inout_outStream << "\n\t\t" << oneKnot;
		    }
		  while (otherKnots);
		  }  // end of knots iteration
		 };  // end of direction loop
		// number of control points in U and V direction
		// and controlpoints
		inout_outStream
			<< "\n# Number of control points along U and V"
			<< '\n' << in_kNurb.GetUnbrCtrlPts() 
			<< "   " << in_kNurb.GetVnbrCtrlPts()
			<< "\n# Control Points (one by line, U increasing first)";
		  { // begin of control points iteration
		  G4NURBS::t_doubleCtrlPt   oneCP;
		  G4NURBS::CtrlPtsIterator  cpI(in_kNurb);
		  G4bool otherCPs;
		  do
		    {
		     otherCPs = cpI.pick(&oneCP);
		     inout_outStream 
			<< "\n\t" << oneCP[G4NURBS::X]
			<< "\t" << oneCP[G4NURBS::Y]
			<< "\t" << oneCP[G4NURBS::Z]
			<< "\t" << oneCP[G4NURBS::W];
		    }
		  while (otherCPs);
		  } // end of control point iteration
		inout_outStream << "\n# That's all!" << endl;	// endl do an \n and a flush
		return inout_outStream;
		}

// the CC compiler issue some "maybe no value returned"
// but everything is ok

G4float	G4NURBS::GetfloatKnot(t_direction in_dir, t_indKnot in_index) const
	{
	in_dir = (t_direction)(in_dir & DMask);
	if ( in_index < m[in_dir].nbrKnots )
	return ((G4float)(m[in_dir].pKnots[in_index])); 
	else
	 {
	 G4cerr << "\nERROR: G4NURBS::GetfloatKnot: index out of range\n"
	      << "\n\t in_dir : " << in_dir << ", in_index : " << in_index
              << "m[in_dir].nbrKnots : " << m[in_dir].nbrKnots << endl;
	 return ((G4float)m[in_dir].pKnots[m[in_dir].nbrKnots-1]); 
	 };
	}
 
G4double	G4NURBS::GetdoubleKnot(t_direction in_dir, t_indKnot in_index) const
		{
		in_dir = (t_direction)(in_dir & DMask);
 		if ( in_index < m[in_dir].nbrKnots )
		return (G4double)(m[in_dir].pKnots[in_index]);
		else
		 {
		 G4cerr << "\nERROR: G4NURBS::GetdoubleKnot: index out of range"
		      << "\n\t in_dir : " << in_dir << ", in_index : " << in_index
	              << "m[in_dir].nbrKnots : " << m[in_dir].nbrKnots << endl;
		 return (G4double)(m[in_dir].pKnots[m[in_dir].nbrKnots-1]); 
		 };
		}

G4NURBS::t_floatCtrlPt*	G4NURBS::GetfloatCtrlPt(t_indCtrlPt in_onedimindex) const
		{
		if (in_onedimindex < mtotnbrCtrlPts)
  		return TofloatCtrlPt(mpCtrlPts[in_onedimindex]); 			
 		else
		 {
		 G4cerr 	<< "\nERROR: G4NURBS::GetfloatCtrlPt: index out of range"
			<< "\n\t in_onedimindex : " << in_onedimindex
			<< " , mtotnbrCtrlPts : " << mtotnbrCtrlPts << endl;
		 return TofloatCtrlPt(mpCtrlPts[mtotnbrCtrlPts-1]);  
		 };		
		}

G4NURBS::t_floatCtrlPt*	G4NURBS::GetfloatCtrlPt(t_inddCtrlPt in_Uindex, t_inddCtrlPt in_Vindex) const
		{
		if (  
			(in_Uindex < m[U].nbrCtrlPts)
  		     && (in_Vindex < m[V].nbrCtrlPts)
		   )
		return TofloatCtrlPt(mpCtrlPts[To1d(in_Uindex, in_Vindex)]);
		else
		 {
		 G4cerr 	<< "\nERROR: G4NURBS::GetfloatCtrlPt: index(s) out of range"
			<< "\n\t in_Uindex : " << in_Uindex
			<< " , in_Vindex : " << in_Vindex
			<< " , UnbrCtrlPts : " << m[U].nbrCtrlPts
			<< " , VnbrCtrlPts : " << m[V].nbrCtrlPts << endl;
		 return TofloatCtrlPt(mpCtrlPts[mtotnbrCtrlPts-1]);  
		 };		
		}

G4NURBS::t_doubleCtrlPt*	G4NURBS::GetdoubleCtrlPt(t_indCtrlPt in_onedimindex) const
		{
		if ( in_onedimindex < mtotnbrCtrlPts )
		return TodoubleCtrlPt(mpCtrlPts[in_onedimindex]);
		else
		 {
		 G4cerr 	<< "\nERROR: G4NURBS::getdoubleCtrlPts: index out of range"
			<< "\n\t in_onedimindex : " << in_onedimindex
			<< " , mtotnbrCtrlPts : " << mtotnbrCtrlPts << endl;
		 return TodoubleCtrlPt(mpCtrlPts[mtotnbrCtrlPts-1]);  
		 };		
		}

G4NURBS::t_doubleCtrlPt*	G4NURBS::GetdoubleCtrlPt(t_inddCtrlPt in_Uindex, t_inddCtrlPt in_Vindex) const
		{
		if (  
			(in_Uindex < m[U].nbrCtrlPts)
  		     && (in_Vindex < m[V].nbrCtrlPts)
		   )
		return TodoubleCtrlPt(mpCtrlPts[To1d(in_Uindex, in_Vindex)]);
		else
		 {
		 G4cerr 	<< "\nERROR: G4NURBS::GetdoubleCtrlPt: index(s) out of range"
			<< "\n\t in_Uindex : " << in_Uindex
			<< " , in_Vindex : " << in_Vindex
			<< " , UnbrCtrlPts : " << m[U].nbrCtrlPts
			<< " , VnbrCtrlPts : " << m[V].nbrCtrlPts << endl;
		 return TodoubleCtrlPt(mpCtrlPts[mtotnbrCtrlPts-1]);  
		 };		
		}
 
// Total copy
G4float *	G4NURBS::GetfloatAllKnots(t_direction in_dir) const
		{
		in_dir = (t_direction)(in_dir & DMask);
		G4float * p = new G4float [m[in_dir].nbrKnots];
		for (t_indKnot i = 0; i < m[in_dir].nbrKnots; i++)
			p[i] = (G4float)m[in_dir].pKnots[i];
		return p;
		}

G4double *	G4NURBS::GetdoubleAllKnots(t_direction in_dir) const
		{
		in_dir = (t_direction)(in_dir & DMask);
		G4double * p = new G4double [m[in_dir].nbrKnots];
		for (t_indKnot i = 0; i < m[in_dir].nbrKnots; i++)
			p[i] = (G4double)m[in_dir].pKnots[i];
		return p;
		}



G4float *	G4NURBS::GetfloatAllCtrlPts() const
		{
		G4float * p = new G4float [mtotnbrCtrlPts*NofC];
		for (t_indKnot i = 0; i < mtotnbrCtrlPts*NofC; i++)
			p[i] = (G4float)(((t_Coord *)mpCtrlPts)[i]);
		return p;
		}

G4double *	G4NURBS::GetdoubleAllCtrlPts() const
		{
		G4double * p = new G4double [mtotnbrCtrlPts*NofC];
		for (t_indKnot i = 0; i < mtotnbrCtrlPts*NofC; i++)
			p[i] = (G4double)(((t_Coord *)mpCtrlPts)[i]);
		return p;
		}



// Iterators

	G4NURBS::KnotsIterator::KnotsIterator(const G4NURBS & in_rNurb, G4NURBS::t_direction in_dir, t_indKnot in_startIndex)
		:	kmdir((G4NURBS::t_direction)(in_dir &  G4NURBS::DMask)),
			kmpMax(in_rNurb.m[kmdir].pKnots + in_rNurb.m[kmdir].nbrKnots)
		
		{
		if (in_startIndex < in_rNurb.m[kmdir].nbrKnots)
		mp = in_rNurb.m[kmdir].pKnots + in_startIndex;
		else
		  {
		  G4cerr 	<< "\nERROR: G4NURBS::KnotsIterator: in_startIndex out of range"
			<< "\n\tin_startIndex : " << in_startIndex
			<< ", nbr of knots : " << in_rNurb.m[kmdir].nbrKnots
			<< "\n\t mp set to NULL, calls to picking functions will fail"
			<< endl;
		  mp = NULL;
		  };
		}

G4bool	G4NURBS::KnotsIterator::pick(G4double * inout_pDbl)
		{
		(*inout_pDbl) = (G4double)(*mp);
		return (G4bool)((++mp)<kmpMax);
		}
		
G4bool	G4NURBS::KnotsIterator::pick(G4float * inout_pFlt)
		{
		(*inout_pFlt) = (G4float)(*mp);
		return (G4bool)((++mp)<kmpMax);
		}



	G4NURBS::CtrlPtsCoordsIterator::CtrlPtsCoordsIterator(const G4NURBS & in_rNurb, t_indCtrlPt in_startCtrlPtIndex)
		:	kmpMax((const t_Coord *)(in_rNurb.mpCtrlPts + in_rNurb.mtotnbrCtrlPts))
		{
 		if (in_startCtrlPtIndex < in_rNurb.mtotnbrCtrlPts )
		mp = (const t_Coord *)(in_rNurb.mpCtrlPts + in_startCtrlPtIndex);
		else
		  {
		  G4cerr 	<< "\nERROR: G4NURBS::CtrlPtsCoordsIterator: in_startCtrlPtIndex out of range"
			<< "\n\tin_startCtrlPtIndex : " << in_startCtrlPtIndex
			<< ", nbr of CtrlPts : " << in_rNurb.mtotnbrCtrlPts 
			<< "\n\t mp set to NULL, calls to picking functions will fail"
			<< endl;
		  mp = NULL;
		  };
		}

G4bool G4NURBS::CtrlPtsCoordsIterator::pick(G4double  * inout_pDbl)
		{
		(*inout_pDbl) = (G4double)((*mp));
		return (G4bool)((++mp)<kmpMax);
		}

G4bool G4NURBS::CtrlPtsCoordsIterator::pick(G4float * inout_pFlt)
		{
		(*inout_pFlt) = (G4float)((*mp));
		return (G4bool)((++mp)<kmpMax);
		}
		
	G4NURBS::CtrlPtsIterator::CtrlPtsIterator(const G4NURBS & in_rNurb, t_indCtrlPt in_startIndex)
		:	kmpMax(in_rNurb.mpCtrlPts + in_rNurb.mtotnbrCtrlPts)
		{
 		if (in_startIndex < in_rNurb.mtotnbrCtrlPts )
		mp = (in_rNurb.mpCtrlPts + in_startIndex);
		else
		  {
		  G4cerr 	<< "\nERROR: G4NURBS::CtrlPtsIterator: in_startIndex out of range"
			<< "\n\tin_startIndex : " << in_startIndex
			<< ", nbr of CtrlPts : " << in_rNurb.mtotnbrCtrlPts 
			<< "\n\t mp set to NULL, calls to picking functions will fail"
			<< endl;
		  mp = NULL;
		  };
		}

G4bool G4NURBS::CtrlPtsIterator::pick(t_doubleCtrlPt * inout_pDblCtrlPt)
                {
                for (t_indCoord i = G4NURBS::X; i < G4NURBS::NofC; i++)
                   (*inout_pDblCtrlPt)[i] = (G4double)((*mp)[i]);
                return (G4bool)((++mp)<kmpMax);
                }

G4bool G4NURBS::CtrlPtsIterator::pick(t_floatCtrlPt * inout_pFltCtrlPt)
                {
                for (t_indCoord i = G4NURBS::X; i < G4NURBS::NofC; i++)
                   (*inout_pFltCtrlPt)[i] = (G4float)((*mp)[i]);
                return (G4bool)((++mp)<kmpMax);
		}		


////////////////////////////////////////////////////////////////////////
// Building functions


G4bool	G4NURBS::MakeKnotVector(t_Dir & io_d, t_KnotVectorGenFlag in_KVGFlag)
	{
	G4bool isgood =    (io_d.order + io_d.nbrCtrlPts == io_d.nbrKnots)
                        && (io_d.pKnots == NULL);
	if ( isgood )
	{
	io_d.pKnots = new t_Knot [io_d.nbrKnots];
	if (in_KVGFlag != UserDefined)
	 	{  // let's do the knots
		t_indKnot indKnot = 0;
		t_index nbrCentralDistinctKnots = io_d.nbrCtrlPts-io_d.order;
		if ( (nbrCentralDistinctKnots % in_KVGFlag) == 0)
		{
		nbrCentralDistinctKnots /= in_KVGFlag; 
		// first and last knots repeated 'order' Times
		for (t_index i=0; i < io_d.order; indKnot++,i++)
		   {
		   io_d.pKnots[indKnot] = 0;
		   io_d.pKnots[indKnot+io_d.nbrCtrlPts] = 1;
		   };

		t_Knot stepKnot = 1.0/(t_Knot)(nbrCentralDistinctKnots+1);
		t_Knot valKnot = stepKnot;

		// central knots
		for (t_indKnot j=0; j < nbrCentralDistinctKnots; valKnot += stepKnot, j++)
		   {
		   for (t_indKnot k=0; k < in_KVGFlag; indKnot++, k++)
		      io_d.pKnots[indKnot] = valKnot;
		   };
		}
		else isgood = false;
		}; // end of knots making
	};	
	return isgood;
	}


ostream &	operator << (ostream & io_ostr, G4NURBS::t_KnotVectorGenFlag in_f)
		{
		switch (in_f) 
			{
			case G4NURBS::UserDefined: io_ostr << "UserDefined"; break;
			case G4NURBS::Regular:     io_ostr << "Regular"; break;
			case G4NURBS::RegularRep:  io_ostr << "RegularRep"; break;
			default:                   io_ostr << (int)in_f;
			};
		return io_ostr;
		}



////////////////////////////////////////////////////////////////////////
// Constructors and co

void	G4NURBS::Conscheck() const
		{
		G4int dummy;
		t_direction dir;
		for (dummy=0; (dummy?(dir=V):(dir=U)),(dummy < NofD); dummy++)
			{ 
			if (m[dir].order<=0)
				{ 
				G4cerr 	<< "\nFATAL ERROR: G4NURBS::G4NURBS: The order in the "
					<< G4NURBS::Tochar(dir) 
					<< " direction must be >= 1" << endl;
				exit(-1);
				};
			if (m[dir].nbrCtrlPts<=0)
				{
				G4cerr 	<< "\nFATAL ERROR: G4NURBS::G4NURBS: The number of control points "
					<< G4NURBS::Tochar(dir) 
					<< " direction must be >= 1" << endl;
				exit(-1);
				};
			};	// end of dummy
		}

	G4NURBS::G4NURBS
			(
			t_order in_Uorder, t_order in_Vorder,
			t_inddCtrlPt in_UnbrCtrlPts, t_inddCtrlPt in_VnbrCtrlPts,
			t_CtrlPt * in_pCtrlPts,
			t_Knot * in_pUKnots, t_Knot * in_pVKnots,
			t_CheckFlag in_CheckFlag
			)
		{

		m[U].order=in_Uorder; m[V].order=in_Vorder;
		m[U].nbrCtrlPts=in_UnbrCtrlPts; m[V].nbrCtrlPts=in_VnbrCtrlPts;

		mtotnbrCtrlPts = m[U].nbrCtrlPts * m[V].nbrCtrlPts;
		m[U].nbrKnots = m[U].order + m[U].nbrCtrlPts;
		m[V].nbrKnots = m[V].order + m[V].nbrCtrlPts;
		
		if (in_CheckFlag)
			Conscheck();

		// CtrlPts
		if (! (mpCtrlPts = in_pCtrlPts) )
			{
			G4cerr 	<< "\nFATAL ERROR: G4NURBS::G4NURBS: A NURBS MUST HAVE CONTROL POINTS!\n\teven if they are defined later, the array must be allocated."
				<< "\n\tgood bye. Have a nice debuging." << endl;
			exit(-1);
			};
		//mnbralias = 0;

		// Knots
		t_direction dir;
		G4int dummy;
		for (dummy=0; (dummy?(dir=V):(dir=U)),(dummy < NofD); dummy++)
		   {
		    if ( !(m[dir].pKnots = (dummy?in_pVKnots:in_pUKnots)) )
		  	{	// make some regular knots between 0 & 1
			if(!MakeKnotVector(m[dir], Regular))
				{
				G4cerr	<< "\nFATAL ERROR: G4NURBS::G4NURBS: Unable to make a Regular knot vector along "
					<< G4NURBS::Tochar(dir)
					<< " direction."
					<< "\n\tgood bye. Have a nice debuging."
					<< endl;
				exit(-1);

				};
			//m[dir].nbralias = 0;
			};	// end of knots-making
		   };// end for dummy

		} // end of G4NURBS::G4NURBS 

		
		// second constructor

	G4NURBS::G4NURBS(
			t_order in_Uorder, t_order in_Vorder,
			t_inddCtrlPt in_UnbrCtrlPts, t_inddCtrlPt in_VnbrCtrlPts,
			t_KnotVectorGenFlag in_UKVGFlag,
			t_KnotVectorGenFlag in_VKVGFlag,
			t_CheckFlag in_CheckFlag
			)
		{

		m[U].order=in_Uorder; m[V].order=in_Vorder;
		m[U].nbrCtrlPts=in_UnbrCtrlPts; m[V].nbrCtrlPts=in_VnbrCtrlPts;

		mtotnbrCtrlPts = m[U].nbrCtrlPts * m[V].nbrCtrlPts;
		m[U].nbrKnots = m[U].order + m[U].nbrCtrlPts;
		m[V].nbrKnots = m[V].order + m[V].nbrCtrlPts;
		
		if (in_CheckFlag)
			Conscheck();

		// Allocate CtrlPts
		mpCtrlPts = new t_CtrlPt [mtotnbrCtrlPts];
		//mnbralias = 0;

		// Knots
		t_direction dir;
		G4int dummy;
		for (dummy=0; (dummy?(dir=V):(dir=U)),(dummy < NofD); dummy++)
		   {
		    t_KnotVectorGenFlag flag = (dummy?in_VKVGFlag:in_UKVGFlag);
		    m[dir].pKnots = NULL;	// (allocation under our control)
		    if (  flag && !MakeKnotVector(m[dir], flag) )
			{
			G4cerr	<< "\nFATAL ERROR: G4NURBS::G4NURBS: Unable to make knot vector along "
				<< G4NURBS::Tochar(dir)
				<< " direction. (" << m[dir].nbrKnots 
				<< " knots requested for a " 
				<< flag 
				<< " knots vector)"
				<< "\n\tgood bye. Have a nice debuging."
				<< endl;
			exit(-1);
			};
		    //m[dir].nbralias = 0;
		   };
		}
			
			





	G4NURBS::G4NURBS(const G4NURBS & in_krNurb)
		{
		// we assume the in nurbs is ok

		// the number of CtrlPts can be copied straightly
		mtotnbrCtrlPts = in_krNurb.mtotnbrCtrlPts;

		// the main datas

		// but as m is an array of t_Dir and as t_Dir
		// is just a structure and not a class with a copy cons
		// whe need to duplicate the knots
		t_direction dir;
		G4int dummy;
		for (dummy=0; (dummy?(dir=V):(dir=U)),(dummy < NofD); dummy++)
			{
			// first we do a 'stupid' copy of m[dir]
			m[dir] = in_krNurb.m[dir];
			// but as m is an array of t_Dir and as t_Dir
			// is just a structure and not a class with a copy cons
			// whe need to duplicate the knots
			m[dir].pKnots = new G4double [m[dir].nbrKnots];
			// we copy the knots with memcpy. This function should be the fastest
			memcpy(m[dir].pKnots, in_krNurb.m[dir].pKnots, m[dir].nbrKnots * sizeof(G4double));
			};	// end of dummy loop

		// the control points
		// once again we need to do the copy
		mpCtrlPts = new t_CtrlPt [mtotnbrCtrlPts];
		memcpy(mpCtrlPts, in_krNurb.mpCtrlPts, mtotnbrCtrlPts*sizeof(t_CtrlPt)); 

		// and as it's very strange to copy a nurbs in G4
		// we issue a warning :
		G4cerr << "\nWARNING: G4NURBS::G4NURBS(const G4NURBS &) used" << endl;
		}	// end of G4NURBS::G4NURBS(const G4NURBS &)


	G4NURBS::~G4NURBS()
		{
		// we must free the two knots vector
		t_direction dir;
		G4int dummy;
		for (dummy=0; (dummy?(dir=V):(dir=U)),(dummy < NofD); dummy++)
			{
			if (m[dir].pKnots)
			  delete m[dir].pKnots;	// [m[dir].nbrKnots] if t_Knot become a class
			m[dir].pKnots = NULL;
			};
		// now we free the CtrlPts array
		if (mpCtrlPts)
		  delete [] mpCtrlPts;		// [mtotnbrCtrlPts] if t_CtrlPt become a class
		mpCtrlPts = NULL;
		}

/************************************************************************
 *                                                                      *
 * Return the current knot the parameter u is less than or equal to.    *
 * Find this "breakpoint" allows the evaluation routines to concentrate *
 * on only those control points actually effecting the curve around u.] *
 *                                                                      *
 *	m   is the number of points on the curve (or surface direction) *
 *	k   is the order of the curve (or surface direction)            *
 *	kv  is the knot vector ([0..m+k-1]) to find the break point in. *
 *                                                                      *
 ************************************************************************/
static int FindBreakPoint(double u, const Float *kv, int m, int k)
{
  int i;
  if (u == kv[m+1]) return m;      /* Special case for closed interval */
  i = m + k;
  while ((u < kv[i]) && (i > 0)) i--;
  return(i);
}

/************************************************************************
 *                                                                      *
 * Compute Bi,k(u), for i = 0..k.                                       *
 *  u        the parameter of the spline to find the basis functions for*
 *  brkPoint the start of the knot interval ("segment")                 *
 *  kv       the knot vector                                            *
 *  k        the order of the curve                                     *
 *  bvals    the array of returned basis values.                        *
 *                                                                      *
 * (From Bartels, Beatty & Barsky, p.387)                               *
 *                                                                      *
 ************************************************************************/
static void BasisFunctions(double u, int brkPoint,
                           const Float *kv, int k, double *bvals)
{
  int r, s, i;
  double omega;

  bvals[0] = 1.0;
  for (r=2; r <= k; r++) {
    i = brkPoint - r + 1;
    bvals[r-1] = 0.0;
    for (s=r-2; s >= 0; s--) {
      i++;
      if (i < 0) {
	omega = 0.0;
      }else{
	omega = (u - kv[i]) / (kv[i+r-1] - kv[i]);
      }
      bvals[s+1] = bvals[s+1] + (1.0-omega) * bvals[s];
      bvals[s]   = omega * bvals[s];
    }
  }
}

/************************************************************************
 *                                                                      *
 * Compute derivatives of the basis functions Bi,k(u)'                  *
 *                                                                      *
 ************************************************************************/
static void BasisDerivatives(double u, int brkPoint,
                             const Float *kv, int k, double *dvals)
{
  int s, i;
  double omega, knotScale;

  BasisFunctions(u, brkPoint, kv, k-1, dvals);

  dvals[k-1] = 0.0;	    /* BasisFunctions misses this */

  knotScale = kv[brkPoint+1] - kv[brkPoint];

  i = brkPoint - k + 1;
  for (s=k-2; s >= 0; s--) {
    i++;
    omega = knotScale * ((double)(k-1)) / (kv[i+k-1] - kv[i]);
    dvals[s+1] += -omega * dvals[s];
    dvals[s] *= omega;
  }
}

/***********************************************************************
 *                                                                     *
 * Calculate a point p on NurbSurface n at a specific u, v             *
 * using the tensor product.                                           *
 *                                                                     *
 * Note the valid parameter range for u and v is                       *
 * (kvU[orderU] <= u < kvU[numU), (kvV[orderV] <= v < kvV[numV])       *
 *                                                                     *
 ***********************************************************************/
void G4NURBS::CalcPoint(double u, double v, G4Point3D &p,
			G4Vector3D &utan, G4Vector3D &vtan) const
{
#define MAXORDER 50
  struct Point4 {
    double x, y, z, w;
  };

  int i, j, ri, rj;
  int ubrkPoint, ufirst;
  double bu[MAXORDER], buprime[MAXORDER];
  int vbrkPoint, vfirst;
  double bv[MAXORDER], bvprime[MAXORDER];
  Point4 r, rutan, rvtan;

  r.x = 0.0; r.y = 0.0; r.z = 0.0; r.w = 0.0;
  rutan = r;   rvtan = r;

  int numU   = GetUnbrCtrlPts();
  int numV   = GetVnbrCtrlPts();
  int orderU = GetUorder();
  int orderV = GetVorder();
  
  /* Evaluate non-uniform basis functions (and derivatives) */
  
  ubrkPoint = FindBreakPoint(u, m[U].pKnots, numU-1, orderU);
  ufirst    = ubrkPoint - orderU + 1;
  BasisFunctions  (u, ubrkPoint, m[U].pKnots, orderU, bu);
  BasisDerivatives(u, ubrkPoint, m[U].pKnots, orderU, buprime);

  vbrkPoint = FindBreakPoint(v, m[V].pKnots, numV-1, orderV);
  vfirst    = vbrkPoint - orderV + 1;
  BasisFunctions  (v, vbrkPoint, m[V].pKnots, orderV, bv);
  BasisDerivatives(v, vbrkPoint, m[V].pKnots, orderV, bvprime);

  /* Weight control points against the basis functions */

  t_doubleCtrlPt *cpoint;
  Point4 cp;
  double tmp;

  for (i=0; i<orderV; i++) {
    for (j=0; j<orderU; j++) {
      ri = orderV - 1 - i;
      rj = orderU - 1 - j;
	
      tmp = bu[rj] * bv[ri];
      cpoint = GetdoubleCtrlPt(j+ufirst, i+vfirst);
      cp.x = *cpoint[G4NURBS::X];
      cp.y = *cpoint[G4NURBS::Y];
      cp.z = *cpoint[G4NURBS::Z];
      cp.w = *cpoint[G4NURBS::W];
      r.x += cp.x * tmp;
      r.y += cp.y * tmp;
      r.z += cp.z * tmp;
      r.w += cp.w * tmp;
      
      tmp = buprime[rj] * bv[ri];
      rutan.x += cp.x * tmp;
      rutan.y += cp.y * tmp;
      rutan.z += cp.z * tmp;
      rutan.w += cp.w * tmp;

      tmp = bu[rj] * bvprime[ri];
      rvtan.x += cp.x * tmp;
      rvtan.y += cp.y * tmp;
      rvtan.z += cp.z * tmp;
      rvtan.w += cp.w * tmp;
    }
  }

  /* Project tangents, using the quotient rule for differentiation */
  
  double wsqrdiv = 1.0 / (r.w * r.w);

  utan.setX((r.w * rutan.x - rutan.w * r.x) * wsqrdiv);
  utan.setY((r.w * rutan.y - rutan.w * r.y) * wsqrdiv);
  utan.setZ((r.w * rutan.z - rutan.w * r.z) * wsqrdiv);

  vtan.setX((r.w * rvtan.x - rvtan.w * r.x) * wsqrdiv);
  vtan.setY((r.w * rvtan.y - rvtan.w * r.y) * wsqrdiv);
  vtan.setZ((r.w * rvtan.z - rvtan.w * r.z) * wsqrdiv);

  p.setX(r.x / r.w);
  p.setY(r.y / r.w);
  p.setZ(r.z / r.w);
}

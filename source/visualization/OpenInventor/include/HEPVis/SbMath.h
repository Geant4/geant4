#ifndef HEPVis_SbMath_h
#define HEPVis_SbMath_h

#include <math.h>
#ifndef M_PI
#define M_PI       3.1415926535897931160E0
#define M_PI_2     1.5707963267948965580E0  
#endif

#define SbMinimum(a,b) ((a)<(b)?a:b)
#define SbMaximum(a,b) ((a)>(b)?a:b)

#define FCOS(x)   ((float)cos((double)(x)))
#define FSIN(x)   ((float)sin((double)(x)))
#define FACOS(x)  ((float)acos((double)(x)))
#define FASIN(x)  ((float)asin((double)(x)))
#define FTAN(x)   ((float)tan((double)(x)))
#define FATAN(x)  ((float)atan((double)(x)))
#define FSQRT(x)  ((float)sqrt((double)(x)))
#define FPOW(x,y) ((float)pow((double)(x),(double)(y)))
#define FLOG(x)   ((float)log((double)(x)))
#define FLOG10(x) ((float)log10((double)(x)))
#define FFLOOR(x) ((float)floor((double)(x)))
#define FFABS(x)  ((float)fabs((double)(x)))
#define FCEIL(x)  ((float)ceil((double)(x)))

#endif

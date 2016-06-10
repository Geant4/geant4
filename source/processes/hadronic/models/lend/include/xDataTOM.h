/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef xDataTOM_h_included
#define xDataTOM_h_included

#include <statusMessageReporting.h>

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

typedef int xDataTOM_Int;

/* Note: xDataTOM_interpolationFlag_flat must be last for current logic in xDataTOM_interpolation.c to work. */
enum xDataTOM_interpolationFlag { xDataTOM_interpolationFlag_invalid, xDataTOM_interpolationFlag_linear, xDataTOM_interpolationFlag_log, 
    xDataTOM_interpolationFlag_byRegion, xDataTOM_interpolationFlag_flat };
enum xDataTOM_interpolationQualifier { xDataTOM_interpolationQualifier_invalid, xDataTOM_interpolationQualifier_dependent,
    xDataTOM_interpolationQualifier_none, xDataTOM_interpolationQualifier_unitBase, xDataTOM_interpolationQualifier_correspondingPoints };
enum xDataTOM_frame { xDataTOM_frame_invalid, xDataTOM_frame_lab, xDataTOM_frame_centerOfMass };
enum xDataTOM_subAxesType { xDataTOM_subAxesType_proxy, xDataTOM_subAxesType_intepolationAxes };
enum xDataTOM_KalbachMannType { xDataTOM_KalbachMannType_fr, xDataTOM_KalbachMannType_fra };

typedef struct xDataTOM_interpolation_s xDataTOM_interpolation;
typedef struct xDataTOM_axis_s xDataTOM_axis;
typedef struct xDataTOM_axes_s xDataTOM_axes;
typedef struct xDataTOM_subAxes_s xDataTOM_subAxes;

typedef struct xDataTOM_XYs_s xDataTOM_XYs;
typedef struct xDataTOM_regionsXYs_s xDataTOM_regionsXYs;
typedef struct xDataTOM_W_XYs_s xDataTOM_W_XYs;
typedef struct xDataTOM_V_W_XYs_s xDataTOM_V_W_XYs;

typedef struct xDataTOM_LegendreSeries_s xDataTOM_LegendreSeries;
typedef struct xDataTOM_W_XYs_LegendreSeries_s xDataTOM_W_XYs_LegendreSeries;
typedef struct xDataTOM_regionsW_XYs_LegendreSeries_s xDataTOM_regionsW_XYs_LegendreSeries;
typedef struct xDataTOM_V_W_XYs_LegendreSeries_s xDataTOM_V_W_XYs_LegendreSeries;
typedef struct xDataTOM_KalbachMannCoefficients_s xDataTOM_KalbachMannCoefficients;
typedef struct xDataTOM_KalbachMann_s xDataTOM_KalbachMann;
typedef struct xDataTOM_polynomial_s xDataTOM_polynomial;

typedef struct xDataTOM_xDataInfo_s xDataTOM_xDataInfo;

typedef struct xDataTOM_attribute_s xDataTOM_attribute;
typedef struct xDataTOM_attributionList_s xDataTOM_attributionList;
typedef struct xDataTOM_elementListItem_s xDataTOM_elementListItem;
typedef struct xDataTOM_elementList_s xDataTOM_elementList;
typedef struct xDataTOM_element_s xDataTOM_element;
typedef struct xDataTOM_TOM_s xDataTOM_TOM;

typedef int (*xDataTOM_sortElementFunc)( const void *, const void * );

struct xDataTOM_interpolation_s {
    enum xDataTOM_interpolationFlag independent, dependent;
    enum xDataTOM_interpolationQualifier qualifier;
};

struct xDataTOM_axis_s {
    int index;
    char *label;
    char *unit;
    xDataTOM_interpolation interpolation;
};

struct xDataTOM_axes_s {
    int numberOfAxes;
    xDataTOM_axis *axis;
};

struct xDataTOM_subAxes_s {
    enum xDataTOM_subAxesType type;
    int offset;
    xDataTOM_axes *axes;
    xDataTOM_interpolation interpolation;
};

struct xDataTOM_XYs_s {
    int index, length;
    double value, accuracy;
    xDataTOM_subAxes subAxes;
    double *data;
};

struct xDataTOM_regionsXYs_s {
    int length;
    xDataTOM_axes *axes;
    xDataTOM_XYs *XYs;
};

struct xDataTOM_W_XYs_s {
    int index, length;
    double value;
    xDataTOM_subAxes subAxes;
    xDataTOM_XYs *XYs;
};

struct xDataTOM_V_W_XYs_s {
    int length;
    xDataTOM_subAxes subAxes;
    xDataTOM_W_XYs *W_XYs;
};

struct xDataTOM_LegendreSeries_s {
    int index, length;
    double value;
    double *LegendreSeries;
};

struct xDataTOM_W_XYs_LegendreSeries_s {
    int index, length;
    double value;
    xDataTOM_subAxes subAxes;
    xDataTOM_LegendreSeries *LegendreSeries;
};

struct xDataTOM_regionsW_XYs_LegendreSeries_s {
    int length;
    xDataTOM_axes *axes;
    xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries;
};

struct xDataTOM_V_W_XYs_LegendreSeries_s {
    int length;
    xDataTOM_subAxes subAxes;
    xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries;
};

struct xDataTOM_KalbachMannCoefficients_s {
    int index, length;
    double value;
    double *coefficients;
};

struct xDataTOM_KalbachMann_s {
    enum xDataTOM_KalbachMannType type;
    int numberOfEnergies;
    xDataTOM_subAxes subAxes;
    xDataTOM_KalbachMannCoefficients *coefficients;
};

struct xDataTOM_polynomial_s {
    int length;
    xDataTOM_subAxes subAxes;
    double *coefficients;
};

struct xDataTOM_xDataInfo_s {
    const char *ID;
    xDataTOM_element *element;
    xDataTOM_axes axes;
    void *data;
};

struct xDataTOM_attribute_s {
    xDataTOM_attribute *next;
    char *name;
    char *value;
};

struct xDataTOM_attributionList_s {
    int number;
    xDataTOM_attribute *attributes;
};

struct xDataTOM_elementListItem_s {
    xDataTOM_element *element;
    const char *sortString;
};

struct xDataTOM_elementList_s {
    int n;
    xDataTOM_elementListItem *items;
};

struct xDataTOM_element_s {
    int ordinal;                                    /* Counting from 0. */
    int index;                                      /* Value from "index" attribute if present or -1. */
    xDataTOM_element *parent;
    xDataTOM_element *next;
    char *name;
    xDataTOM_attributionList attributes;
    int numberOfChildren;
    xDataTOM_element *children;
    xDataTOM_xDataInfo xDataInfo;
};

struct xDataTOM_TOM_s {
    char *fileName;
    char *realFileName;
    xDataTOM_element root;
};

/*
* Stuff in common/xDataTOM.c
*/
xDataTOM_TOM *xDataTOM_importFile( statusMessageReporting *smr, const char *fileName );
xDataTOM_TOM *xDataTOM_mallocTOM( statusMessageReporting *smr );
int xDataTOM_initializeTOM( statusMessageReporting *smr, xDataTOM_TOM *doc );
void *xDataTOM_freeTOM( statusMessageReporting *smr, xDataTOM_TOM **TOM );
int xDataTOM_setFileNameTOM( statusMessageReporting *smr, xDataTOM_TOM *doc, const char *fileName );
void xDataTOM_displayTree( statusMessageReporting *smr, xDataTOM_TOM *TOM, int printAttributes );

xDataTOM_element *xDataTOM_mallocElement( statusMessageReporting *smr, xDataTOM_element *parent, int ordinal, int index, char const *name );
void xDataTOM_freeElement( xDataTOM_element **element );
void xDataTOM_releaseElement( xDataTOM_element *element );
xDataTOM_element *xDataTOM_addElementInElement( statusMessageReporting *smr, xDataTOM_element *parent, int index, char const *name );
xDataTOM_element *xDataTOM_getDocumentsElement( xDataTOM_TOM *TOM );
xDataTOM_element *xDataTOME_getFirstElement( xDataTOM_element *element );
xDataTOM_element *xDataTOME_getNextElement( xDataTOM_element *element );
xDataTOM_element *xDataTOME_getOneElementByName( statusMessageReporting *smr, xDataTOM_element *element, char const *name, int required );
int xDataTOM_numberOfElementsByName( statusMessageReporting *smr, xDataTOM_element *element, char const *name );
int xDataTOME_addAttribute( statusMessageReporting *smr, xDataTOM_element *element, char const *name, char const *value );
char const *xDataTOM_getAttributesValueInElement( xDataTOM_element *element, char const *name );
int xDataTOME_copyAttributionList( statusMessageReporting *smr, xDataTOM_attributionList *desc, xDataTOM_element *element );
int xDataTOME_convertAttributeToInteger( statusMessageReporting *smr, xDataTOM_element *element, char const *name, int *n );
int xDataTOME_convertAttributeToDouble( statusMessageReporting *smr, xDataTOM_element *element, char const *name, double *d );
int xDataTOME_convertAttributeToDoubleWithUnit( statusMessageReporting *smr, xDataTOM_element *element, char const *name, double *d, char *unit );
int xDataTOME_getInterpolation( statusMessageReporting *smr, xDataTOM_element *element, int index, 
    enum xDataTOM_interpolationFlag *independent, enum xDataTOM_interpolationFlag *dependent, enum xDataTOM_interpolationQualifier *qualifier );

void xDataTOMAL_initial( statusMessageReporting *smr, xDataTOM_attributionList *attributes );
void xDataTOMAL_release( xDataTOM_attributionList *attributes );
int xDataTOMAL_addAttribute( statusMessageReporting *smr, xDataTOM_attributionList *attributes, char const *name, char const *value );
char const *xDataTOMAL_getAttributesValue( xDataTOM_attributionList *attributes, char const *name );
int xDataTOMAL_copyAttributionList( statusMessageReporting *smr, xDataTOM_attributionList *desc, xDataTOM_attributionList *src );
int xDataTOMAL_convertAttributeToInteger( statusMessageReporting *smr, xDataTOM_attributionList *attributes, char const *name, int *n );
int xDataTOMAL_convertAttributeToDouble( statusMessageReporting *smr, xDataTOM_attributionList *attributes, char const *name, double *d );

void *xData_initializeData( statusMessageReporting *smr, xDataTOM_element *TE, char const *ID, size_t size );
int xDataTOM_isXDataID( xDataTOM_element *TE, char const *ID );

/*
* Stuff in common/xDataTOMMisc.c
*/
char *xDataTOMMisc_getAbsPath( statusMessageReporting *smr, const char *fileName );
int xDataTOM_setMessageError_ReturnInt( int value, statusMessageReporting *smr, void *userData, const char *file, int line, int code, const char *fmt, ... );
xDataTOM_element *xDataTOM_getLinksElement( statusMessageReporting *smr, xDataTOM_element *element, char const *link );

#define xDataTOMMisc_allocateCopyString2( smr, s, forItem ) xDataTOMMisc_allocateCopyString( smr, s, forItem, __FILE__, __LINE__ )

/*
* Stuff in common/xDataTOM_interpolation.c
*/
int xDataTOM_interpolation_set( statusMessageReporting *smr, xDataTOM_interpolation *interpolation, enum xDataTOM_interpolationFlag independent, 
    enum xDataTOM_interpolationFlag dependent, enum xDataTOM_interpolationQualifier qualifier );
int xDataTOM_interpolation_setFromString( statusMessageReporting *smr, xDataTOM_interpolation *interpolation, char const *str );
int xDataTOM_interpolation_copy( statusMessageReporting *smr, xDataTOM_interpolation *desc, xDataTOM_interpolation *src );

/*
* Stuff in common/xDataTOM_axes.c
*/
int xDataTOM_axes_initialize( statusMessageReporting *smr, xDataTOM_axes *axes, int numberOfAxes );
int xDataTOM_axes_release( xDataTOM_axes *axes );
char const *xDataTOM_axes_getLabel( statusMessageReporting *smr, xDataTOM_axes *axes, int index );
char const *xDataTOM_axes_getUnit( statusMessageReporting *smr, xDataTOM_axes *axes, int index );
int xDataTOM_axes_getInterpolation( statusMessageReporting *smr, xDataTOM_axes *axes, int index, enum xDataTOM_interpolationFlag *independent,
        enum xDataTOM_interpolationFlag *dependent, enum xDataTOM_interpolationQualifier *qualifier );

int xDataTOM_subAxes_initialize( statusMessageReporting *smr, xDataTOM_subAxes *subAxes, enum xDataTOM_subAxesType type, int offset, 
    xDataTOM_axes *axes, xDataTOM_interpolation *interpolation );
int xDataTOM_subAxes_release( xDataTOM_subAxes *subAxes );
char const *xDataTOM_subAxes_getLabel( statusMessageReporting *smr, xDataTOM_subAxes *subAxes, int index );
char const *xDataTOM_subAxes_getUnit( statusMessageReporting *smr, xDataTOM_subAxes *subAxes, int index );

xDataTOM_axis *xDataTOM_axis_new( statusMessageReporting *smr, int index, char const *label, char const *unit, xDataTOM_interpolation *interpolation );
int xDataTOM_axis_initialize( statusMessageReporting *smr, xDataTOM_axis *axis, int index, char const *label, char const *unit, 
    xDataTOM_interpolation *interpolation );
xDataTOM_axis *xDataTOM_axis_release( statusMessageReporting *smr, xDataTOM_axis *axis );
enum xDataTOM_frame xDataTOM_axis_stringToFrame( statusMessageReporting *smr, char const *frame );
char const *xDataTOM_axis_frameToString( statusMessageReporting *smr, enum xDataTOM_frame frame );

/*
* Stuff in common/xDataTOM_XYs.c
*/
int xDataTOM_XYs_free( xDataTOM_xDataInfo *xDI );
int xDataTOM_XYs_release( xDataTOM_XYs *XYs );
int xDataTOM_XYs_getData( xDataTOM_XYs *XYs, double **data );
int xDataTOM_XYs_getDataFromXDataInfo( xDataTOM_xDataInfo *xDI, double **data );

/*
* Stuff in common/xDataTOM_regionsXYs.c
*/
int xDataTOM_regionsXYs_free( xDataTOM_xDataInfo *xDI );

/*
* Stuff in common/xDataTOM_W_XYs.c
*/
xDataTOM_W_XYs *xDataTOM_W_XYs_new( statusMessageReporting *smr, int index, int length, double value, xDataTOM_axes *axes, int axesOffset );
int xDataTOM_W_XYs_initialize( statusMessageReporting *smr, xDataTOM_W_XYs *W_XYs, int index, int length, double value, xDataTOM_axes *axes,
    int axesOffset );
xDataTOM_W_XYs *xDataTOM_W_XYs_free( xDataTOM_W_XYs *W_XYs );
int xDataTOM_W_XYs_freeFrom_xDataInfo( xDataTOM_xDataInfo *xDI );
int xDataTOM_W_XYs_release( xDataTOM_W_XYs *W_XYs );
xDataTOM_xDataInfo *xDataTOME_getXData( xDataTOM_element *TE );
void *xDataTOME_getXDataIfID( statusMessageReporting *smr, xDataTOM_element *TE, char const *ID );

/*
* Stuff in common/xDataTOM_V_W_XYs.c
*/
int xDataTOM_V_W_XYs_initialize( statusMessageReporting *smr, xDataTOM_V_W_XYs *V_W_XYs, int length, xDataTOM_axes *axes );
int xDataTOM_V_W_XYs_free( xDataTOM_xDataInfo *xDI );

/*
* Stuff in common/xDataTOM_LegendreSeries.c
*/
int xDataTOM_LegendreSeries_initialize( statusMessageReporting *smr, xDataTOM_LegendreSeries *LegendreSeries, int index, int length, double value );
int xDataTOM_LegendreSeries_release( xDataTOM_LegendreSeries *LegendreSeries );

/*
* Stuff in common/xDataTOM_W_XYs_LegendreSeries.c
*/
int xDataTOM_W_XYs_LegendreSeries_initialize( statusMessageReporting *smr, xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries, int index, 
        int length, double value, enum xDataTOM_subAxesType subAxesType, xDataTOM_axes *axes, xDataTOM_interpolation *interpolation );
int xDataTOM_W_XYs_LegendreSeries_free( xDataTOM_xDataInfo *xDI );
int xDataTOM_W_XYs_LegendreSeries_release( xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries );

/*
* Stuff in common/xDataTOM_regionsW_XYs_LegendreSeries.c
*/
int xDataTOM_regionsW_XYs_LegendreSeries_initialize( statusMessageReporting *smr, xDataTOM_regionsW_XYs_LegendreSeries *regionsW_XYs_LegendreSeries,
        int length, xDataTOM_axes *axes );
int xDataTOM_regionsW_XYs_LegendreSeries_free( xDataTOM_xDataInfo *xDI );
int xDataTOM_regionsW_XYs_LegendreSeries_release( xDataTOM_regionsW_XYs_LegendreSeries *regionsW_XYs_LegendreSeries );

/*
* Stuff in common/xDataTOM_V_W_XYs_LegendreSeries.c
*/
int xDataTOM_V_W_XYs_LegendreSeries_initialize( statusMessageReporting *smr, xDataTOM_V_W_XYs_LegendreSeries *V_W_XYs_LegendreSeries, 
        int length, xDataTOM_axes *axes );
int xDataTOM_V_W_XYs_LegendreSeries_free( xDataTOM_xDataInfo *xDI );

/*
* Stuff in common/xDataTOM_KalbachMann.c
*/
int xDataTOM_KalbachMann_initialize( statusMessageReporting *smr, xDataTOM_KalbachMann *KalbachMann, int length, xDataTOM_axes *axes );
int xDataTOM_KalbachMann_free( xDataTOM_xDataInfo *xDI );
int xDataTOM_KalbachMann_release( xDataTOM_KalbachMann *KalbachMann );

/*
* Stuff in common/xDataTOM_polynomial.c
*/
int xDataTOM_polynomial_initialize( statusMessageReporting *smr, xDataTOM_polynomial *polynomial, int length, xDataTOM_axes *axes );
int xDataTOM_polynomial_free( xDataTOM_xDataInfo *xDI );
int xDataTOM_polynomial_release( xDataTOM_polynomial *polynomial );
int xDataTOM_polynomial_getData( xDataTOM_polynomial *polynomial, double **data );
int xDataTOM_polynomial_getDataFromXDataInfo( xDataTOM_xDataInfo *xDI, double **data );

#if defined __cplusplus
    }
    }
#endif

#endif              /* End of xDataTOM_h_included. */

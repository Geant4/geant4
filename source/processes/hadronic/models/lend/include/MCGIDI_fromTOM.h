/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#ifndef MCGIDI_fromTOM_h_included
#define MCGIDI_fromTOM_h_included

#include <xDataTOM_importXML_private.h>
#include "MCGIDI.h"

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

ptwXYPoints *MCGIDI_fromTOM_XYs_to_ptwXYPoints_linear( statusMessageReporting *smr, xDataTOM_XYs *XYs, enum ptwXY_interpolation_e interpolation );
int MCGIDI_fromTOM_pdfsOfXGivenW( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_pdfsOfXGivenW *dists, ptwXYPoints *norms,
    char const *toUnits[3] );
int MCGIDI_fromTOM_pdfOfX( statusMessageReporting *smr, ptwXYPoints *pdfXY, MCGIDI_pdfOfX *dist, double *norm );
int MCGIDI_fromTOM_interpolation( statusMessageReporting *smr, xDataTOM_element *element, int index, enum ptwXY_interpolation_e *interpolation );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of MCGIDI_fromTOM_h_included. */

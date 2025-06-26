/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"

namespace GIDI {

static Vector collapseVector( Vector const &a_vector, std::vector<int> const &a_collapseIndices, 
        std::vector<double> const &a_weight, bool a_normalize );
static void multiGroupSetup( Transporting::MultiGroup const &a_boundaries, ptwXPoints **a_boundaries_xs,
                Transporting::Flux const &a_flux, ptwXYPoints **a_fluxes_xys, ptwXPoints **a_multiGroupFlux );

/* *********************************************************************************************************//**
 * Collapses a multi-group vector.
 *
 * @param a_vector                  [in]    The Vector to collapse.
 * @param a_settings                [in]    Specifies the uncollapsed and collapsed multi-group information and the flux.
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_temperature             [in]    The temperature of the flux to use when collapsing.
 *
 * @return                                  Returns the collapsed Vector.
 ***********************************************************************************************************/

Vector collapse( Vector const &a_vector, Transporting::Settings const &a_settings, Transporting::Particles const &a_particles, double a_temperature ) {

    Transporting::Particle const *projectile( a_particles.particle( a_settings.projectileID( ) ) );
    Transporting::ProcessedFlux const *flux( projectile->nearestProcessedFluxToTemperature( a_temperature ) );
    std::vector<double> const &multiGroupFlux( flux->multiGroupFlux( ) );
    std::vector<int> const &collapseIndices( projectile->collapseIndices( ) );

    return( collapseVector( a_vector, collapseIndices, multiGroupFlux, true ) );
}

/* *********************************************************************************************************//**
 * Collapses a multi-group vector.
 *
 * @param a_vector                  [in]    The Vector to collapse.
 * @param a_collapseIndices         [in]    Maps uncollapsed indices to collapsed indices.
 * @param a_weight                  [in]    The uncollapsed flux weighting.
 * @param a_normalize               [in]    If true, divide each collapsed value by it corresponding collapsed weight value.
 * @return                                  Returns the collapsed Vector.
 ***********************************************************************************************************/

static Vector collapseVector( Vector const &a_vector, std::vector<int> const &a_collapseIndices, 
        std::vector<double> const &a_weight, bool a_normalize ) {

    std::size_t n1( a_collapseIndices.size( ) - 1 );
    std::size_t index1( a_collapseIndices[0] );
    Vector vectorCollapsed( n1 );

    if( a_vector.size( ) > 0 ) {
        for( std::size_t i1 = 0; i1 < n1; ++i1 ) {
            std::size_t index2( a_collapseIndices[i1+1] );
            double fluxSum = 0;
            double valueSum = 0;

            for( std::size_t i2 = index1; i2 < index2; ++i2 ) {
                fluxSum += a_weight[i2];
                valueSum += a_weight[i2] * a_vector[i2];
            }
            if( a_normalize && ( fluxSum != 0 ) ) valueSum /= fluxSum;
            vectorCollapsed[i1] = valueSum;
            index1 = index2;
        }
    }

    return( vectorCollapsed );
}

/* *********************************************************************************************************//**
 * Collapses a multi-group matrix.
 *
 * @param a_matrix                  [in]    The Matrix to collapse.
 * @param a_settings                [in]    Specifies the uncollapsed and collapsed multi-group information and the flux.
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_temperature             [in]    The temperature of the flux to use when collapsing.
 * @param a_productID               [in]    Particle id of the outgoing particle.
 * @return                                  Returns the collapsed Matrix.
 ***********************************************************************************************************/

Matrix collapse( Matrix const &a_matrix, Transporting::Settings const &a_settings, Transporting::Particles const &a_particles, double a_temperature, std::string const &a_productID ) {

    if( a_matrix.size( ) == 0 ) return( a_settings.multiGroupZeroMatrix( a_particles, a_productID, true ) );

    Transporting::Particle const *projectile( a_particles.particle( a_settings.projectileID( ) ) );
    Transporting::ProcessedFlux const *flux( projectile->nearestProcessedFluxToTemperature( a_temperature ) );
    std::vector<double> const &multiGroupFlux( flux->multiGroupFlux( ) );
    std::vector<int> const &projectileCollapseIndices( projectile->collapseIndices( ) );

    Transporting::Particle const *product( a_particles.particle( a_productID ) );
    std::size_t n2 = product->numberOfGroups( );

    std::vector<int> productCollapseIndices( product->collapseIndices( ) );
    productCollapseIndices[0] = 0;
    productCollapseIndices[n2] = a_matrix[0].size( );

    std::vector<double> conservationWeight( a_matrix[0].size( ), 1. );
    if( product->conserve() == Transporting::Conserve::energyOut ) {
        std::vector<double> boundaries = product->fineMultiGroup().boundaries();
        for( std::size_t i1 = 0; i1 < boundaries.size() - 1; ++i1 ) {
            conservationWeight[i1] = 0.5 * (boundaries[i1] + boundaries[i1+1]);
        }
    }

    Matrix productCollapsed( 0, 0 );
    for( std::size_t i1 = 0; i1 < a_matrix.size( ); ++i1 ) {
        productCollapsed.push_back( collapseVector( a_matrix[i1], productCollapseIndices, conservationWeight, false ) );
    }

    Matrix productCollapsedTranspose = productCollapsed.transpose( );
    Matrix collapsedTranspose( 0, 0 );
    for( std::size_t i2 = 0; i2 < n2; ++i2 ) {
        collapsedTranspose.push_back( collapseVector( productCollapsedTranspose[i2], projectileCollapseIndices, multiGroupFlux, true ) );
    }

    if( product->conserve() == Transporting::Conserve::energyOut ) {
        double denominator = 1;
        std::vector<double> boundaries = product->multiGroup().boundaries();
        for( std::size_t i1 = 0; i1 < boundaries.size() - 1; ++i1 ) {
            denominator = 0.5 * (boundaries[i1] + boundaries[i1+1]);
            collapsedTranspose[i1] /= denominator;
        }
    }

    return( collapsedTranspose.transpose( ) );
}

/* *********************************************************************************************************//**
 * Transport correct a vector.
 *
 * @param a_vector                  [in]    The Vector to transport correct.
 * @param a_transportCorrection     [in]    The Vector that has the transport correction terms.
 * @return                                  Returns the collapsed Matrix.
 ***********************************************************************************************************/

Vector transportCorrect( Vector const &a_vector, Vector const &a_transportCorrection ) {

    return( a_vector - a_transportCorrection );
}

/* *********************************************************************************************************//**
 * Transport correct a Matrix.
 *
 * @param a_matrix                  [in]    The Matrix to transport correct.
 * @param a_transportCorrection     [in]    The Vector that has the transport correction terms.
 * @return                                  Returns the collapsed Matrix.
 ***********************************************************************************************************/

Matrix transportCorrect( Matrix const &a_matrix, Vector const &a_transportCorrection ) {

    std::size_t size = a_transportCorrection.size( );
    Matrix corrected( a_matrix );

    if( size == 0 ) return( corrected );
    if( a_matrix.size( ) == 0 ) {
        corrected = Matrix( size, size ); }
    else {
        if( size != a_matrix.size( ) ) throw Exception( "transportCorrect: matrix rows different than vector size." );
    }

    for( std::size_t index = 0; index < size; ++index ) corrected[index][index] -= a_transportCorrection[index];
    return( corrected );
}

/*! \class MultiGroupCalulationInformation
 * This class stores data as needed to multi-group data. Since the flux may be temperature dependent, an instance of this
 * should only be used for one temperature.
 */

/* *********************************************************************************************************//**
 * Constructor.
 *
 * @param   a_multiGroup                [in]    The multi-group boundaries.
 * @param   a_flux                      [in]    The flux weighting.
 ***********************************************************************************************************/

MultiGroupCalulationInformation::MultiGroupCalulationInformation( Transporting::MultiGroup const &a_multiGroup, Transporting::Flux const &a_flux ) :
        m_multiGroup( a_multiGroup ),
        m_flux( a_flux ),
        m_boundaries_xs( nullptr ),
        m_fluxes_xys( nullptr ),
        m_multiGroupFlux( nullptr ) {

    multiGroupSetup( m_multiGroup, &m_boundaries_xs, m_flux, &m_fluxes_xys, &m_multiGroupFlux );
}

/* *********************************************************************************************************//**
 * Destructor that frees allocated memory.
 ***********************************************************************************************************/

MultiGroupCalulationInformation::~MultiGroupCalulationInformation( ) {

    ptwX_free( m_boundaries_xs );
    ptwXY_free( m_fluxes_xys );
    ptwX_free( m_multiGroupFlux );
}

/* *********************************************************************************************************//**
 * Returns a flux weighted multi-group version of the function *a_function*.
 *
 * @param a_boundaries              [in]    List of multi-group boundaries.
 * @param a_function                [in]    Function to multi-group.
 * @param a_flux                    [in]    Flux to use for weighting.
 *
 * @return                                  Returns the multi-grouped Vector of *a_function*.
 ***********************************************************************************************************/

Vector multiGroupXYs1d( Transporting::MultiGroup const &a_boundaries, Functions::XYs1d const &a_function, Transporting::Flux const &a_flux ) {

    std::vector<double> const &boundaries = a_boundaries.boundaries( );
    ptwXPoints *boundaries_xs = ptwX_create( nullptr, boundaries.size( ), boundaries.size( ), &(boundaries[0]) );
    if( boundaries_xs == nullptr ) throw Exception( "GIDI::multiGroup: ptwX_create failed." );

    Transporting::Flux_order const &flux_order_0 = a_flux[0];
    double const *energies = flux_order_0.energies( );
    double const *fluxes = flux_order_0.fluxes( );
    ptwXYPoints *fluxes_xys = ptwXY_createFrom_Xs_Ys( nullptr, ptwXY_interpolationLinLin, ptwXY_interpolationToString( ptwXY_interpolationLinLin ), 
        12, 1e-3, flux_order_0.size( ), 10, flux_order_0.size( ), energies, fluxes, 0 );
    if( fluxes_xys == nullptr ) {
        ptwX_free( boundaries_xs );
        throw Exception( "GIDI::multiGroup: ptwXY_createFrom_Xs_Ys failed." );
    }

    ptwXPoints *multiGroupFlux = ptwXY_groupOneFunction( nullptr, fluxes_xys, boundaries_xs, ptwXY_group_normType_none, nullptr );
    if( multiGroupFlux == nullptr ) {
        ptwX_free( boundaries_xs );
        ptwXY_free( fluxes_xys );
        throw Exception( "GIDI::multiGroup: ptwXY_groupOneFunction failed." );
    }

    ptwXYPoints *ptwXY = ptwXY_clone2( nullptr, a_function.ptwXY( ) );
    ptwXPoints *groups = nullptr;
    if( ptwXY != nullptr ) {
        ptwXY_mutualifyDomains( nullptr, ptwXY, 1e-12, 1e-12, 1, fluxes_xys, 1e-12, 1e-12, 1 );
        groups = ptwXY_groupTwoFunctions( nullptr, ptwXY, fluxes_xys, boundaries_xs, ptwXY_group_normType_norm, multiGroupFlux );
    }
    ptwX_free( boundaries_xs );
    ptwXY_free( fluxes_xys );
    ptwX_free( multiGroupFlux );
    ptwXY_free( ptwXY );
    if( groups == nullptr ) throw Exception( "GIDI::multiGroup: ptwXY_groupTwoFunctions failed." );

    Vector vector( ptwX_length( nullptr, groups ), ptwX_getPointAtIndex( nullptr, groups, 0 ) );
    ptwX_free( groups );

    return( vector );
}

/* *********************************************************************************************************//**
 * This function returns a flux weighted multi-group version of the product of *a_function1* times * *a_function2*.
 * The caller owns the returned instance and is respondible for deleting it (i.e., freeing its memory when no longer needed).
 *
 * @param   a_multiGroupCalulationInformation   [in]    Store multi-group boundary and flux data used for multi-grouping.
 * @param   a_function1                         [in]    First function of the product.
 * @param   a_function2                         [in]    Second function of the product, generally a reaction's cross section.
 *
 * @return                                      Returns the multi-grouped Vector.
 ***********************************************************************************************************/

Vector *multiGroupTwoXYs1ds( MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, Functions::XYs1d const &a_function1, 
                Functions::XYs1d const &a_function2 ) {

    ptwXPoints *boundaries_xs = const_cast<ptwXPoints *>( a_multiGroupCalulationInformation.m_boundaries_xs );
    ptwXPoints *multiGroupFlux = const_cast<ptwXPoints * >( a_multiGroupCalulationInformation.m_multiGroupFlux );
    ptwXYPoints *ptwXY1 = nullptr, *ptwXY2 = nullptr, *fluxes_xys = nullptr;
    ptwXPoints *groups = nullptr;
    std::string errorMessage( "GIDI::multiGroupTwoXYs1ds: ptwXY_clone2 for a_function1 failed." );
    statusMessageReporting *smr = nullptr;

    ptwXY1 = ptwXY_clone2( smr, a_function1.ptwXY( ) );
    if( ptwXY1 != nullptr ) {
        ptwXY2 = ptwXY_clone2( smr, a_function2.ptwXY( ) );
        if( ptwXY2 == nullptr ) {
            errorMessage = "GIDI::multiGroupTwoXYs1ds: ptwXY_clone2 for a_function2 failed."; }
        else {
            double min1, min2, max1, max2;

            ptwXY_domainMin( smr, ptwXY1, &min1 );
            ptwXY_domainMax( smr, ptwXY1, &max1 );
            ptwXY_domainMin( smr, ptwXY2, &min2 );
            ptwXY_domainMax( smr, ptwXY2, &max2 );
            fluxes_xys = ptwXY_domainSlice( smr, const_cast<ptwXYPoints *>( a_multiGroupCalulationInformation.m_fluxes_xys ), 
                    std::max( min1, min2 ), std::min( max1, max2), 10, 1 );

            if( fluxes_xys == nullptr ) {
                errorMessage = "GIDI::multiGroupTwoXYs1ds: ptwXY_domainSlice for flux failed."; }
            else {
                ptwXY_mutualifyDomains( smr, ptwXY1,     1e-12, 1e-12, 1, ptwXY2,     1e-12, 1e-12, 1 );
                ptwXY_mutualifyDomains( smr, ptwXY1,     1e-12, 1e-12, 1, fluxes_xys, 1e-12, 1e-12, 1 );
                ptwXY_mutualifyDomains( smr, fluxes_xys, 1e-12, 1e-12, 1, ptwXY2,     1e-12, 1e-12, 1 );
                groups = ptwXY_groupThreeFunctions( smr, ptwXY1, ptwXY2, fluxes_xys, boundaries_xs, ptwXY_group_normType_norm, multiGroupFlux );
            }
        }
    }
    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    ptwXY_free( fluxes_xys );

    if( groups == nullptr ) {
        throw Exception( errorMessage );
    }

    Vector *vector = new Vector( ptwX_length( smr, groups ), ptwX_getPointAtIndex( smr, groups, 0 ) );
    ptwX_free( groups );

    return( vector );
}

/* *********************************************************************************************************//**
 * Setups *a_boundaries_xs*, *a_fluxes_xys* and *a_multiGroupFlux* as needed by multi-grouping functions. Calling code are
 * responsible for free-ing *a_boundaries_xs*, *a_fluxes_xys* and *a_multiGroupFlux*.
 *
 * @param a_boundaries              [in]    List of multi-group boundaries.
 * @param a_boundaries_xs           [out]   A ptwXPoints representation of *a_boundaries*.
 * @param a_flux                    [in]    Flux to use for weighting.
 * @param a_fluxes_xys              [out]   A ptwXYPoints representation of *a_flux*.
 * @param a_multiGroupFlux          [out]   A ptwYPoints multi-grouped representation of *a_fluxes_xys*.
 ***********************************************************************************************************/

static void multiGroupSetup( Transporting::MultiGroup const &a_boundaries, ptwXPoints **a_boundaries_xs,
                Transporting::Flux const &a_flux, ptwXYPoints **a_fluxes_xys, ptwXPoints **a_multiGroupFlux ) {

    std::vector<double> const &boundaries = a_boundaries.boundaries( );
    *a_boundaries_xs = ptwX_create( nullptr, boundaries.size( ), boundaries.size( ), &(boundaries[0]) );
    if( *a_boundaries_xs == nullptr ) throw Exception( "GIDI::multiGroup: ptwX_create failed." );

    Transporting::Flux_order const &flux_order_0 = a_flux[0];
    double const *energies = flux_order_0.energies( );
    double const *fluxes = flux_order_0.fluxes( );
    *a_fluxes_xys = ptwXY_createFrom_Xs_Ys( nullptr, ptwXY_interpolationLinLin, ptwXY_interpolationToString( ptwXY_interpolationLinLin ),
        12, 1e-3, flux_order_0.size( ), 10, flux_order_0.size( ), energies, fluxes, 0 );
    if( *a_fluxes_xys == nullptr ) {
        *a_boundaries_xs = ptwX_free( *a_boundaries_xs );
        throw Exception( "GIDI::multiGroup: ptwXY_createFrom_Xs_Ys failed." );
    }

    *a_multiGroupFlux = ptwXY_groupOneFunction( nullptr, *a_fluxes_xys, *a_boundaries_xs, ptwXY_group_normType_none, nullptr );
    if( *a_multiGroupFlux == nullptr ) {
        *a_boundaries_xs = ptwX_free( *a_boundaries_xs );
        *a_fluxes_xys = ptwXY_free( *a_fluxes_xys );
        throw Exception( "GIDI::multiGroup: ptwXY_groupOneFunction failed." );
    }
}

/* *********************************************************************************************************//**
 * This function finds a component's data that can be multi-grouped (generally a **Functions::XYs1d** instance, multi-groups it with
 * the boundaries and flux data in *a_multiGroupCalulationInformation* with weight *a_crossSection* and adds/replaces with style
 * label *a_heatedMultiGroupLabel*.
 *
 * @param   a_heatedMultiGroupLabel             [in]    The label of the style for the multi-group data being added.
 * @param   a_multiGroupCalulationInformation   [in]    Store multi-group boundary and flux data used for multi-grouping.
 * @param   a_component                         [in]    The Component whose data will be multi-grouped.
 * @param   a_weight                            [in]    An additional function to weight the data with. This is generally a reactions cross section.
 ***********************************************************************************************************/

void calculate1dMultiGroupDataInComponent( LUPI_maybeUnused ProtareSingle const *a_protare, std::string const &a_heatedMultiGroupLabel,
                MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, Component &a_component, Functions::XYs1d const &a_weight ) {

    Functions::XYs1d const *xys1d = static_cast<Functions::XYs1d const *>( a_component.findInstanceOfTypeInLineage( a_heatedMultiGroupLabel, GIDI_XYs1dChars ) );
    Functions::XYs1d const *xys1d2{ nullptr };

    if( xys1d == nullptr ) {
        Functions::Constant1d const *constand1d = 
                static_cast<Functions::Constant1d const *>( a_component.findInstanceOfTypeInLineage( a_heatedMultiGroupLabel, GIDI_constant1dChars ) );
        if( constand1d != nullptr ) {
            xys1d2 = GIDI::Functions::XYs1d::makeConstantXYs1d( GIDI::Axes( ), constand1d->domainMin( ), constand1d->domainMax( ), constand1d->value( ) ); }
        else {
            Functions::Regions1d const *regions1d = 
                    static_cast<Functions::Regions1d const *>( a_component.findInstanceOfTypeInLineage( a_heatedMultiGroupLabel, GIDI_regions1dChars ) );
            if( regions1d != nullptr ) {
                xys1d2 = regions1d->asXYs1d( true, 1e-4, 1e-6, 1e-6 ); }
            else {
                Functions::Branching1d const *branching1d = 
                        static_cast<Functions::Branching1d const *>( a_component.findInstanceOfTypeInLineage( a_heatedMultiGroupLabel, GIDI_branching1dChars ) );
                if( branching1d != nullptr ) {
                    xys1d2 = Functions::XYs1d::makeConstantXYs1d( GIDI::Axes( ), a_weight.domainMin( ), a_weight.domainMax( ), branching1d->multiplicity( ) ); }
                else {
                    Functions::Polynomial1d const *polynomial1d =
                            static_cast<Functions::Polynomial1d const *>( a_component.findInstanceOfTypeInLineage( a_heatedMultiGroupLabel, GIDI_polynomial1dChars ) );
                    if( polynomial1d != nullptr ) {
                        xys1d2 = polynomial1d->asXYs1d( true, 1e-4, 1e-6, 1e-6 ); }
                    else {
                        throw Exception( "calculate1dMultiGroupDataInComponent: from findInstanceOfTypeInLineage, no XYs1d, Constant1d, Regions1d, Branching1d or Polynomial1d form found in "
                                + a_component.toXLink( ) + "." );
                    }
                }
            }
        }
        xys1d = xys1d2;
    }

    Vector *vector = multiGroupTwoXYs1ds( a_multiGroupCalulationInformation, *xys1d, a_weight );
    Functions::Gridded1d *gridded1d = a_component.get<Functions::Gridded1d>( a_heatedMultiGroupLabel );
    gridded1d->setData( *vector );
    delete vector;
    delete xys1d2;
}

/* *********************************************************************************************************//**
 * This function 
 *
 * @param   a_multiGroupCalulationInformation   [in]    Store multi-group boundary and flux data used for multi-grouping.
 * @param   a_weight                            [in]    An additional function to weight the data with. This is generally a reactions cross section.
 * @param   a_evaluated                         [in]    This is the evaluated form of the fission energy released.
 * @param   a_gridded1d                         [in]    This is the current multi-grouped form whose data will be replace.
 ***********************************************************************************************************/

void calculate1dMultiGroupFissionEnergyRelease( MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, Functions::XYs1d const &a_weight,
                Functions::Function1dForm const *a_evaluated, Functions::Function1dForm *a_gridded1d ) {

    if( a_gridded1d == nullptr ) return;

    Functions::XYs1d const *xys1d = nullptr;
    Functions::Gridded1d *gridded1d = static_cast<Functions::Gridded1d *>( a_gridded1d );

    if( a_evaluated->moniker( ) == GIDI_XYs1dChars ) {
        xys1d = static_cast<Functions::XYs1d const *>( a_evaluated ); }
    else if( a_evaluated->moniker( ) == GIDI_polynomial1dChars ) {
        Functions::Polynomial1d const *polynomial1d = static_cast<Functions::Polynomial1d const *>( a_evaluated );
        xys1d = polynomial1d->asXYs1d( true, 1e-4, 1e-6, 1e-6 ); }
    else {
        throw Exception( "calculate1dMultiGroupFissionEnergyRelease: form not XYs1d or Polynomial1d: " + a_evaluated->toXLink( ) + "." );
    }

    Vector *vector = multiGroupTwoXYs1ds( a_multiGroupCalulationInformation, *xys1d, a_weight );
    gridded1d->setData( *vector );
    delete vector;

    if( a_evaluated != xys1d ) delete xys1d;
}

}

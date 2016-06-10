/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#define GIDI_USE_BDFLS 0

#ifndef GIDI_settings_hpp_included
#define GIDI_settings_hpp_included 1

#include <string>
#include <vector>
#include <map>

#include <ptwX.h>
#include <ptwXY.h>
#include <statusMessageReporting.h>

/* Disable Effective C++ warnings in GIDI header files. */
#if defined( __INTEL_COMPILER )
#pragma warning( push )

#if   __INTEL_COMPILER > 1399
#pragma warning( disable:2021 )
#elif __INTEL_COMPILER > 1199
#pragma warning( disable:2304 )
#endif

#endif

#define GIDI_settings_projectileEnergyMode_continuousEnergy 1
#define GIDI_settings_projectileEnergyMode_grouped ( 1 << 1 )
#define GIDI_settings_projectileEnergyMode_fixedGrid ( 1 << 2 )

class GIDI_settings_group {

    private:
        std::string mLabel;
        std::vector<double> mBoundaries;

    public:
        GIDI_settings_group( std::string const &label = "empty", int size = 0 );
        GIDI_settings_group( std::string const &label, int length, double const *values );
        GIDI_settings_group( std::string const &label, std::vector<double> const &boundaries );
        GIDI_settings_group( GIDI_settings_group const &group );
        ~GIDI_settings_group( );

        inline double operator[]( int const index ) const { return( mBoundaries[index] ); }
        inline int size( void ) const { return( (int) mBoundaries.size( ) ); }
        inline int getNumberOfGroups( void ) const { return( (int) ( mBoundaries.size( ) - 1 ) ); }
        inline double const *pointer( void ) const { return( &(mBoundaries[0]) ); }

        void setFromCDoubleArray( int length, double *values );
        inline std::string getLabel( ) const { return( mLabel ); }
        int getGroupIndexFromEnergy( double energy, bool encloseOutOfRange ) const;
        inline bool isLabel( std::string &label ) const { return( label == mLabel ); }
        void print( bool outline = false, int valuesPerLine = 10 ) const;

    private:
        void initialize( std::string const &label, int size, int length, double const *values );
};

#if GIDI_USE_BDFLS
#include <cbdfls.h>

class GIDI_settings_groups_from_bdfls {

    private:
        std::vector<GIDI_settings_group> mGroups;

    public:
        GIDI_settings_groups_from_bdfls( std::string const &fileName );
        GIDI_settings_groups_from_bdfls( char const *fileName );
        GIDI_settings_groups_from_bdfls( cbdfls_file const *bdfls );
        ~GIDI_settings_groups_from_bdfls( );

        GIDI_settings_group getViaGID( int gid ) const;
        std::vector<std::string> getLabels( void ) const;
        std::vector<int> getGIDs( void ) const;
        void print( bool outline = true, int valuesPerLine = 10 ) const;

    private:
        void initialize( char const *fileName );
        void initialize2( cbdfls_file const *bdfls );
};
#endif

/**
    This class stores the flux for one Legendre order (see class GIDI_settings_flux).
*/
class GIDI_settings_flux_order {

    private:
        int mOrder;                 /**< The Legendre order of the flux. */
        std::vector<double> mEnergies;   /**< List of flux energies. */
        std::vector<double> mFluxes;     /**< List of flux values - one for each element of mEnergies. */

    public:
        GIDI_settings_flux_order( int order     /**< The Legendre order for this flux data. */ );
        GIDI_settings_flux_order( int           order       /**< The Legendre order for this flux data. */,
                                  int           length      /**< The number or values in energies and fluxes. */,
                                  double const  *energies   /**< List of energies where flux is given. */,
                                  double const  *fluxes     /**< List of flux value for each energies value. */ );
        GIDI_settings_flux_order( int                   order       /**< The Legendre order for this flux data. */,
                                  std::vector<double> const  &energies   /**< List of energies where flux is given. */,
                                  std::vector<double> const  &fluxes     /**< List of flux value for each energies value. */ );
        GIDI_settings_flux_order( GIDI_settings_flux_order const &fluxOrder /**< Legendre flux order to copy. */ );
        ~GIDI_settings_flux_order( );

        inline int getOrder( void ) const { return( mOrder ); }
        inline int size( void ) const { return( (int) mEnergies.size( ) ); }
        inline double const *getEnergies( void ) const { return( &(mEnergies[0]) ); }
        inline double const *getFluxes( void ) const { return( &(mFluxes[0]) ); }
        void print( int valuesPerLine = 10 ) const;

    private:
        void initialize( int order, int length, double const *energies, double const *fluxes );
};

class GIDI_settings_flux {

    private:
        std::string mLabel;                                  /**< Label for the flux. */
        double mTemperature;
        std::vector<GIDI_settings_flux_order> mFluxOrders;   /**< List of fluxes for each Legendre order, l, sorted by Legendre order starting with l = 0. */

    public:
        GIDI_settings_flux( std::string const &label, double temperature_MeV );
        GIDI_settings_flux( char const *label, double temperature_MeV );
        GIDI_settings_flux( GIDI_settings_flux const &flux );
        ~GIDI_settings_flux( );

        GIDI_settings_flux_order const *operator[]( int order ) const;
        inline int getMaxOrder( void ) const { return( (int) mFluxOrders.size( ) - 1 ); }
        inline int size( void ) const { return( (int) mFluxOrders.size( ) ); }

        inline std::string getLabel( ) const { return( mLabel ); }
        inline bool isLabel( std::string const &label ) const { return( label == mLabel ); }
        inline bool isLabel( char const *label ) const { return( label == mLabel ); }
        inline double getTemperature( ) const { return( mTemperature ); }
        void addFluxOrder( GIDI_settings_flux_order const &fluxOrder );
        void print( bool outline = true, int valuesPerLine = 10 ) const;
};

#if GIDI_USE_BDFLS
class GIDI_settings_fluxes_from_bdfls {

    private:
        std::vector<GIDI_settings_flux> mFluxes;

    public:
        GIDI_settings_fluxes_from_bdfls( std::string const &fileName, double temperature_MeV );
        GIDI_settings_fluxes_from_bdfls( char const *fileName, double temperature_MeV );
        GIDI_settings_fluxes_from_bdfls( cbdfls_file const *bdfls, double temperature_MeV );
        ~GIDI_settings_fluxes_from_bdfls( );

        GIDI_settings_flux getViaFID( int fid );
        std::vector<std::string> getLabels( void );
        std::vector<int> getFIDs( void );
        void print( bool outline = true, int valuesPerLine = 10 );

    private:
        void initialize( char const *fileName, double temperature_MeV );
        void initialize2( cbdfls_file const *bdfls, double temperature_MeV );
};
#endif

class GIDI_settings_processedFlux {

    private:
        GIDI_settings_flux mFlux;
        std::vector<GIDI::ptwXYPoints *> mFluxXY;          /* Same as mFlux but stored as ptwXYPoints for each l-order. */
        std::vector<GIDI::ptwXPoints *> mGroupedFlux;      /* mFlux grouped using mGroupX, and stored as ptwXPoints for each l-order. */

    public:
        GIDI_settings_processedFlux( GIDI_settings_flux const &flux, GIDI::ptwXPoints *groupX );
        GIDI_settings_processedFlux( GIDI_settings_processedFlux const &flux );
        ~GIDI_settings_processedFlux( );

        inline double getTemperature( ) const { return( mFlux.getTemperature( ) ); }
        GIDI::ptwXPoints *groupFunction( GIDI::statusMessageReporting *smr, GIDI::ptwXPoints *groupX, GIDI::ptwXYPoints *ptwXY1, int order ) const;
};

class GIDI_settings_particle {

    private:
        int mPoPId;
        bool mTransporting;
        int mEnergyMode;
        GIDI_settings_group mGroup;
        GIDI::ptwXPoints *mGroupX;                    /* Same as mGroup but stored as ptwXPoints. */
        std::vector<GIDI_settings_processedFlux> mProcessedFluxes;

    public:
        GIDI_settings_particle( int PoPId, bool transporting, int energyMode );
        GIDI_settings_particle( GIDI_settings_particle const &particle );
        int initialize( int PoPId, bool transporting, int energyMode );
        ~GIDI_settings_particle( );

        int addFlux( GIDI::statusMessageReporting *smr, GIDI_settings_flux const &flux );
        GIDI_settings_processedFlux const *nearestFluxToTemperature( double temperature ) const;
        inline int getGroupIndexFromEnergy( double e_in, bool encloseOutOfRange ) const { return( mGroup.getGroupIndexFromEnergy( e_in, encloseOutOfRange ) ); };
        inline int getNumberOfGroups( void ) const { return( mGroup.getNumberOfGroups( ) ); };
        inline int getPoPId( void ) const { return( mPoPId ); }
        inline int getEnergyMode( void ) const { return( mEnergyMode ); }
        inline bool getTransporting( void ) const { return( mTransporting ); }
        inline GIDI_settings_group getGroup( void ) const { return( mGroup ); }
        GIDI_settings_flux const *getFlux( double temperature ) const;
        GIDI::ptwXPoints *groupFunction( GIDI::statusMessageReporting *smr, GIDI::ptwXYPoints *ptwXY1, double temperature, int order ) const;
        void setGroup( GIDI_settings_group const &group );

        inline bool isEnergyMode_continuous( void ) const { return( this->mEnergyMode & GIDI_settings_projectileEnergyMode_continuousEnergy ); }
        inline bool isEnergyMode_grouped( void ) const { return( this->mEnergyMode & GIDI_settings_projectileEnergyMode_grouped ); }
        inline bool isEnergyMode_fixedGrid( void ) const { return( this->mEnergyMode & GIDI_settings_projectileEnergyMode_fixedGrid ); }

    private:
        GIDI_settings_flux const *getProcessedFlux( double temperature ) const;
};

class GIDI_settings {

    private:
        std::map<int, GIDI_settings_particle> mParticles;

    public:
        GIDI_settings( );
        ~GIDI_settings( );

        int addParticle( GIDI_settings_particle const &particle );
        GIDI_settings_particle const *getParticle( int PoPId ) const;
        int eraseParticle( int PoPId );
        void releaseMemory( ) { mParticles.clear( ); }
};

#if defined( __INTEL_COMPILER )
#pragma warning( pop )
#endif

#endif      // End of GIDI_settings_hpp_included

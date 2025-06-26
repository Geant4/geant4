/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef RISI_hpp_included
#define RISI_hpp_included 1

#include <map>
#include <set>

#include <LUPI.hpp>

namespace GIDI {

namespace RISI {

class Projectile;

class Reaction {

    public:
        double m_effectiveThreshold;                        /**< The effective threshold for the reaction. */
        std::vector<std::string> m_products;                /**< The list of final products for the reaction. */
        std::vector<int> m_multiplicities;                  /**< The multiplicities for each product in *m_products*. */
        std::vector<std::string> m_intermediates;           /**< The list of intermediates products for the reaction. */
        std::string m_process;                              /**< The process for the reaction. */
        std::string m_reactionLabel;                        /**< The label of the reaction. */
        std::string m_convarianceFlag;                      /**< A flag indicating if covariance data are present for the reaction. */

    public:
        Reaction( double a_effectiveThreshold, std::vector<std::string> const &a_products, std::vector<int> const &a_multiplicities, 
                std::vector<std::string> const &a_intermediates, std::string const &a_process, std::string const &reactionLabel,
                std::string const &convarianceFlag );

        void products( double a_energyMax, std::set<std::string> &a_products ) const ;
};

class Protare {

    private:
        int m_addMode;                                      /**< Indicates which method **add** calls. */
        std::string m_projectile;                           /**< The PoPs id for the projectile. */
        std::string m_target;                               /**< The PoPs id for the target. */
        std::string m_evaluation;                           /**< The evaluation for the protare. */
        double m_energyConversionFactor;                    /**< Factor to convert from file energy units to user energy units. */

        std::map<std::string, std::string> m_aliases;       /**< The list of meta-stable aliases in the protare. */
        std::vector<Reaction *> m_reactions;                /**< The list of **Reaction** instances for the protare. */

    public:
        Protare( std::string const &a_projectile, std::string const &a_target, std::string const &a_evaluation, 
                std::string const &a_protareEnergyUnit, std::string const &a_requestedEnergyUnit );
        ~Protare();

        std::string const &projectile( ) { return( m_projectile ); }
        std::string const &target( ) { return( m_target ); }
        std::string const &evaluation( ) { return( m_evaluation ); }

        void Oops( std::vector<std::string> const &a_elements );
        void addAlias( std::vector<std::string> const &a_elements );
        void setAddingAliases( ) { m_addMode = 1; }         /**< Tells **add** method to call the **addAlias** method. */
        void addReaction( std::vector<std::string> const &a_elements );
        void setAddingReactions( ) { m_addMode = 2; }       /**< Tells **add** method to call the **addReaction** method. */
        void add( std::vector<std::string> const &a_elements );

        void products( Projectile const *a_projectile, int a_level, int a_maxLevel, double a_energyMax, std::map<std::string, int> &a_products ) const ;
};

class Target {

    private:
        std::string m_id;
        std::vector<Protare *> m_protares;

    public:
        Target( std::string const &a_id ) :
                m_id( a_id ) {
        }
        ~Target( );

        void add( Protare *a_protare );
        void products( Projectile const *a_projectile, int a_level, int a_maxLevel, double a_energyMax, std::map<std::string, int> &a_products ) const ;
        void print( std::string const &a_indent = "" ) const ;
};

class Projectile {

    private:
        std::string m_id;
        std::map<std::string, Target *> m_targets;

    public:
        Projectile( std::string const &a_id ) :
                m_id( a_id ) {
        }
        ~Projectile( );

        void add( Protare *a_protare );
        void products( std::string const &a_target, int a_level, int a_maxLevel, double a_energyMax, std::map<std::string, int> &a_products ) const ;
        void print( std::string const &a_indent = "" ) const ;
};

class Projectiles {

    private:
        std::map<std::string, Projectile *> m_projectiles;

    public:
        Projectiles( ) {}
        ~Projectiles( );

        void add( Protare *a_protare );
        void clear( );
        std::vector<std::string> products( std::string const &a_projectile, std::vector<std::string> const &a_seedTargets, int a_maxLevel, 
                double a_energyMax ) const ;
        void print( std::string const &a_indent = "" ) const ;
};

void readRIS( std::string const &a_fileName, std::string const &a_energyUnit, Projectiles &a_projectiles );

}           // End of namespace RISI.

}           // End of namespace GIDI.

#endif      // End of RISI_hpp_included

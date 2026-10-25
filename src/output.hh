// This file is part of MUSIC
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2010-2024 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __OUTPUT_HH
#define __OUTPUT_HH

#include <algorithm>
#include <stdexcept>
#include <string>
#include <map>

#include "general.hh"
#include "mesh.hh"

/*!
 * @class particle_load_mask
 * @brief masked simple cubic particle load on the finest level (as 'ParticleLoad = masked' in monofonIC)
 *
 * Every 2x2x2 block of finest-level cells carries one gas particle at the block centre (the cell
 * opposite the corner) and DM at the corner, plus the edge centres (ParticleMaskType bit 0) and/or
 * the face centres (bit 1). Remaining cells are left empty. Without baryons the load is plain SC.
 * Indices are local to the finest grid, which the refinement mask always fills in aligned 2x2x2 blocks.
 */
class particle_load_mask
{
	int mask_[8]; // 0 = DM, 1 = gas, -1 = empty

	int value(size_t i, size_t j, size_t k) const
	{
		return mask_[((i % 2) * 2 + j % 2) * 2 + k % 2];
	}

public:
	bool requested; //!< ParticleLoad = masked
	bool enabled;	//!< requested and baryons enabled

	explicit particle_load_mask(config_file &cf)
	{
		const std::string load = cf.get_value_safe<std::string>("setup", "ParticleLoad", "sc");
		// other formats never read ParticleLoad, so monofonIC-style values are only rejected for swift
		if (load != "sc" && load != "masked" && cf.get_value<std::string>("output", "format") == "swift")
			throw std::runtime_error("ParticleLoad must be 'sc' or 'masked' for swift output, got '" + load + "'");
		requested = load == "masked";
		enabled = requested && cf.get_value_safe<bool>("setup", "baryons", false);

		unsigned mask_type = requested ? cf.get_value_safe<unsigned>("setup", "ParticleMaskType", 3) : 3;
		if (mask_type > 3)
			throw std::runtime_error("ParticleMaskType must be 0, 1, 2 or 3");

		for (auto &m : mask_)
			m = -1;
		mask_[0] = 0; // DM at corner
		if (mask_type & 1)
			mask_[1] = mask_[2] = mask_[4] = 0; // DM at edge centres
		if (mask_type & 2)
			mask_[3] = mask_[5] = mask_[6] = 0; // DM at face centres
		mask_[7] = 1; // gas at block centre
	}

	bool is_dm(size_t i, size_t j, size_t k) const { return !enabled || value(i, j, k) == 0; }
	bool is_gas(size_t i, size_t j, size_t k) const { return !enabled || value(i, j, k) == 1; }

	//! fraction of finest-level cells holding a DM / gas particle
	double dm_fraction() const { return enabled ? std::count(mask_, mask_ + 8, 0) / 8.0 : 1.0; }
	double gas_fraction() const { return enabled ? 1.0 / 8.0 : 1.0; }
};


/*!
 * @class output_plugin
 * @brief abstract base class for output plug-ins
 *
 * This class provides the abstract base class for all output plug-ins.
 * All output plug-ins need to derive from it and implement the purely
 * virtual member functions.
 */
class output_plugin
{
protected:
	
	//! reference to the config_file object that holds all configuration options
	config_file& cf_;
	
	//! output file or directory name
	std::string fname_;
	
	//! minimum refinement level
	unsigned levelmin_;
	
	//! maximum refinement level
	unsigned levelmax_;
	
	std::vector<unsigned> 
		offx_,		//!< vector describing the x-offset of each level
		offy_,		//!< vector describing the y-offset of each level
		offz_,		//!< vector describing the z-offset of each level
		sizex_,		//!< vector describing the extent in x of each level
		sizey_,		//!< vector describing the extent in y of each level
		sizez_;		//!< vector describing the extent in z of each level
	
	//! quick access function to query properties of the refinement grid from the configuration options
	/*! @param name	name of the config property
	 *  @param icomp component index (0=x, 1=y, 2=z)
	 *  @param oit output iterator (e.g. std::back_inserter for vectors)
	 */
	template< typename output_iterator >
	void query_grid_prop( std::string name, int icomp, output_iterator oit )
	{
		char str[128];
		for( unsigned i=levelmin_; i<=levelmax_; ++i )
		{
			snprintf( str, 128, "%s(%u,%d)", name.c_str(), i, icomp );
			*oit = cf_.get_value<unsigned>( "setup", str );
			++oit;
		}
	}
	
public:
	
	//! constructor
	explicit output_plugin( config_file& cf )
	: cf_(cf)
	{ 
		fname_		= cf.get_value<std::string>("output","filename");
		levelmin_	= cf.get_value<unsigned>( "setup", "levelmin" );
		levelmax_	= cf.get_value<unsigned>( "setup", "levelmax" );

		query_grid_prop( "offset", 0, std::back_inserter(offx_) );
		query_grid_prop( "offset", 1, std::back_inserter(offy_) );
		query_grid_prop( "offset", 2, std::back_inserter(offz_) );
		
		query_grid_prop( "size", 0, std::back_inserter(sizex_) );
		query_grid_prop( "size", 1, std::back_inserter(sizey_) );
		query_grid_prop( "size", 2, std::back_inserter(sizez_) );
	}
	
	//! destructor
	virtual ~output_plugin()
	{ }
	
	//! purely virtual prototype to write the masses for each dark matter particle
	virtual void write_dm_mass( const grid_hierarchy& gh ) = 0;
	
	//! purely virtual prototype to write the dark matter density field 
	virtual void write_dm_density( const grid_hierarchy& gh ) = 0;
	
	//! purely virtual prototype to write the dark matter gravitational potential (from which displacements are computed in 1LPT)
	virtual void write_dm_potential( const grid_hierarchy& gh ) = 0;
	
	//! purely virtual prototype to write dark matter particle velocities
	virtual void write_dm_velocity( int coord, const grid_hierarchy& gh ) = 0;
	
	//! purely virtual prototype to write dark matter particle positions
	virtual void write_dm_position( int coord, const grid_hierarchy& gh ) = 0;

	//! purely virtual prototype to write the baryon velocities
	virtual void write_gas_velocity( int coord, const grid_hierarchy& gh ) = 0;

	//! purely virtual prototype to write the baryon coordinates
	virtual void write_gas_position( int coord, const grid_hierarchy& gh ) = 0;
	
	//! purely virtual prototype to write the baryon density field
	virtual void write_gas_density( const grid_hierarchy& gh )  = 0;
	
	//! purely virtual prototype to write the baryon gravitational potential (from which displacements are computed in 1LPT)
	virtual void write_gas_potential( const grid_hierarchy& gh )  = 0;
	
	//! purely virtual prototype for all things to be done at the very end
	virtual void finalize( void ) = 0;
};

/*!
 * @brief implements abstract factory design pattern for output plug-ins
 */
struct output_plugin_creator
{
	//! create an instance of a plug-in
	virtual output_plugin * create( config_file& cf ) const = 0;
	
	//! destroy an instance of a plug-in
	virtual ~output_plugin_creator() { }
};

//! maps the name of a plug-in to a pointer of the factory pattern 
std::map< std::string, output_plugin_creator *>& get_output_plugin_map();

//! print a list of all registered output plug-ins
void print_output_plugins();

/*!
 * @brief concrete factory pattern for output plug-ins
 */
template< class Derived >
struct output_plugin_creator_concrete : public output_plugin_creator
{
	//! register the plug-in by its name
	output_plugin_creator_concrete( const std::string& plugin_name )
	{
		get_output_plugin_map()[ plugin_name ] = this;
	}
	
	//! create an instance of the plug-in
	output_plugin * create( config_file& cf ) const
	{
		return new Derived( cf );
	}
};

//! failsafe version to select the output plug-in
output_plugin *select_output_plugin( config_file& cf );

#endif // __OUTPUT_HH

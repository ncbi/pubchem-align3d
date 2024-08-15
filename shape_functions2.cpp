/*  $Id: shape_functions2.cpp 685605 2024-07-26 12:29:33Z thiessen $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Authors:  Evan Bolton, Leonid Zaslavsky, Paul Thiessen
*
* ===========================================================================
*/

// C++ Includes

#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <iostream>
#include <algorithm>


#include "shape_exp_lookup.hpp"
#include "shape_exp_lookup_trunc.hpp"
#include "shape_math_functions.hpp"
#include "shape_functions.hpp"
#include "shape_constants.hpp"
#include "shape_debug.hpp"


namespace Align3D {

void getJointColorTypeSet( const unsigned* ref_atom_types,
				     const unsigned ref_natom, 
				     const unsigned* fit_atom_types,
				     const unsigned fit_natom,
				     std::set<unsigned>& joint_color_atom_type_set )
{
  std::set<unsigned> ref_color_atom_type_set;
  for( unsigned k=0; k < ref_natom; ++k ){
    unsigned atom_type=ref_atom_types[k]; 
    if( atom_type > 0 ){
      ref_color_atom_type_set.insert( atom_type );
    }
  }

  joint_color_atom_type_set.clear();
  for( unsigned k=0; k < fit_natom; ++k ){
    unsigned atom_type = fit_atom_types[k];
    if( atom_type > 0 ){
      if( ref_color_atom_type_set.find( atom_type ) != ref_color_atom_type_set.end() ){
	joint_color_atom_type_set.insert( atom_type );
      }
    }
  }
}


void getVolumeAtomIndexVector( const unsigned* atom_types,
					 const unsigned natom,
					 std::vector< unsigned >& atomIndexVector )
{
  for( unsigned k=0; k<natom; ++k ){
    unsigned atom_type = atom_types[ k ];
    if( atom_type == 0 ){
      atomIndexVector.push_back( k );
    }
  }
}




void getColorAtomType2IndexVectorMap( const unsigned* atom_types,
						const unsigned natom,
						std::map< unsigned, std::vector< unsigned > >& atom_type_to_index_vector_map )
{
  for( unsigned k=0; k<natom; ++k ){
    unsigned atom_type = atom_types[ k ];
    if( atom_type > 0 ){
      std::map< unsigned, std::vector< unsigned > >::iterator ptr = atom_type_to_index_vector_map.find( atom_type );
      if( ptr == atom_type_to_index_vector_map.end() ){
	atom_type_to_index_vector_map[ atom_type] = std::vector<unsigned>();
	ptr = atom_type_to_index_vector_map.find( atom_type );
      }
      ptr->second.push_back( k );
    }
  }
}

void restrictColorAtomType2IndexVectorMap( 
						    std::map< unsigned, std::vector< unsigned > >& atom_type_to_index_vector_map,
						    const std::set<unsigned>& joint_color_atom_type_set ){

  std::set< unsigned > delete_set;
  for( std::map< unsigned, std::vector< unsigned > >::const_iterator iter = atom_type_to_index_vector_map.begin(); iter != atom_type_to_index_vector_map.end(); ++iter ){
    unsigned atom_type = iter->first;
    if( joint_color_atom_type_set.find( atom_type ) == joint_color_atom_type_set.end() ){
      delete_set.insert( atom_type );
    }
  }

  for( unsigned atom_type : delete_set ){
    atom_type_to_index_vector_map.erase( atom_type );
  }
}





void ComputeOverlapAndAnalyticGradient(
						 double& o,                                                  // overlap
						 double* g,                                                  // gradient[ 7 ]
						 const double* m,                                            // rotation/translation matrix [ 12 ]
						 const float* ref_coord,                                     // per ref atom, XYZ coordinates [3*nref]
						 const std::vector<double>& alpha_ref_vector,
						 const std::vector< unsigned >& ref_volume_atom_index_vector,
						 const float* fit_coord,                                     // original per fit atom, XYZ coordinates [3*nfit]
						 const std::vector<double>& alpha_fit_vector,
						 const std::vector< unsigned >& fit_volume_atom_index_vector,
						 bool needGrad
						 )
{
  // Initialize shape volume overlap and its gradient
  o = 0.0;
  if( needGrad ){
    memset( g, 0, 7 * sizeof( double ) );
  }
  
  // Loop over atoms
  for( unsigned i : fit_volume_atom_index_vector ){
    // Transform original coordinates
    const double  old_x = fit_coord[ ( i * 3 )     ];
    const double  old_y = fit_coord[ ( i * 3 ) + 1 ];
    const double  old_z = fit_coord[ ( i * 3 ) + 2 ];

    double x, y, z;
    if( m == NULL ){
      x = old_x;
      y = old_y;
      z = old_z;
    }else{ 
      x = m[  9 ] + m[ 0 ] * old_x + m[ 1 ] * old_y + m[ 2 ] * old_z;
      y = m[ 10 ] + m[ 3 ] * old_x + m[ 4 ] * old_y + m[ 5 ] * old_z;
      z = m[ 11 ] + m[ 6 ] * old_x + m[ 7 ] * old_y + m[ 8 ] * old_z;
    }
    
    // Since we are applying the derivative on transformed coordinates, quaternion is (1, 0, 0, 0).
    
      
    for( unsigned j : ref_volume_atom_index_vector ){
      double x_diff = ref_coord[ ( j * 3 )     ] - x;
      double y_diff = ref_coord[ ( j * 3 ) + 1 ] - y;
      double z_diff = ref_coord[ ( j * 3 ) + 2 ] - z;
	
      double d2 = x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;
      if ( d2 < getD2CutOff() ) {
	// Overlap

	double alpha1 = alpha_fit_vector[i];
	double alpha2 = alpha_ref_vector[j];
	double a_ak = getA_ak( alpha1, alpha2 );
	double a_oe = getOverlap0( alpha1, alpha2 );
	double expMinus_ = exp( - d2 * a_ak );
	o += a_oe * expMinus_;
	
	// Gradient terms
	if (needGrad ){
	
	  double derivative0_ = -2.0 * a_oe * a_ak;
	
	  x_diff *= derivative0_ * expMinus_;
	  y_diff *= derivative0_ * expMinus_;
	  z_diff *= derivative0_ * expMinus_;
	  
	  const double  xp2 = x + x;
	  const double  yp2 = y + y;
	  const double  zp2 = z + z;
	  
	  // Gradient
	  g[ 0 ] += x_diff * xp2;
	  //            g[ 1 ] += x_diff * dxdx;
	  g[ 2 ] += x_diff * zp2;
	  g[ 3 ] -= x_diff * yp2;  // -
	  g[ 4 ] += x_diff;  // Tx
	  
	  g[ 0 ] += y_diff * yp2;
	  g[ 1 ] -= y_diff * zp2;  // -
	  //            g[ 2 ] += y_diff * dydy;
	  g[ 3 ] += y_diff * xp2;
	  g[ 5 ] += y_diff;  // Ty
	  
	  g[ 0 ] += z_diff * zp2;
	  g[ 1 ] += z_diff * yp2;
	  g[ 2 ] -= z_diff * xp2;  // -
	  //            g[ 3 ] += z_diff * dzdz;
	  g[ 6 ] += z_diff;  // Tz
	}
      }
    }
  }
}


double ComputeShapeOverlap(
				     const float* ref_coord,
				     const std::vector<double>& alpha_ref_vector,
				     const std::vector< unsigned >& ref_volume_atom_index_vector,
				     const float* fit_coord,
				     const std::vector<double>& alpha_fit_vector,
				     const std::vector< unsigned >& fit_volume_atom_index_vector
				     //////	     const GrantPickupParameters& gpparams
				     )
{
  double overlap = 0.0;

  ComputeOverlapAndAnalyticGradient(
					      overlap,
					      NULL,      
					      ////   NULL, 
					      NULL,
					      ref_coord,
					      alpha_ref_vector,
					      ref_volume_atom_index_vector,
					      fit_coord,
					      alpha_fit_vector,
					      fit_volume_atom_index_vector,
					      //////     gpparams,
					      false
					      );

  return overlap;
}


  

double ComputeFeatureOverlap( const float* ref_coord,                                                    // per ref atom, XYZ coordinates
					const std::vector<double>& alpha_ref_vector,
					const std::map< unsigned, std::vector< unsigned > >& ref_atom_type_to_index_vector_map,
					const float* fit_coord,                                                    // per fit atom, XYZ coordinates
					const std::vector<double>& alpha_fit_vector,
					const std::map< unsigned, std::vector< unsigned > >& fit_atom_type_to_index_vector_map,
					const double* m
					/////	const GrantPickupParameters& gpparams
					)                                                       // Rot/trans matrix to rotate fit to ref
{
  double overlap = 0.0;

  if( debug_print2 ) std::cerr << std::endl << "^^^^ComputeFeatureOverlap : m=" << m << std::endl;

  // Iteration over atom types
  for( std::map< unsigned, std::vector< unsigned > >::const_iterator ref_ptr = ref_atom_type_to_index_vector_map.begin();
       ref_ptr != ref_atom_type_to_index_vector_map.end();
       ++ref_ptr ){
    const unsigned atom_type = ref_ptr->first;

    if( debug_print2 ) std::cerr << "^^^^ComputeFeatureOverlap : atom_type: " << atom_type << std::endl;
    
    const std::vector< unsigned >& ref_index_vector = ref_ptr->second;
    std::map< unsigned, std::vector< unsigned > >::const_iterator fit_ptr = fit_atom_type_to_index_vector_map.find( atom_type );
    if( fit_ptr != fit_atom_type_to_index_vector_map.end() ){
      const std::vector< unsigned >& fit_index_vector = fit_ptr->second;
      
      // Loop over fit atoms
      for( unsigned i : fit_index_vector ){

	
	
	double old_x = fit_coord[ ( i * 3 )     ];
	double old_y = fit_coord[ ( i * 3 ) + 1 ];
	double old_z = fit_coord[ ( i * 3 ) + 2 ];

	double x, y, z;
	if( m == NULL ){
	  x = old_x;
	  y = old_y;
	  z = old_z;
	}else{
	  x = m[  9 ] + m[ 0 ] * old_x + m[ 1 ] * old_y + m[ 2 ] * old_z;
	  y = m[ 10 ] + m[ 3 ] * old_x + m[ 4 ] * old_y + m[ 5 ] * old_z;
	  z = m[ 11 ] + m[ 6 ] * old_x + m[ 7 ] * old_y + m[ 8 ] * old_z;
	}


	if( debug_print2 ) std::cerr << "\t" << "i: " << i << "\t" << x << "\t" << y << "\t" << z << std::endl;

	
	// Loop over fit atoms
	for( unsigned j : ref_index_vector ){

	  if( debug_print2 ) std::cerr << "\t\t" << " j: " << j << "\t" << ref_coord[ ( j * 3 ) ] << "\t" << ref_coord[ ( j * 3 + 1 ) ] << "\t" << ref_coord[ ( j * 3 + 2 )] << std::endl;
	  
	  double x_diff = x - ref_coord[ ( j * 3 )     ];
	  double y_diff = y - ref_coord[ ( j * 3 ) + 1 ];
	  double z_diff = z - ref_coord[ ( j * 3 ) + 2 ];
	  
	  double d2 = x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;

	  if( debug_print2 ) std::cerr << "\t\t\t" << "i: " << i << " j: " << j << "\t" << x_diff << "\t" << y_diff << "\t" << z_diff << std::endl;
	  if( debug_print2 ) std::cerr << "\t\t\t" << "i: " << i << " j: " << j << "\td2: " << d2 << "\tthreshold: " << getD2CutOff() << std::endl;


	  
	  // Overlap
	  if ( d2 < getD2CutOff() ) {
	    double alpha1 = alpha_fit_vector[i];
	    double alpha2 = alpha_ref_vector[j];
	    double a_ak = getA_ak( alpha1, alpha2 );
	    double a_oe = getOverlap0( alpha1, alpha2 );
	    double expMinus_ = exp( -d2 * a_ak );;
	    double add_overlap = a_oe * expMinus_;
	    overlap += add_overlap;

	    if( debug_print2 ){
	    std::cerr << "\t\t\t\t ^^^^^ i: "
		      << i << " j: " << j
		      << " a_ak: " << a_ak
		      << " a_oe: " << a_oe
		      << " expMinus_: " << expMinus_
		      << " add_overlap: " << add_overlap
		      << " overlap: " << overlap
		      << std::endl;
	    }
	    
	  }									  }
      }
    }
  }
  return  overlap;
}




double ComputeFeatureOverlap(
				       const float* ref_coord,
				       const std::vector<double>& alpha_ref_vector,
				       const std::map< unsigned, std::vector< unsigned > >& ref_atom_type_to_index_vector_map,
				       const float* fit_coord,
				       const std::vector<double>& alpha_fit_vector,
				       const std::map< unsigned, std::vector< unsigned > >& fit_atom_type_to_index_vector_map
				       )
{
  return ComputeFeatureOverlap(
					 ref_coord,
					 alpha_ref_vector,
					 ref_atom_type_to_index_vector_map,
					 fit_coord,
					 alpha_fit_vector,
					 fit_atom_type_to_index_vector_map,
					 NULL
					 );
}

    
double dCutOff = std::numeric_limits<double>::max();
double d2CutOff = std::numeric_limits<double>::max(); 

// Use Cut Off?
void setUseCutOff( bool useCutOff = false )
{
  if( useCutOff ){
    dCutOff = default_dCutOff;
    d2CutOff = dCutOff * dCutOff;
  }else{
    dCutOff = std::numeric_limits<double>::max();
    d2CutOff = std::numeric_limits<double>::max(); 
  }
}


  
double getD2CutOff(){
  return d2CutOff;
}


double getAlpha( double r ){

  double alpha( static_cast<double>( 0 ) );
  if( r > static_cast<double>( 0 ) ){

    alpha = kappa/( r*r );

    if( debug_print2 ){
      std::cerr << "\t\t\t\t\tkappa: " << kappa
		<< "<< r:" << r
		<< "alpha: " << alpha
		<< std::endl;
    }
  }
  return alpha;
}

double getA_ak( double alpha1, double alpha2 ){
  double a_ak = static_cast< double >( 0 );
  double sum = alpha1 + alpha2;
  if( sum > static_cast<double>( 0 ) ){
    a_ak = alpha1*alpha2 / sum;
  }

  if( debug_print2 ){
    std::cerr << "\t\t\t\t\t "
	      << "^^^^^^ getA_ak( unsigned i, unsigned j ):"
	      << " alpha1: " << alpha1
	      << " alpha2: " << alpha2
	      << std::endl;
  }
  
  return a_ak;
}



double getOverlap0( double alpha1, double alpha2 ){
  double a_pi = static_cast< double>( 0 );
  double sum = alpha1 + alpha2;
  if( sum > static_cast<double>( 0 ) ){
    a_pi = pi / sum;
  }
  double a_oe = _p2_ * a_pi * std::sqrt( a_pi );
  return a_oe;
}


void setAlpha( const double* rad,
			 const unsigned natom,
			 std::vector< double >& alpha_vec ){

  if( debug_print2 ) std::cerr << "--- natom: " << natom << std::endl;
  
  alpha_vec.clear();
  for( unsigned k=0; k < natom; ++k ){
    double r = rad[ k ];
    double alpha = getAlpha( r );
    alpha_vec.push_back( alpha );

    if( debug_print2 ) std::cerr << "---" << r << "\t" << alpha << "\t" << alpha_vec[k] << std::endl;
  }
}

} // namespace Align3D

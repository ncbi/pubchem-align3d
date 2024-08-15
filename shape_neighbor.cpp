/*  $Id: shape_neighbor.cpp 685605 2024-07-26 12:29:33Z thiessen $
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
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <limits>


#include "shape_functions.hpp"
#include "shape_constants.hpp"
#include "shape_debug.hpp"


namespace Align3D {

void Neighbor_Conformers(
				      const float* ref_coord,
				      const std::vector<double>& alpha_ref_vector,
				      const std::vector< unsigned >& ref_volumeAtomIndexVector,
				      const std::map< unsigned, std::vector< unsigned > >& ref_colorAtomType2IndexVectorMap,
				      const double ref_sov,
				      const double ref_sof,
				      const float* fit_coord,
				      const std::vector<double>& alpha_fit_vector,
				      const std::vector< unsigned >& fit_volumeAtomIndexVector,
				      const std::map< unsigned, std::vector< unsigned > >& fit_colorAtomType2IndexVectorMap,
				      const double fit_sov,
				      const double fit_sof,
				      const bool considerColorAtoms,
				      const unsigned max_preiters,
				      const unsigned max_postiters,
				      const double opt_param,
				      float* matrix,
				      double& nbr_st,
				      double& nbr_ct
				      ){
  
  if (debug_print0) {
    std::cerr << "debug_print1: " << debug_print1 << std::endl;
    std::cerr << "debug_print2: " << debug_print2 << std::endl;
  }
  
  // Sums needed for further use
  const double  st_sum_sov = ( ref_sov + fit_sov );  
  const double ct_sum_sov = ( ref_sof + fit_sof );

  double st_ratio = static_cast<double>(0);
  {
    double minSov = std::min( ref_sov, fit_sov );
    if( st_sum_sov > minSov ){
      st_ratio = minSov/( st_sum_sov - minSov );
    }
  }

  double ct_ratio = static_cast<double>(0);
  {
    double minSof = std::min( ref_sof, fit_sof );
    if( ct_sum_sov > minSof ){
      ct_ratio = minSof/( ct_sum_sov - minSof );
    }
  }
       
  // Qrat stuff
  const float qrat_threshold = 0.7225; // 0.85*0.85;
    
  unsigned ref_qrat = 1000;
  {
    std::vector< float > ref_coord_work_vec(alpha_ref_vector.size() * 3);
    memcpy(&(ref_coord_work_vec[0]), ref_coord, (alpha_ref_vector.size() * 3 * sizeof(float)));
    TransformCoordToStericFrameAndCalculateQrat( alpha_ref_vector.size(), qrat_threshold, &( ref_coord_work_vec[0] ), ref_qrat );
  }
  unsigned fit_qrat = 1000;
  {
    std::vector< float > fit_coord_work_vec(alpha_fit_vector.size() * 3);
    memcpy(&(fit_coord_work_vec[0]), fit_coord, (alpha_fit_vector.size() * 3 * sizeof(float)));
    TransformCoordToStericFrameAndCalculateQrat( alpha_fit_vector.size(), qrat_threshold, &( fit_coord_work_vec[0] ), fit_qrat );
  }

  bool is_qrat_nonzero = (ref_qrat > 0 || fit_qrat > 0);
  
  unsigned  npose = 0;
  if( is_qrat_nonzero ){
    npose = 16;
  }else{
    npose = 4;
  }

  if( debug_print1 ) std::cerr << "*** st_sum_sov: " << st_sum_sov << std::endl;
  if( debug_print1 ) std::cerr << "*** ct_sum_sov: " << ct_sum_sov << std::endl;
  if( debug_print1 ) std::cerr << "*** st_ratio: " << st_ratio << std::endl;
  if( debug_print1 ) std::cerr << "*** ct_ratio: " << ct_ratio << std::endl;

  // Globally best
  
  nbr_st = static_cast<double>(0);
  nbr_ct = static_cast<double>(0);
  double nbr_comb_prelim = static_cast<double>(0);
  double nbr_comb = static_cast<double>(0);
  double  best_qom[ 7 ];  // Best quaternion/trans matrix (optimized neighbor pose)
    
  double  best_qom_pose[ 16 ][ 7 ];
  double  best_mom_pose[ 16 ][ 12 ];
  
  if( debug_print1 ) std::cerr << "^^^*** max_preiters: " << max_preiters << std::endl;
  if( debug_print1 ) std::cerr << "^^^*** max_postiters: " << max_postiters << std::endl;
  if( debug_print1 ) std::cerr << "^^^*** st_sum_sov: " << st_sum_sov << std::endl;
  if( debug_print1 ) std::cerr << "^^^*** ct_sum_sov: " << ct_sum_sov << std::endl;
  if( debug_print1 ) std::cerr << "^^^*** considerColorAtoms: " << considerColorAtoms << std::endl;
  if( debug_print1 ) std::cerr << "*** npose: " << npose << std::endl;
  if( debug_print1 ) std::cerr << std::endl;

  // Preparing initial poses
  
  for ( unsigned  ibest_pose = 0;  ibest_pose < npose;  ibest_pose++ ) {
    
    if( debug_print1 ) std::cerr << "Processing initial pose: " << ibest_pose << "\t is_qrat_nonzero: " << is_qrat_nonzero << "\t  npose: " << npose << std::endl;
      
    // Point to appropriate transformation matrix for steric pose
      
    double  qom[ 7 ];  // Grid-based pose 4x1 quaternion vector and 3x1 translation vector
    const double  *mom;  // Grid-based pose 3x3 rotation matrix and 3x1 translation vector
    double  st_tot = (double ) 0;  // AB volume overlap
    double  g[ 7 ];  // Analytic gradient for 4x1 quaternion vector and 3x1 translation vector (really 6x1 since [0] element is zero)

    if( debug_print1 ) std::cerr << "\tBefore qrat-based selection of poses, is_qrat_nonzero: " << is_qrat_nonzero << " st_tot: " << st_tot << std::endl;
    
    if ( is_qrat_nonzero ) {
      
      if( debug_print1 ) std::cerr << "\t\tis_qrat_nonzero: " << is_qrat_nonzero << " -- Select mom for ibest_pose, where ibest_pose=" << ibest_pose << std::endl;

      memcpy( best_qom_pose[ ibest_pose ], qmatrix_typed[ ibest_pose ], 7 * sizeof( double ) );    
      memcpy( best_mom_pose[ ibest_pose ], matrix_typed[ ibest_pose ], 12 * sizeof( double ) );
    } else {
      if( debug_print1 ) std::cerr << "\t\tis_qrat_nonzero: " << is_qrat_nonzero << " -- Select 8 mom's for ibest_pose, where ibest_pose=" << ibest_pose << std::endl;
      
      // Find best m_loc_ from mloc_ = 0, 1, ..., 6, set it to m_loc, and apply insyed of the best pose from the grid
      unsigned m_loc = 0;
      double F = (double ) -1;
      for ( unsigned  m_loc_ = 0; m_loc_ < 7; ++m_loc_ ){
	double  mc[ 12 ];  // Scratch array for rot/trans matrix	  
	memcpy( qom, pqmatrix_typed[ ibest_pose ][ m_loc_ ], 7 * sizeof( double ) );
	mom = pmatrix_typed[ ibest_pose ][ m_loc_ ];
	
	ComputeOverlapAndAnalyticGradient(
					 st_tot,   
					 g,   
					 mom,   
					 ref_coord,
					 alpha_ref_vector,
					 ref_volumeAtomIndexVector,
					 fit_coord,
					 alpha_fit_vector,
					 fit_volumeAtomIndexVector
					  );

	// Compute shape Tanimoto
	double fst = 0.0; 
	if( st_sum_sov > 0 && st_sum_sov > st_tot ){
	  fst = st_tot / ( st_sum_sov - st_tot );
	}

	// Compute feature Tanimoto
	double ct_tot = 0.0;
	double fct = 0.0;
	if( considerColorAtoms ){
	  if( ct_sum_sov > 0 ){
	    VQuaternionTransArray2RotTransMatrix( mc, qom, false );  // Already normalized
	    ct_tot = ComputeFeatureOverlap(
					   ref_coord,
					   alpha_ref_vector,
					   ref_colorAtomType2IndexVectorMap,
					   fit_coord,
					   alpha_fit_vector,
					   fit_colorAtomType2IndexVectorMap,
					   mc
					   );


	    
	    if( ct_sum_sov > ct_tot ){
	      fct = ct_tot / ( ct_sum_sov - ct_tot );
	    }
	  }
	}
	  	
	double F_ = (double) 0;
	if( considerColorAtoms ){
	  F_ = opt_param*(fst-fct) + fct;
	}else{
	  F_ = fst;
	}

	if( debug_print1 ) std:: cerr << "\t\tibest_pose: " << ibest_pose << "\tm_loc_: " << m_loc_ << "\tst_tot: " << st_tot << "\tst_sum_sov: " << st_sum_sov << "\tct_tot: " << ct_tot << "\tct_sum_sov: " << ct_sum_sov << "\tfst: " << fst << "\t fct: " << fct << "\t F_: " << F_ << std::endl;
	
	// Find best m_loc
	if( F_ > F ){
	  m_loc = m_loc_;
	  F = F_;
	}
      }

      if( debug_print1 ) std:: cerr << "ibest_pose: " << ibest_pose << " FOUND m_loc: " << m_loc << "\tF: " <<  F << "\n";

	
      // Setting "best" m_loc for the current best_pose ibest_pose
      memcpy( best_qom_pose[ ibest_pose ], pqmatrix_typed[ ibest_pose ][ m_loc ], 7 * sizeof( double ) );
      memcpy( best_mom_pose[ ibest_pose ], pmatrix_typed[ ibest_pose ][ m_loc ], 12 * sizeof( double ) );
    }  // end if "qrat"
  } // end loop ibest_pose


  const double  max_quaternion_step  = 0.075;   // Maximum step size for quaternion
  const double  max_translation_step = 0.500;   // Maximum step size for translation
  const double  min_quaternion_step  = 0.0002;  // Convergence criteria for quaternion
  const double  min_translation_step = 0.0020;  // Convergence criteria for translation
  const double  min_volume_overlap   = 0.0020;  // Convergence criteria for ST AB volume overlap
	
  
  // Previous step quaternion/translation and respective 1st derivative gradient
  double  old_quattrans_pose[ 16 ][ 7 ];
  double old_gradient_pose[ 16 ][ 7 ];
  double qstep_size_pose[16];
  double tstep_size_pose[ 16 ];
  bool allowed_pose[16];
  double st_tot_pose[ 16 ];
  double fcomb_best_pose[ 16 ];
  
  for ( unsigned  ibest_pose = 0;  ibest_pose < npose;  ibest_pose++ ) {
    allowed_pose[ ibest_pose ] = true;
    fcomb_best_pose[ ibest_pose ] = static_cast<double>( 0 );
    qstep_size_pose[ ibest_pose ]  = -0.0010;  // Quaternion step size using gradient for direction
    tstep_size_pose[ ibest_pose ]  = -0.0100;  // Translation step size using gradient for direction
    
    // Evaluate the gradient at the starting point

    ComputeOverlapAndAnalyticGradient(
				      st_tot_pose[ ibest_pose ],
				      old_gradient_pose[ ibest_pose ],
				      best_mom_pose[ ibest_pose ],
				      ref_coord,
				      alpha_ref_vector,
				      ref_volumeAtomIndexVector,
				      fit_coord,
				      alpha_fit_vector,
				      fit_volumeAtomIndexVector
				      );
   
    if( debug_print1 ){
      float fst = st_tot_pose[ ibest_pose ] / ( st_sum_sov - st_tot_pose[ ibest_pose ]  );
      std::cerr << "\tibest_pose: " << ibest_pose << " Beginning: st_tot: " << st_tot_pose[ ibest_pose ] << "; fst: " << fst << std::endl;
    }
    
  }
  
  // Loops over cycles, allowed poses and iterations
  for( unsigned int icycle=0; icycle < 2; ++icycle ){
    for ( unsigned  ibest_pose = 0;  ibest_pose < npose;  ibest_pose++ ) {
      
      // Beginning of the cycle -- start with old gradient
      double  g[ 7 ];
      memcpy( g, old_gradient_pose[ ibest_pose ], 7 * sizeof( double ) );
      double& qstep_size = qstep_size_pose[ ibest_pose ];
      double& tstep_size = tstep_size_pose[ ibest_pose ];
      double *old_gradient = old_gradient_pose[ ibest_pose ];
      double* old_quattrans = old_quattrans_pose[ ibest_pose ];
      double& st_tot = st_tot_pose[ ibest_pose ];
      double  *qom = best_qom_pose[ ibest_pose ];
      //	double  *mom = best_mom_pose[ ibest_pose ];

	  
      if( icycle == 1 ){ // Beginning of the cycle with icycle=0

	
	if( debug_print1 ) std::cerr << "*** icycle: " << icycle
				    << " ibest_pose: " << ibest_pose
				    << " nbr_comb_prelim: " << nbr_comb_prelim
				    << " fcomb_best_pose[ ibest_pose ]: " << fcomb_best_pose[ ibest_pose ]
				    << " fcomb_best_pose[ ibest_pose ]/nbr_comb_prelim: " << fcomb_best_pose[ ibest_pose ]/nbr_comb_prelim
				    << std::endl;
	
	if( fcomb_best_pose[ ibest_pose ] <= 0.7 * nbr_comb_prelim ){
	  allowed_pose[ ibest_pose ] = false;
	}

	if( debug_print1 ){
	  if( allowed_pose[ ibest_pose ] ){
	    std::cerr << "*** pose: " << ibest_pose << " allowed" << std::endl;
	  }else{
	    std::cerr << "*** pose: " << ibest_pose << " not allowed" << std::endl;
	  }
	}
      }

      if( allowed_pose[ ibest_pose ] ){
	if( debug_print1 ) std::cerr << "\tibest_pose: " << ibest_pose << " Compute (non-grid-based) shape volume overlap for the current pose ibest_pose " << ibest_pose << " and its \"best\" starting qom" << std::endl;
 

	
	// Iterations
	unsigned start_iter = 0;
	unsigned end_iter = 0;
	if( icycle == 0 ){
	  end_iter = max_preiters;
	}else{
	  start_iter = max_preiters;
	    end_iter = max_preiters + max_postiters;
	}
	
	for ( unsigned iter = start_iter;  iter < end_iter;  iter++ ) {
	  // Determine the step to take
	  double  quat_squared;  // Scratch for deriving fourth quaternion value [0] from other three [1-3]
	  double  step[ 7 ];  // Actual step taken
	  //                    unsigned  ustep[ 6 ];  // Standard step size

	  // If not first iteration and if we bracketed the maximum, use old step
	  bool  converged = false;
	  if ( 0 == iter ) {
	    // 1st iteration through, use the default step size
	    //                        double  
	    step[ 1 ] = qstep_size * g[ 1 ];
	    step[ 2 ] = qstep_size * g[ 2 ];
	    step[ 3 ] = qstep_size * g[ 3 ];
	    step[ 4 ] = tstep_size * g[ 4 ];
	    step[ 5 ] = tstep_size * g[ 5 ];
	    step[ 6 ] = tstep_size * g[ 6 ];
	  } else {
	    if ( std::signbit( g[ 1 ] ) != std::signbit( old_gradient[ 1 ] ) ) {
	      step[ 1 ] = ( ( ( qom[ 1 ] * fabs( old_gradient[ 1 ] ) ) + ( old_quattrans[ 1 ] * fabs( g[ 1 ] ) ) )
			    / ( fabs( old_gradient[ 1 ] ) + fabs( g[ 1 ] ) ) ) - qom[ 1 ];
	      double  new_step = qstep_size * g[ 1 ];
	      if ( fabs( step[ 1 ] ) > fabs( new_step ) ) {
		step[ 1 ] *= fabs( new_step / step[ 1 ] );
	      }
	    } else {
	      step[ 1 ] = qstep_size * g[ 1 ];
	    }

	    if ( std::signbit( g[ 2 ] ) != std::signbit( old_gradient[ 2 ] ) ) {
	      step[ 2 ] = ( ( ( qom[ 2 ] * fabs( old_gradient[ 2 ] ) ) + ( old_quattrans[ 2 ] * fabs( g[ 2 ] ) ) )
			    / ( fabs( old_gradient[ 2 ] ) + fabs( g[ 2 ] ) ) ) - qom[ 2 ];
	      double  new_step = qstep_size * g[ 2 ];
	      if ( fabs( step[ 2 ] ) > fabs( new_step ) ) {
		step[ 2 ] *= fabs( new_step / step[ 2 ] );
	      }
	    } else {
	      step[ 2 ] = qstep_size * g[ 2 ];
	    }

	    if ( std::signbit( g[ 3 ] ) != std::signbit( old_gradient[ 3 ] ) ) {
	      step[ 3 ] = ( ( ( qom[ 3 ] * fabs( old_gradient[ 3 ] ) ) + ( old_quattrans[ 3 ] * fabs( g[ 3 ] ) ) )
			    / ( fabs( old_gradient[ 3 ] ) + fabs( g[ 3 ] ) ) ) - qom[ 3 ];
	      double  new_step = qstep_size * g[ 3 ];
	      if ( fabs( step[ 3 ] ) > fabs( new_step ) ) {
		step[ 3 ] *= fabs( new_step / step[ 3 ] );
	      }
	    } else {
	      step[ 3 ] = qstep_size * g[ 3 ];
	    }

	    if ( std::signbit( g[ 4 ] ) != std::signbit( old_gradient[ 4 ] ) ) {
	      step[ 4 ] = ( ( ( qom[ 4 ] * fabs( old_gradient[ 4 ] ) ) + ( old_quattrans[ 4 ] * fabs( g[ 4 ] ) ) )
			    / ( fabs( old_gradient[ 4 ] ) + fabs( g[ 4 ] ) ) ) - qom[ 4 ];
	      double  new_step = tstep_size * g[ 4 ];
	      if ( fabs( step[ 4 ] ) > fabs( new_step ) ) {
		step[ 4 ] *= fabs( new_step / step[ 4 ] );
	      }
	    } else {
	      step[ 4 ] = tstep_size * g[ 4 ];
	    }

	    if ( std::signbit( g[ 5 ] ) != std::signbit( old_gradient[ 5 ] ) ) {
	      step[ 5 ] = ( ( ( qom[ 5 ] * fabs( old_gradient[ 5 ] ) ) + ( old_quattrans[ 5 ] * fabs( g[ 5 ] ) ) )
			    / ( fabs( old_gradient[ 5 ] ) + fabs( g[ 5 ] ) ) ) - qom[ 5 ];
	      double  new_step = tstep_size * g[ 5 ];
	      if ( fabs( step[ 5 ] ) > fabs( new_step ) ) {
		step[ 5 ] *= fabs( new_step / step[ 5 ] );
	      }
	    } else {
	      step[ 5 ] = tstep_size * g[ 5 ];
	    }

	    if ( std::signbit( g[ 6 ] ) != std::signbit( old_gradient[ 6 ] ) ) {
	      step[ 6 ] = ( ( ( qom[ 6 ] * fabs( old_gradient[ 6 ] ) ) + ( old_quattrans[ 6 ] * fabs( g[ 6 ] ) ) )
			    / ( fabs( old_gradient[ 6 ] ) + fabs( g[ 6 ] ) ) ) - qom[ 6 ];
	      double  new_step = tstep_size * g[ 6 ];
	      if ( fabs( step[ 6 ] ) > fabs( new_step ) ) {
		step[ 6 ] *= fabs( new_step / step[ 6 ] );
	      }
	    } else {
	      step[ 6 ] = tstep_size * g[ 6 ];
	    }
	  }

	  // Retain last step information in case we are not happy with "line search"
	  double  old_overlap = st_tot;
	  memcpy( old_quattrans, qom, 7 * sizeof( double ) );
	  memcpy( old_gradient, g, 7 * sizeof( double ) );

	  if( debug_print1 ) std::cerr << "\t\tibest_pose: " << ibest_pose << " iter: " << iter << "\tOld overlap: " << old_overlap << std::endl;

	  //
	  // Line search (sort of) loop
	  //
	  for ( unsigned  line_iter = 0;  false == converged;  line_iter++ ) {
	    if( debug_print1 ) std::cerr << "\t\t\tibest_pose: " << ibest_pose << " iter: " << iter << "\tline_iter: " << line_iter << " -- started" << std::endl;
	
	    // Check that the absolute max step size does not go beyond some reasonable size
	    double  mtstep = fmax( fmax( fabs( step[ 4 ] ), fabs( step[ 5 ] ) ), fabs( step[ 6 ] ) );
	    if ( mtstep > max_translation_step ) {
	  
	      double  tscale_factor = max_translation_step / mtstep;

	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose: " << ibest_pose << " iter: " << iter << "\tline_iter: " << line_iter << " too large max step; applying scaling factor tscale_factor: " << tscale_factor << std::endl;
       	  
	      if ( fabs( step[ 4 ] ) > max_translation_step ) step[ 4 ] *= tscale_factor;
	      if ( fabs( step[ 5 ] ) > max_translation_step ) step[ 5 ] *= tscale_factor;
	      if ( fabs( step[ 6 ] ) > max_translation_step ) step[ 6 ] *= tscale_factor;
	    }

	    double  mqstep = fmax( fmax( fabs( step[ 1 ] ), fabs( step[ 2 ] ) ), fabs( step[ 3 ] ) );
	    if ( mqstep > max_quaternion_step ) {
	      double  qscale_factor = max_quaternion_step / mqstep;

	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose:" << ibest_pose << " iter: " << iter << "\tline_iter: " << line_iter << " too large max step; applying scaling factor qscale_factor: " << qscale_factor << std::endl;
	  
	      if ( fabs( step[ 1 ] ) > max_quaternion_step ) step[ 1 ] *= qscale_factor;
	      if ( fabs( step[ 2 ] ) > max_quaternion_step ) step[ 2 ] *= qscale_factor;
	      if ( fabs( step[ 3 ] ) > max_quaternion_step ) step[ 3 ] *= qscale_factor;
	    }


	
	    if ( min_quaternion_step > mqstep && min_translation_step > mtstep ) {  // Did we converge (so to speak)?
	      st_tot = 0.0;  // We didn't compute the new overlap yet, ensure we use the old one.
	      converged = true;

	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose: "
					  << ibest_pose
					  << " iter: "
					  << iter
					  << " line_iter: "
					  << line_iter
					  << "-- line search converged because"
					  << " " << min_quaternion_step
					  << " " << mqstep
					  << " " << min_translation_step
					  << " " << mtstep
					  << ". BREAK " << std::endl;
	      break;
	    }else{
	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose: " << ibest_pose << " iter: " << iter << " line_iter: " << line_iter << "-- line search is not converged yet because " << min_quaternion_step << " " << mqstep << " " << min_translation_step << " " << mtstep << std::endl;
	    }
	
	
	    // Derive the 4th component of the quaternion
	    quat_squared = step[ 1 ] * step[ 1 ] + step[ 2 ] * step[ 2 ] + step[ 3 ] * step[ 3 ];
	    step[ 0 ] = sqrt( 1.0 - quat_squared );

	    // Update quaternion with step size (getting new rot/trans matrix to use in the process)

	    double  mc[ 12 ], qc[ 7 ], tm1[ 12 ], tm2[ 12 ];  // Scratch array for rot/trans matrix

	    // Functions originally from PubChem::Data
	    VQuaternionTransArray2RotTransMatrix( tm1, qom, false );  // Already normalized
	    VQuaternionTransArray2RotTransMatrix( tm2, step, false );  // Already normalized
	    VCombineRotTransMatrixMatrix( mc, tm2, tm1 );  // New rot/trans matrix
	    VRotTransMatrix2QuaternionTransArray( qc, mc );  // New quaternion

	    // QC is not used in below function since always assumed to be (1,0,0,0) ... need to remove it...
	    
	    ComputeOverlapAndAnalyticGradient(
					      st_tot,  
					      g,
					      mc,   
					      ref_coord,
					      alpha_ref_vector,
					      ref_volumeAtomIndexVector,
					      fit_coord,
					      alpha_fit_vector,
					      fit_volumeAtomIndexVector
					      );
				       
	    // Did we make a good step?

	    if ( st_tot > old_overlap ) {
	      memcpy( qom, qc, 7 * sizeof( double ) );


	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose: " << ibest_pose << " iter: " << iter << " line_iter: " << line_iter << "-- we made a good step step because st_tot increased. old_overlap: " << old_overlap << " st_tot: " << st_tot << ". BREAK" << std::endl;
	    
	      break;  // Line search is complete
	    }else{
	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose: " << ibest_pose << " iter: " << iter << " line_iter: " << line_iter << "-- we DID NOT make a good step step because st_tot did not increase. old_overlap: " << old_overlap << " st_tot: " << st_tot << std::endl;

	    }
       

	    // Overlap at new point is worse, cut back

	    if ( 2 == line_iter ) {

	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose: "<< ibest_pose << " iter: " << iter << " line_iter: " << line_iter << "-- we line_iter=2 and we did not make a good step step because st_tot did not increase. Decrease steps." << std::endl;

	      qstep_size *= 0.1;
	      tstep_size *= 0.1;

	      step[ 1 ] = qstep_size * old_gradient[ 1 ];
	      step[ 2 ] = qstep_size * old_gradient[ 2 ];
	      step[ 3 ] = qstep_size * old_gradient[ 3 ];
	      step[ 4 ] = tstep_size * old_gradient[ 4 ];
	      step[ 5 ] = tstep_size * old_gradient[ 5 ];
	      step[ 6 ] = tstep_size * old_gradient[ 6 ];

	      qstep_size *= 5.0;
	      tstep_size *= 5.0;
	  
	    } else if ( 3 < line_iter ) {
	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose: " << ibest_pose << " iter: " << iter << " line_iter: " << line_iter << "-- we line_iter>3 and we did not make a good step step because st_tot did not increase. Set converged. BREAK." << std::endl;

	  
	      converged = true;
	      memcpy( qom, qc, 7 * sizeof( double ) );
	      break;
	    
	    } else {
	      if( debug_print1 ) std::cerr << "\t\t\tibest_pose: " << ibest_pose << " iter: " << iter << " line_iter: " << line_iter << "-- we line_iter = 1,3 and we did not make a good step step because st_tot did not increase. Changed steps." << std::endl;


	  
	      if ( std::signbit( g[ 1 ] ) != std::signbit( old_gradient[ 1 ] ) ) {
		step[ 1 ] = ( ( ( qc[ 1 ] * fabs( old_gradient[ 1 ] ) ) + ( old_quattrans[ 1 ] * fabs( g[ 1 ] ) ) )
			      / ( fabs( old_gradient[ 1 ] ) + fabs( g[ 1 ] ) ) ) - qc[ 1 ];
		double  new_step = qstep_size * g[ 1 ];
		if ( fabs( step[ 1 ] ) > fabs( new_step ) ) {
		  step[ 1 ] *= fabs( new_step / step[ 1 ] );
		}
	      } else if ( 1.0 > fabs( g[ 1 ] ) ) {
		step[ 1 ] = qstep_size * g[ 1 ];
	      } else if ( fabs( g[ 1 ] ) > fabs( old_gradient[ 1 ] ) ) {  // Going wrong way (relative to other components)?
		step[ 1 ] += qstep_size * g[ 1 ];
	      } else {
		double  delta = ( g[ 1 ] * ( step[ 1 ] / ( old_gradient[ 1 ] - g[ 1 ] ) ) );
		if ( fabs( delta ) > fabs( step[ 1 ] * 0.1 ) && fabs( delta ) > 0.001 ) {
		  delta *= 0.0005 / fabs( delta );
		}
		step[ 1 ] += delta;
	      }

	      if ( std::signbit( g[ 2 ] ) != std::signbit( old_gradient[ 2 ] ) ) {
		step[ 2 ] = ( ( ( qc[ 2 ] * fabs( old_gradient[ 2 ] ) ) + ( old_quattrans[ 2 ] * fabs( g[ 2 ] ) ) )
			      / ( fabs( old_gradient[ 2 ] ) + fabs( g[ 2 ] ) ) ) - qc[ 2 ];
		double  new_step = qstep_size * g[ 2 ];
		if ( fabs( step[ 2 ] ) > fabs( new_step ) ) {
		  step[ 2 ] *= fabs( new_step / step[ 2 ] );
		}
	      } else if ( 1.0 > fabs( g[ 2 ] ) ) {
		step[ 2 ] = qstep_size * g[ 2 ];
	      } else if ( fabs( g[ 2 ] ) > fabs( old_gradient[ 2 ] ) ) {  // Going wrong way (relative to other components)?
		step[ 2 ] += qstep_size * g[ 2 ];
	      } else {
		double  delta = ( g[ 2 ] * ( step[ 2 ] / ( old_gradient[ 2 ] - g[ 2 ] ) ) );
		if ( fabs( delta ) > fabs( step[ 2 ] * 0.1 ) && fabs( delta ) > 0.001 ) {
		  delta *= 0.0005 / fabs( delta );
		}
		step[ 2 ] += delta;
	      }

	      if ( std::signbit( g[ 3 ] ) != std::signbit( old_gradient[ 3 ] ) ) {
		step[ 3 ] = ( ( ( qc[ 3 ] * fabs( old_gradient[ 3 ] ) ) + ( old_quattrans[ 3 ] * fabs( g[ 3 ] ) ) )
			      / ( fabs( old_gradient[ 3 ] ) + fabs( g[ 3 ] ) ) ) - qc[ 3 ];
		double  new_step = qstep_size * g[ 3 ];
		if ( fabs( step[ 3 ] ) > fabs( new_step ) ) {
		  step[ 3 ] *= fabs( new_step / step[ 3 ] );
		}
	      } else if ( 1.0 > fabs( g[ 3 ] ) ) {
		step[ 3 ] = qstep_size * g[ 3 ];
	      } else if ( fabs( g[ 3 ] ) > fabs( old_gradient[ 3 ] ) ) {  // Going wrong way (relative to other components)?
		step[ 3 ] += qstep_size * g[ 3 ];
	      } else {
		double  delta = ( g[ 3 ] * ( step[ 3 ] / ( old_gradient[ 3 ] - g[ 3 ] ) ) );
		if ( fabs( delta ) > fabs( step[ 3 ] * 0.1 ) && fabs( delta ) > 0.001 ) {
		  delta *= 0.0005 / fabs( delta );
		}
		step[ 3 ] += delta;
	      }

	      if ( std::signbit( g[ 4 ] ) != std::signbit( old_gradient[ 4 ] ) ) {
		step[ 4 ] = ( ( ( qc[ 4 ] * fabs( old_gradient[ 4 ] ) ) + ( old_quattrans[ 4 ] * fabs( g[ 4 ] ) ) )
			      / ( fabs( old_gradient[ 4 ] ) + fabs( g[ 4 ] ) ) ) - qc[ 4 ];
		double  new_step = qstep_size * g[ 4 ];
		if ( fabs( step[ 4 ] ) > fabs( new_step ) ) {
		  step[ 4 ] *= fabs( new_step / step[ 4 ] );
		}
	      } else if ( 1.0 > fabs( g[ 4 ] ) ) {
		step[ 4 ] = qstep_size * g[ 4 ];
	      } else if ( fabs( g[ 4 ] ) > fabs( old_gradient[ 4 ] ) ) {  // Going wrong way (relative to other components)?
		step[ 4 ] += qstep_size * g[ 4 ];
	      } else {
		double  delta = ( g[ 4 ] * ( step[ 4 ] / ( old_gradient[ 4 ] - g[ 4 ] ) ) );
		if ( fabs( delta ) > fabs( step[ 4 ] * 0.5 ) && fabs( delta ) > 0.01 ) {
		  delta *= 0.005 / fabs( delta );
		}
		step[ 4 ] += delta;
	      }

	      if ( std::signbit( g[ 5 ] ) != std::signbit( old_gradient[ 5 ] ) ) {
		step[ 5 ] = ( ( ( qc[ 5 ] * fabs( old_gradient[ 5 ] ) ) + ( old_quattrans[ 5 ] * fabs( g[ 5 ] ) ) )
			      / ( fabs( old_gradient[ 5 ] ) + fabs( g[ 5 ] ) ) ) - qc[ 5 ];
		double  new_step = qstep_size * g[ 5 ];
		if ( fabs( step[ 5 ] ) > fabs( new_step ) ) {
		  step[ 5 ] *= fabs( new_step / step[ 5 ] );
		}
	      } else if ( 1.0 > fabs( g[ 5 ] ) ) {
		step[ 5 ] = qstep_size * g[ 5 ];
	      } else if ( fabs( g[ 5 ] ) > fabs( old_gradient[ 5 ] ) ) {  // Going wrong way (relative to other components)?
		step[ 5 ] += qstep_size * g[ 5 ];
	      } else {
		double  delta = ( g[ 5 ] * ( step[ 5 ] / ( old_gradient[ 5 ] - g[ 5 ] ) ) );
		if ( fabs( delta ) > fabs( step[ 5 ] * 0.5 ) && fabs( delta ) > 0.01 ) {
		  delta *= 0.005 / fabs( delta );
		}
		step[ 5 ] += delta;
	      }

	      if ( std::signbit( g[ 6 ] ) != std::signbit( old_gradient[ 6 ] ) ) {
		step[ 6 ] = ( ( ( qc[ 6 ] * fabs( old_gradient[ 6 ] ) ) + ( old_quattrans[ 6 ] * fabs( g[ 6 ] ) ) )
			      / ( fabs( old_gradient[ 6 ] ) + fabs( g[ 6 ] ) ) ) - qc[ 6 ];
		double  new_step = qstep_size * g[ 6 ];
		if ( fabs( step[ 6 ] ) > fabs( new_step ) ) {
		  step[ 6 ] *= fabs( new_step / step[ 6 ] );
		}
	      } else if ( 1.0 > fabs( g[ 6 ] ) ) {
		step[ 6 ] = qstep_size * g[ 6 ];
	      } else if ( fabs( g[ 6 ] ) > fabs( old_gradient[ 6 ] ) ) {  // Going wrong way (relative to other components)?
		step[ 6 ] += qstep_size * g[ 6 ];
	      } else {
		double  delta = ( g[ 6 ] * ( step[ 6 ] / ( old_gradient[ 6 ] - g[ 6 ] ) ) );
		if ( fabs( delta ) > fabs( step[ 6 ] * 0.5 ) && fabs( delta ) > 0.01 ) {
		  delta *= 0.005 / fabs( delta );
		}
		step[ 6 ] += delta;
	      }
	    }
	    memcpy( qom, old_quattrans, 7 * sizeof( double ) );
	  }  // End of line search

	  ///// if( //// DebugOptions::debug_print1 ) std::cerr << "***** End of line search. converged: " << converged << " min_volume_overlap > ( st_tot - old_overlap ) : " << ( min_volume_overlap > ( st_tot - old_overlap ) ) << std::endl;

	  // Did we converge?
	  if ( true == converged || min_volume_overlap > ( st_tot - old_overlap ) ) {
	    if( debug_print1 ) std::cerr << "\t\tibest_pose_pose: " << ibest_pose << " iter: " << iter << " converged. converged: " << converged << " min_volume_overlap: " << min_volume_overlap << " st_tot - old_overlap: " << ( st_tot - old_overlap ) << " WILL BREAK" << std::endl;
	
	    if ( old_overlap > st_tot ) {  // Previous successful step was better... retain better step
	      st_tot = old_overlap;
	      memcpy( qom, old_quattrans, 7 * sizeof( float ) );

	      if( debug_print1 ) std::cerr << "\t\titer: " << iter << " keep previous value of st_tot. "; 
	    }else{
	      if( debug_print1 ) std::cerr << "\t\titer: " << iter << " keep previous value of st_tot. "; 
	    }

	
	    if( debug_print1 ) std::cerr << " st_tot: " << st_tot << std::endl;

	    break;  // Converged!
	  }      
	}  // End of optimization iteration

	// Compute shape Tanimoto
	double  fst = 0.0;
	if( st_sum_sov > 0 && st_sum_sov > st_tot ){
	  fst = st_tot / ( st_sum_sov - st_tot );
	}

	// Compute feature Tanimoto
	double ct_tot = 0.0;
	double fct = 0.0;
	if ( considerColorAtoms ) {
	  if( ct_sum_sov > 0 ){
	    double  mc[ 12 ];  // Scratch array for rot/trans matrix
	    VQuaternionTransArray2RotTransMatrix( mc, qom, false );  // Already normalized
	    ct_tot = ComputeFeatureOverlap(
					   ref_coord,
					   alpha_ref_vector,
					   ref_colorAtomType2IndexVectorMap,
					   fit_coord,
					   alpha_fit_vector,
					   fit_colorAtomType2IndexVectorMap,
					   mc
					   );


	    
	    if( ct_sum_sov > ct_tot ){
	      fct = ct_tot / ( ct_sum_sov - ct_tot );
	    }
	  }
	}
      


	if( debug_print1 ) std::cerr << "\tibest_pose: " << ibest_pose << " End of optimization iterations for ibest_pose. st_tot: " << st_tot << "\tfst: " << fst << "\t nbr_st: " << nbr_st << "\tnbr_ct: " << nbr_ct << "\tfst: " << fst << "\tfct: " <<  fct << "\n";

	
   
	double fcomb = 0.0;
	if( considerColorAtoms ){
	  fcomb = opt_param*(fst-fct) + fct;
	}else{
	  fcomb = fst;
	}

	// Found better fcomb_best_pose[ ibest_pose ]
	if( fcomb_best_pose[ ibest_pose ] < fcomb ){
	  fcomb_best_pose[ ibest_pose ] = fcomb;
	}

	// Found a better pose

	if( icycle == 0 ){
	  if ( nbr_comb_prelim < fcomb ) {
	    nbr_comb_prelim = fcomb;
	  }
	}else{
	  if ( nbr_comb < fcomb ) {
	
	    memcpy( best_qom, qom, 7 * sizeof( double ) );
	    nbr_st = fst;
	    nbr_ct = fct;
	    nbr_comb = fcomb;

	    if( debug_print1 ) std::cerr << "\tibest_pose: " << ibest_pose << ":: FOUND BETTER POSE fst: " << fst << "\tfct: " <<  fct << "\tfcomb\t" << fcomb << std::endl;
	  }
	}
      }
    }
  }
      
  // A function originally from PubChem::Data
  VQuaternionTransArray2RotTransMatrix( matrix, best_qom, false );
      
  if( debug_print1 ) std::cerr << "FINAL RESULT: fst: " << nbr_st << " fct: " << nbr_ct << " fcomb: " << nbr_comb << "\n";
}

} // namespace Align3D

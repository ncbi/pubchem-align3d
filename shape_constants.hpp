/*  $Id: shape_constants.hpp 685605 2024-07-26 12:29:33Z thiessen $
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

#ifndef ALIGN3D_SHAPE_CONSTANTS__HPP
#define ALIGN3D_SHAPE_CONSTANTS__HPP

namespace  Align3D {

// Grant-Pickup Parameters

  // Based on Sunghwan Kim's research, a cut-off of 3.6A gives 0.0087,
  //    @ 3.8A gives 0.00095, @ 4.0A gives 0.00051, @ 4.4A gives 0.00013
  //    ST AB overlap contribution between two carbon atoms...
  //    this suggests that a cutoff of 4.5A may be overkill and may be reduced
  //    to a smaller value

  const double default_dCutOff = 4.5;
  extern double dCutOff;
  extern double d2CutOff; 

  // Shape method parameters from Grant, Gallardo and Pickup paper
  // A Fast Method of Molecular Shape Comparison: 
  // A Simple Application of A Gaussian Description of Molecular Shape
  // const double  p = 2.0 * sqrt( 2.0 );  // p = 2.828427..

  const double pi = 3.14159265358979323846;
  const double _p_ = 2.7;
  const double _p2_ = _p_ * _p_;
  const double lambda = ( 4.0 * pi ) / ( 3.0 * _p_ );      //const double lambda = 1.5514037795505151794877251275454; // 4*PI/(3*p);
  const double kappa = pi / pow( lambda, ( 2.0 / 3.0 ) );  //const double  kappa = 2.344228685405120246556971866528f; // PI / ( (lambda)^(2/3) );


 // Grid matrix rotation extent titrated to maximize neighbor count using validation suite
  //        const float  rot1 = 0.980198f, rot2 = 0.19802f;  const double  qrot = 0.100;  // 0.10  3163642
  //        const float  rot1 = 0.976089f, rot2 = 0.21737f;  // 0.11  3171031
  //        const float  rot1 = 0.95599f, rot2 = 0.293399f;  // 0.15  3193121
  //        const float  rot1 = 0.94057f, rot2 = 0.3396f;    // 0.175  3200725
  //        const float  rot1 = 0.923077f, rot2 = 0.384615f; // 0.200  3205001
  //        const float  rot1 = 0.915525f, rot2 = 0.402260f; // 0.210
  const float  rot1 = 0.911634686611f, rot2 = 0.411001457621f;  // 0.215  
  const double  drot1 = 0.911634686611, drot2 = 0.411001457621;  // 0.215  
  //        const float  rot1 = 0.907669f, rot2 = 0.419687f; // 0.220  
  //        const float  rot1 = 0.903629f, rot2 = 0.428316f; // 0.225  3205614
  //        const float  rot1 = 0.899516f, rot2 = 0.436889f; // 0.230
  //        const float  rot1 = 0.895330f, rot2 = 0.445403f; // 0.235

  
  const float  fqrot1 = 0.977659114061f, fqrot = 0.210196709523f;  // 0.215 (un-normalized)
  const float  pqmatrix_type[4][7][12] = {
    // Rot/Trans matrices for 1st four starting points with +/- rotation permutations
    //   correspond to the three zeroed quaternions with +/-0.215 change (later normalized)
    {
      {   1.0,    0.0,    0.0,    0.0, 0.0, 0.0, 0.0 },  //  0   X,  Y,  Z
      { fqrot1,   fqrot,    0.0,    0.0, 0.0, 0.0, 0.0 },
      { fqrot1,  -fqrot,    0.0,    0.0, 0.0, 0.0, 0.0 },
      { fqrot1,    0.0,   fqrot,    0.0, 0.0, 0.0, 0.0 },
      { fqrot1,    0.0,  -fqrot,    0.0, 0.0, 0.0, 0.0 },
      { fqrot1,    0.0,    0.0,   fqrot, 0.0, 0.0, 0.0 },
      { fqrot1,    0.0,    0.0,  -fqrot, 0.0, 0.0, 0.0 }
    }, {
      {   0.0,    1.0,    0.0,    0.0, 0.0, 0.0, 0.0 },  //  1   X, -Y, -Z
      {  fqrot,  fqrot1,    0.0,    0.0, 0.0, 0.0, 0.0 },
      {  fqrot, -fqrot1,    0.0,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,  fqrot1,   fqrot,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,  fqrot1,  -fqrot,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,  fqrot1,    0.0,   fqrot, 0.0, 0.0, 0.0 },
      {   0.0,  fqrot1,    0.0,  -fqrot, 0.0, 0.0, 0.0 }
    }, {
      {   0.0,    0.0,    0.0,    1.0, 0.0, 0.0, 0.0 },  //  2  -X, -Y,  Z
      {  fqrot,    0.0,    0.0,  fqrot1, 0.0, 0.0, 0.0 },
      {  fqrot,    0.0,    0.0, -fqrot1, 0.0, 0.0, 0.0 },
      {   0.0,   fqrot,    0.0,  fqrot1, 0.0, 0.0, 0.0 },
      {   0.0,  -fqrot,    0.0,  fqrot1, 0.0, 0.0, 0.0 },
      {   0.0,    0.0,   fqrot,  fqrot1, 0.0, 0.0, 0.0 },
      {   0.0,    0.0,  -fqrot,  fqrot1, 0.0, 0.0, 0.0 }
    }, {
      {   0.0,    0.0,    1.0,    0.0, 0.0, 0.0, 0.0 },  //  3  -X,  Y, -Z
      {  fqrot,    0.0,  fqrot1,    0.0, 0.0, 0.0, 0.0 },
      {  fqrot,    0.0, -fqrot1,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,   fqrot,  fqrot1,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,  -fqrot,  fqrot1,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,    0.0,  fqrot1,   fqrot, 0.0, 0.0, 0.0 },
      {   0.0,    0.0,  fqrot1,  -fqrot, 0.0, 0.0, 0.0 }
    }
  };

    // Quaternion form of the "matrix_type" matrices
  const float  qmvalf = 0.707106781187f;
  const float  qmatrix_type[16][7] = {
    { 1.0,       0.0,       0.0,       0.0,      0, 0, 0 },  //  0   X,  Y,  Z
    { 0.0,       1.0,       0.0,       0.0,      0, 0, 0 },  //  1   X, -Y, -Z
    { 0.0,       0.0,       0.0,       1.0,      0, 0, 0 },  //  2  -X, -Y,  Z
    { 0.0,       0.0,       1.0,       0.0,      0, 0, 0 },  //  3  -X,  Y, -Z
    { qmvalf,   -qmvalf,    0.0,       0.0,      0, 0, 0 },  //  4   X,  Z, -Y
    { qmvalf,    qmvalf,    0.0,       0.0,      0, 0, 0 },  //  5   X, -Z,  Y
    { 0.0,       0.0,      -qmvalf,    qmvalf,   0, 0, 0 },  //  6  -X, -Z, -Y
    { 0.0,       0.0,       qmvalf,    qmvalf,   0, 0, 0 },  //  7  -X,  Z,  Y
    { qmvalf,    0.0,       0.0,      -qmvalf,   0, 0, 0 },  //  8   Y, -X,  Z
    { 0.0,       qmvalf,    qmvalf,    0.0,      0, 0, 0 },  //  9   Y,  X, -Z
    { qmvalf,    0.0,       0.0,       qmvalf,   0, 0, 0 },  // 10  -Y,  X,  Z
    { 0.0,      -qmvalf,    qmvalf,    0.0,      0, 0, 0 },  // 11  -Y, -X, -Z
    { qmvalf,    0.0,       qmvalf,    0.0,      0, 0, 0 },  // 12   Z,  Y, -X
    { 0.0,       qmvalf,    0.0,       qmvalf,   0, 0, 0 },  // 13   Z, -Y,  X
    { 0.0,      -qmvalf,    0.0,       qmvalf,   0, 0, 0 },  // 14  -Z, -Y, -X
    { qmvalf,    0.0,      -qmvalf,    0.0,      0, 0, 0 }   // 15  -Z,  Y,  X
  };
  const double  qmval = 0.707106781187;
  const double  qmatrix_typed[16][7] = {
    { 1.0,       0.0,       0.0,       0.0,      0, 0, 0 },  //  0   X,  Y,  Z
    { 0.0,       1.0,       0.0,       0.0,      0, 0, 0 },  //  1   X, -Y, -Z
    { 0.0,       0.0,       0.0,       1.0,      0, 0, 0 },  //  2  -X, -Y,  Z
    { 0.0,       0.0,       1.0,       0.0,      0, 0, 0 },  //  3  -X,  Y, -Z
    { qmval,    -qmval,     0.0,       0.0,      0, 0, 0 },  //  4   X,  Z, -Y
    { qmval,     qmval,     0.0,       0.0,      0, 0, 0 },  //  5   X, -Z,  Y
    { 0.0,       0.0,      -qmval,     qmval,    0, 0, 0 },  //  6  -X, -Z, -Y
    { 0.0,       0.0,       qmval,     qmval,    0, 0, 0 },  //  7  -X,  Z,  Y
    { qmval,     0.0,       0.0,      -qmval,    0, 0, 0 },  //  8   Y, -X,  Z
    { 0.0,       qmval,     qmval,     0.0,      0, 0, 0 },  //  9   Y,  X, -Z
    { qmval,     0.0,       0.0,       qmval,    0, 0, 0 },  // 10  -Y,  X,  Z
    { 0.0,      -qmval,     qmval,     0.0,      0, 0, 0 },  // 11  -Y, -X, -Z
    { qmval,     0.0,       qmval,     0.0,      0, 0, 0 },  // 12   Z,  Y, -X
    { 0.0,       qmval,     0.0,       qmval,    0, 0, 0 },  // 13   Z, -Y,  X
    { 0.0,      -qmval,     0.0,       qmval,    0, 0, 0 },  // 14  -Z, -Y, -X
    { qmval,     0.0,      -qmval,     0.0,      0, 0, 0 }   // 15  -Z,  Y,  X
  };

    // Grid symmetry transformations for starting points
  const float  matrix_type[16][12] = {
    //  Rxx   Rxy   Rxz   Ryx   Ryy   Ryz   Rzx   Rzy   Rzz   Tx    Ty    Tz
    { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f },  //  0   X,  Y,  Z
    { 1.0f, 0.0f, 0.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f },  //  1   X, -Y, -Z
    {-1.0f, 0.0f, 0.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f },  //  2  -X, -Y,  Z
    {-1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f },  //  3  -X,  Y, -Z

    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 0.0f },  //  4   X,  Z, -Y
    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,-1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f },  //  5   X, -Z,  Y
    {-1.0f, 0.0f, 0.0f, 0.0f, 0.0f,-1.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 0.0f },  //  6  -X, -Z, -Y
    {-1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f },  //  7  -X,  Z,  Y

    { 0.0f, 1.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f },  //  8   Y, -X,  Z
    { 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f },  //  9   Y,  X, -Z
    { 0.0f,-1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f },  // 10  -Y,  X,  Z
    { 0.0f,-1.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f },  // 11  -Y, -X, -Z

    { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f },  // 12   Z,  Y, -X
    { 0.0f, 0.0f, 1.0f, 0.0f,-1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f },  // 13   Z, -Y,  X
    { 0.0f, 0.0f,-1.0f, 0.0f,-1.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f },  // 14  -Z, -Y, -X
    { 0.0f, 0.0f,-1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f }   // 15  -Z,  Y,  X
  };
  const double  matrix_typed[16][12] = {
    //  Rxx  Rxy  Rxz  Ryx  Ryy  Ryz  Rzx  Rzy  Rzz  Tx   Ty   Tz
    { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },  //  0   X,  Y,  Z
    { 1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0 },  //  1   X, -Y, -Z
    {-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },  //  2  -X, -Y,  Z
    {-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0 },  //  3  -X,  Y, -Z

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0 },  //  4   X,  Z, -Y
    { 1.0, 0.0, 0.0, 0.0, 0.0,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },  //  5   X, -Z,  Y
    {-1.0, 0.0, 0.0, 0.0, 0.0,-1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0 },  //  6  -X, -Z, -Y
    {-1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },  //  7  -X,  Z,  Y

    { 0.0, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },  //  8   Y, -X,  Z
    { 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0 },  //  9   Y,  X, -Z
    { 0.0,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },  // 10  -Y,  X,  Z
    { 0.0,-1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0 },  // 11  -Y, -X, -Z

    { 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },  // 12   Z,  Y, -X
    { 0.0, 0.0, 1.0, 0.0,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },  // 13   Z, -Y,  X
    { 0.0, 0.0,-1.0, 0.0,-1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },  // 14  -Z, -Y, -X
    { 0.0, 0.0,-1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 }   // 15  -Z,  Y,  X
  };

  const double  qrot1 = 0.977659114061, qrot = 0.210196709523;  // 0.215 (un-normalized)
  const double  pqmatrix_typed[4][7][12] = {
    // Rot/Trans matrices for 1st four starting points with +/- rotation permutations
    //   correspond to the three zeroed quaternions with +/-0.215 change (later normalized)
    {
      {   1.0,    0.0,    0.0,    0.0, 0.0, 0.0, 0.0 },  //  0   X,  Y,  Z
      { qrot1,   qrot,    0.0,    0.0, 0.0, 0.0, 0.0 },
      { qrot1,  -qrot,    0.0,    0.0, 0.0, 0.0, 0.0 },
      { qrot1,    0.0,   qrot,    0.0, 0.0, 0.0, 0.0 },
      { qrot1,    0.0,  -qrot,    0.0, 0.0, 0.0, 0.0 },
      { qrot1,    0.0,    0.0,   qrot, 0.0, 0.0, 0.0 },
      { qrot1,    0.0,    0.0,  -qrot, 0.0, 0.0, 0.0 }
    }, {
      {   0.0,    1.0,    0.0,    0.0, 0.0, 0.0, 0.0 },  //  1   X, -Y, -Z
      {  qrot,  qrot1,    0.0,    0.0, 0.0, 0.0, 0.0 },
      {  qrot, -qrot1,    0.0,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,  qrot1,   qrot,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,  qrot1,  -qrot,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,  qrot1,    0.0,   qrot, 0.0, 0.0, 0.0 },
      {   0.0,  qrot1,    0.0,  -qrot, 0.0, 0.0, 0.0 }
    }, {
      {   0.0,    0.0,    0.0,    1.0, 0.0, 0.0, 0.0 },  //  2  -X, -Y,  Z
      {  qrot,    0.0,    0.0,  qrot1, 0.0, 0.0, 0.0 },
      {  qrot,    0.0,    0.0, -qrot1, 0.0, 0.0, 0.0 },
      {   0.0,   qrot,    0.0,  qrot1, 0.0, 0.0, 0.0 },
      {   0.0,  -qrot,    0.0,  qrot1, 0.0, 0.0, 0.0 },
      {   0.0,    0.0,   qrot,  qrot1, 0.0, 0.0, 0.0 },
      {   0.0,    0.0,  -qrot,  qrot1, 0.0, 0.0, 0.0 }
    }, {
      {   0.0,    0.0,    1.0,    0.0, 0.0, 0.0, 0.0 },  //  3  -X,  Y, -Z
      {  qrot,    0.0,  qrot1,    0.0, 0.0, 0.0, 0.0 },
      {  qrot,    0.0, -qrot1,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,   qrot,  qrot1,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,  -qrot,  qrot1,    0.0, 0.0, 0.0, 0.0 },
      {   0.0,    0.0,  qrot1,   qrot, 0.0, 0.0, 0.0 },
      {   0.0,    0.0,  qrot1,  -qrot, 0.0, 0.0, 0.0 }
    }
  };

    const double  pmatrix_typed[4][7][12] = {
    // Rot/Trans matrices for 1st four starting points with +/- rotation permutations
    //   correspond to the three zeroed quaternions with +/-0.225 change (later normalized)
    {
      { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,0.0,0.0 },  //  0   X,  Y,  Z
      { 1.0, 0.0, 0.0, 0.0, drot1,-drot2, 0.0, drot2, drot1, 0.0,0.0,0.0 },   // 1, -.225,   0,   0  {Y-Z Z+Y}
      { 1.0, 0.0, 0.0, 0.0, drot1, drot2, 0.0,-drot2, drot1, 0.0,0.0,0.0 }, // 1, +.225,   0,   0  {Y+Z  Z-Y}
      { drot1, 0.0, drot2, 0.0, 1.0, 0.0,-drot2, 0.0, drot1, 0.0,0.0,0.0 },   // 1,   0, -.225,   0  {X+Z  Z-X}
      { drot1, 0.0,-drot2, 0.0, 1.0, 0.0, drot2, 0.0, drot1, 0.0,0.0,0.0 }, // 1,   0, +.225,   0  {X-Z  Z+X}
      { drot1,-drot2, 0.0, drot2, drot1, 0.0, 0.0, 0.0, 1.0, 0.0,0.0,0.0 },   // 1,   0,   0, -.225  {X-Y  Y+X}
      { drot1, drot2, 0.0,-drot2, drot1, 0.0, 0.0, 0.0, 1.0, 0.0,0.0,0.0 }  // 1,   0,   0, +.225  {X+Y  Y-X}
    }, {
      { 1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0,-1.0, 0.0,0.0,0.0 },  //  1   X, -Y, -Z
      { 1.0, 0.0, 0.0, 0.0,-drot1,-drot2, 0.0, drot2,-drot1, 0.0,0.0,0.0 },
      { 1.0, 0.0, 0.0, 0.0,-drot1, drot2, 0.0,-drot2,-drot1, 0.0,0.0,0.0 },
      { drot1, drot2, 0.0, drot2,-drot1, 0.0, 0.0, 0.0,-1.0, 0.0,0.0,0.0 },
      { drot1,-drot2, 0.0,-drot2,-drot1, 0.0, 0.0, 0.0,-1.0, 0.0,0.0,0.0 },
      { drot1, 0.0, drot2, 0.0,-1.0, 0.0, drot2, 0.0,-drot1, 0.0,0.0,0.0 },
      { drot1, 0.0,-drot2, 0.0,-1.0, 0.0,-drot2, 0.0,-drot1, 0.0,0.0,0.0 }
    }, {
      {-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0, 1.0, 0.0,0.0,0.0 },  //  2  -X, -Y,  Z
      {-drot1,-drot2, 0.0, drot2,-drot1, 0.0, 0.0, 0.0, 1.0, 0.0,0.0,0.0 },
      {-drot1, drot2, 0.0,-drot2,-drot1, 0.0, 0.0, 0.0, 1.0, 0.0,0.0,0.0 },
      {-drot1, 0.0, drot2, 0.0,-1.0, 0.0, drot2, 0.0, drot1, 0.0,0.0,0.0 },
      {-drot1, 0.0,-drot2, 0.0,-1.0, 0.0,-drot2, 0.0, drot1, 0.0,0.0,0.0 },
      {-1.0, 0.0, 0.0, 0.0,-drot1, drot2, 0.0, drot2, drot1, 0.0,0.0,0.0 },
      {-1.0, 0.0, 0.0, 0.0,-drot1,-drot2, 0.0,-drot2, drot1, 0.0,0.0,0.0 }
    }, {
      {-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,-1.0, 0.0,0.0,0.0 },  //  3  -X,  Y, -Z
      {-drot1, 0.0, drot2, 0.0, 1.0, 0.0,-drot2, 0.0,-drot1, 0.0,0.0,0.0 },
      {-drot1, 0.0,-drot2, 0.0, 1.0, 0.0, drot2, 0.0,-drot1, 0.0,0.0,0.0 },
      {-drot1, drot2, 0.0, drot2, drot1, 0.0, 0.0, 0.0,-1.0, 0.0,0.0,0.0 },
      {-drot1,-drot2, 0.0,-drot2, drot1, 0.0, 0.0, 0.0,-1.0, 0.0,0.0,0.0 },
      {-1.0, 0.0, 0.0, 0.0, drot1, drot2, 0.0, drot2,-drot1, 0.0,0.0,0.0 },
      {-1.0, 0.0, 0.0, 0.0, drot1,-drot2, 0.0,-drot2,-drot1, 0.0,0.0,0.0 }
    }
  };

}

#endif // ALIGN3D_SHAPE_CONSTANTS__HPP

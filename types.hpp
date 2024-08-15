/*  $Id: types.hpp 685605 2024-07-26 12:29:33Z thiessen $
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

#ifndef ALIGN3D_TYPES__HPP
#define ALIGN3D_TYPES__HPP

#if (__SIZEOF_LONG__ == 8)
    typedef signed long Int8;    
    typedef unsigned long Uint8;
#elif (__SIZEOF_LONG_LONG__ == 8)
    typedef signed long long Int8;    
    typedef unsigned long long Uint8;
#else
    #error "This platform does not support 8-byte integer"
#endif

#endif // ALIGN3D_TYPES__HPP

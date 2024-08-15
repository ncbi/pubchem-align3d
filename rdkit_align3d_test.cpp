/*  $Id: rdkit_align3d_test.cpp 685605 2024-07-26 12:29:33Z thiessen $
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

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>

#include "shape_functions.hpp"

using namespace std;

#define ERRORTHROW(msg_stream) do { \
    ostringstream os; \
    os << msg_stream; \
    throw runtime_error(os.str()); \
} while (0)

#define DEBUG_MSG(msg_stream) cout << msg_stream << '\n'
//#define DEBUG_MSG(msg_stream)

using namespace RDKit;

// Bondi radii
const double radius_carbon      = 1.70;
const double radius_nitrogen    = 1.55;
const double radius_oxygen      = 1.52;
const double radius_fluorine    = 1.47;
const double radius_silicon     = 2.10;
const double radius_phosphorous = 1.80;
const double radius_sulfur      = 1.80;
const double radius_chlorine    = 1.75;
const double radius_bromine     = 1.85;
const double radius_iodine      = 1.98;
const double radius_color       = 1.08265;  // same radius for all feature/color "atoms"

void PrepareConformer(
    RWMol& mol,
    vector < float >& coord,
    vector < double >& alpha_vector,
    vector < unsigned int >& atom_type_vector,
    vector < unsigned int >& volumeAtomIndexVector,
    map < unsigned int, vector < unsigned int > >& colorAtomType2IndexVectorMap,
    double& sov,
    double& sof)
{
    coord.clear();
    alpha_vector.clear();
    atom_type_vector.clear();
    volumeAtomIndexVector.clear();
    colorAtomType2IndexVectorMap.clear();
    sov = 0.0;
    sof = 0.0;
    
    // unpack features (PubChem-specific property from SDF)
    // NOTE: this unpacking assumes that RWMol-atom-index = SDF-atom-number - 1 
    //   e.g. RWMol uses [0..N-1] and SDF uses [1..N], with atoms in the same order
    
    vector < pair < vector < unsigned int >, unsigned int > > feature_idx_type;
    
    if (mol.hasProp("PUBCHEM_PHARMACOPHORE_FEATURES")) {
    
        // regular atoms have type 0; feature "atoms" (features represented by a single point+radius) must have type > 0
        static const map < string, unsigned int > atomTypes = {
            { "acceptor", 1 },
            { "anion", 2 },
            { "cation", 3 },
            { "donor", 4 },
            { "hydrophobe", 5 },
            { "rings", 6 },
        };

        string features;
        mol.getProp("PUBCHEM_PHARMACOPHORE_FEATURES", features);
        istringstream iss(features);
        string line;
        unsigned int n = 0;
        while (getline(iss, line)) {
            
            if (n == 0) {
                feature_idx_type.resize(stoul(line));
            } 
            
            else {
                unsigned int f = n - 1;
                if (f >= feature_idx_type.size())
                    ERRORTHROW("Too many features");
                    
                istringstream iss2(line);
                string token;
                unsigned int t = 0;
                while (getline(iss2, token, ' ')) {
                    if (t == 0) {
                        feature_idx_type[f].first.resize(stoul(token));
                    } else if (t <= feature_idx_type[f].first.size()) {
                        feature_idx_type[f].first[t - 1] = stoul(token) - 1;
                    } else {
                        map < string, unsigned int >::const_iterator type = atomTypes.find(token);
                        if (type == atomTypes.end())
                            ERRORTHROW("Invalid feature type");
                        feature_idx_type[f].second = type->second;
                    }
                    ++t;
                }
                if (t != (feature_idx_type[f].first.size() + 2))
                    ERRORTHROW("Wrong number of tokens in feature");
            }
            
            ++n;
        }
        if (n != (feature_idx_type.size() + 1))
            ERRORTHROW("Wrong number of features");
        
        DEBUG_MSG("# features: " << feature_idx_type.size());
    }

    // unpack atoms
    
    if (mol.getNumConformers() != 1)
        ERRORTHROW("Must have exactly one conformer");
    Conformer& conformer = mol.getConformer();
    if (!conformer.is3D())
        ERRORTHROW("Conformer must be 3D");
        
    unsigned int nAtoms = mol.getNumAtoms();
    //DEBUG_MSG("num atoms: " << nAtoms);
    
    unsigned int nHeavyAtoms = 0;
    unsigned int i;
    for (i=0; i<nAtoms; ++i)
        if (mol.getAtomWithIdx(i)->getAtomicNum() > 1)
            ++nHeavyAtoms;
    DEBUG_MSG("num heavy atoms: " << nHeavyAtoms);
    
    unsigned int nAlignmentAtoms = nHeavyAtoms + feature_idx_type.size();
    
    vector < double > rad_vector(nAlignmentAtoms);
    atom_type_vector.resize(nAlignmentAtoms);
    
    unsigned int array_idx = 0;
    float ave_x = 0.0, ave_y = 0.0, ave_z = 0.0; 
    for (i=0; i<nAtoms; ++i) {
        unsigned int Z = mol.getAtomWithIdx(i)->getAtomicNum();
        if (Z > 1) {
        
            atom_type_vector[array_idx] = 0;
            
            const RDGeom::Point3D& pos = conformer.getAtomPos(i);
            ave_x += pos.x;
            ave_y += pos.y;
            ave_z += pos.z;
        
            switch (Z) {
                case 6:  rad_vector[array_idx] = radius_carbon;         break;
                case 7:  rad_vector[array_idx] = radius_nitrogen;       break;
                case 8:  rad_vector[array_idx] = radius_oxygen;         break;
                case 9:  rad_vector[array_idx] = radius_fluorine;       break;
                case 14: rad_vector[array_idx] = radius_silicon;        break;
                case 15: rad_vector[array_idx] = radius_phosphorous;    break;
                case 16: rad_vector[array_idx] = radius_sulfur;         break;
                case 17: rad_vector[array_idx] = radius_chlorine;       break;
                case 35: rad_vector[array_idx] = radius_bromine;        break;
                case 53: rad_vector[array_idx] = radius_iodine;         break;
                default:
                    ERRORTHROW("Can't use molecules with element Z=" << Z);
            }
            
            ++array_idx;
        }
    }

    // translate steric center to origin
    ave_x /= nHeavyAtoms;
    ave_y /= nHeavyAtoms;
    ave_z /= nHeavyAtoms;
    DEBUG_MSG("steric center: (" << ave_x << ", " << ave_y << ", " << ave_z << ")");

    coord.resize(nAlignmentAtoms * 3);

    array_idx = 0;
    for (i=0; i<nAtoms; ++i) {
    
        // translate all atoms
        RDGeom::Point3D& pos = conformer.getAtomPos(i);
        pos.x -= ave_x;
        pos.y -= ave_y;
        pos.z -= ave_z;
        
        // but use only non-H for alignment
        if (mol.getAtomWithIdx(i)->getAtomicNum() > 1) {        
            coord[array_idx * 3] = pos.x;
            coord[(array_idx * 3) + 1] = pos.y;
            coord[(array_idx * 3) + 2] = pos.z;
            ++array_idx;
        }
    }

    // get feature coordinates - simply the average of coords of all atoms in the feature
    for (i=0; i<feature_idx_type.size(); ++i) {
    
        double x = 0.0, y = 0.0, z = 0.0;
        for (unsigned int j=0; j<feature_idx_type[i].first.size(); ++j) {
            unsigned int idx = feature_idx_type[i].first[j];
            if (idx >= nAtoms || mol.getAtomWithIdx(idx)->getAtomicNum() <= 1)
                ERRORTHROW("Invalid feature atom index");
            const RDGeom::Point3D& pos = conformer.getAtomPos(idx);
            x += pos.x;
            y += pos.y;
            z += pos.z;
        }
        x /= feature_idx_type[i].first.size();
        y /= feature_idx_type[i].first.size();
        z /= feature_idx_type[i].first.size();
        DEBUG_MSG("feature type " << feature_idx_type[i].second << " (" << x << ", " << y << ", " << z << ")");
        
        array_idx = nHeavyAtoms + i;
        coord[array_idx * 3] = x;
        coord[(array_idx * 3) + 1] = y;
        coord[(array_idx * 3) + 2] = z;
        rad_vector[array_idx] = radius_color;
        atom_type_vector[array_idx] = feature_idx_type[i].second;
    }
    
    Align3D::setAlpha(&(rad_vector[0]), rad_vector.size(), alpha_vector);
        
    // regular atom self overlap
    Align3D::getVolumeAtomIndexVector(&(atom_type_vector[0]), atom_type_vector.size(), volumeAtomIndexVector);
    sov = Align3D::ComputeShapeOverlap(
        &(coord[0]),
        alpha_vector,
        volumeAtomIndexVector,
        &(coord[0]),
        alpha_vector,
        volumeAtomIndexVector);
    DEBUG_MSG("sov: " << sov);
        
    // feature self overlap
    if (feature_idx_type.size() > 0) {
        Align3D::getColorAtomType2IndexVectorMap(&(atom_type_vector[0]), atom_type_vector.size(), colorAtomType2IndexVectorMap);
        sof = Align3D::ComputeFeatureOverlap(
	    &(coord[0]),
	    alpha_vector,
	    colorAtomType2IndexVectorMap,
	    &(coord[0]),
	    alpha_vector,
	    colorAtomType2IndexVectorMap);
        DEBUG_MSG("sof: " << sof);
    }
}

void Neighbor(
    // inputs
    RWMol& ref,
    RWMol& fit,
    const double opt_param,
    const unsigned int max_preiters,
    const unsigned int max_postiters,
    // outputs
    float* matrix,
    double& nbr_st,
    double& nbr_ct)
{
    vector < float > ref_coord;
    vector < double > alpha_ref_vector;
    vector < unsigned int > ref_volumeAtomIndexVector;
    map < unsigned int, vector < unsigned int > > ref_colorAtomType2IndexVectorMap;
    double ref_sov;
    double ref_sof;
    vector < unsigned int > ref_atom_type_vector;

    vector < float > fit_coord;
    vector < double > alpha_fit_vector;
    vector < unsigned int > fit_volumeAtomIndexVector;
    map < unsigned int, vector < unsigned int > > fit_colorAtomType2IndexVectorMap;
    double fit_sov;
    double fit_sof;
    vector < unsigned int > fit_atom_type_vector;

    Align3D::setUseCutOff(true);

    DEBUG_MSG("Reference details:");
    PrepareConformer(
        ref, 
        ref_coord, 
        alpha_ref_vector, 
        ref_atom_type_vector,
        ref_volumeAtomIndexVector, 
        ref_colorAtomType2IndexVectorMap, 
        ref_sov,
        ref_sof);
        
    DEBUG_MSG("Fit details:");
    PrepareConformer(
        fit, 
        fit_coord, 
        alpha_fit_vector, 
        fit_atom_type_vector,
        fit_volumeAtomIndexVector, 
        fit_colorAtomType2IndexVectorMap, 
        fit_sov,
        fit_sof);
        
    // write out ref and fit, both translated to origin but not yet aligned 
    {{
        SDWriter ref_sdf("ref.sdf");
        ref_sdf.write(ref);
        SDWriter fit_sdf("fit.sdf");
        fit_sdf.write(fit);
    }}

    set < unsigned int > jointColorAtomTypeSet;
    Align3D::getJointColorTypeSet( 
        &(ref_atom_type_vector[0]),
	ref_atom_type_vector.size(),
	&(fit_atom_type_vector[0]),
	fit_atom_type_vector.size(),
	jointColorAtomTypeSet);
    Align3D::restrictColorAtomType2IndexVectorMap(ref_colorAtomType2IndexVectorMap, jointColorAtomTypeSet);
    Align3D::restrictColorAtomType2IndexVectorMap(fit_colorAtomType2IndexVectorMap, jointColorAtomTypeSet);

    DEBUG_MSG("Running alignment...");
    Align3D::Neighbor_Conformers(
        &(ref_coord[0]),
        alpha_ref_vector,
        ref_volumeAtomIndexVector,
        ref_colorAtomType2IndexVectorMap,
        ref_sov,
        ref_sof,
        &(fit_coord[0]),
        alpha_fit_vector,
        fit_volumeAtomIndexVector,
        fit_colorAtomType2IndexVectorMap,
        fit_sov,
        fit_sof,
        !jointColorAtomTypeSet.empty(),
        max_preiters,
        max_postiters,
        opt_param,
        matrix,
        nbr_st,
        nbr_ct);
        
    DEBUG_MSG("Done!");
    DEBUG_MSG("nbr_st: " << nbr_st);
    DEBUG_MSG("nbr_ct: " << nbr_ct);
    
    // transform fit coords
    Conformer& fit_conformer = fit.getConformer();
    vector < float > orig(fit.getNumAtoms() * 3), transformed(fit.getNumAtoms() * 3);
    unsigned int i;
    for (i=0; i<fit.getNumAtoms(); ++i) {
        const RDGeom::Point3D& pos = fit_conformer.getAtomPos(i);
        orig[i * 3] = pos.x;
        orig[(i * 3) + 1] = pos.y;
        orig[(i * 3) + 2] = pos.z;
    }
    
    Align3D::VApplyRotTransMatrix(&(transformed[0]), &(orig[0]), fit.getNumAtoms(), matrix);
    
    for (i=0; i<fit.getNumAtoms(); ++i) {
        RDGeom::Point3D& pos = fit_conformer.getAtomPos(i);
        pos.x = transformed[i * 3];
        pos.y = transformed[(i * 3) + 1];
        pos.z = transformed[(i * 3) + 2];
    }

    SDWriter fit_sdf("fit-aligned.sdf");
    fit_sdf.write(fit);
}

int main(int argc, char **argv)
{
    int status = -1;
    
    try {
    
        if (argc != 6)
            ERRORTHROW("Usage: rdkit_example <ref_conformer.sdf> <fit_conformer.sdf> opt_param max_preiters max_postiters");
        
        SDMolSupplier s_ref(argv[1], false, false, true);
        unique_ptr < ROMol > ro_ref(s_ref.next());
        if (!ro_ref.get())
            ERRORTHROW("Failed to read ref conformer");
        RWMol ref(*ro_ref);
        
        SDMolSupplier s_fit(argv[2], false, false, true);
        unique_ptr < ROMol > ro_fit(s_fit.next());
        if (!ro_fit.get())
            ERRORTHROW("Failed to read fit conformer");
        RWMol fit(*ro_fit);
        
        vector < float > matrix(12, 0.0);
        double nbr_st = 0.0;
        double nbr_ct = 0.0;
        Neighbor(ref, fit, stod(argv[3]), stoul(argv[4]), stoul(argv[5]), &(matrix[0]), nbr_st, nbr_ct);
        
        status = 0;
        
    } catch (std::exception& e) {
        cerr << "Caught std::exception: " << e.what() << '\n';
    } catch (...) {
        cerr << "Caught unknown exception\n";
    }
    
    return status;
}

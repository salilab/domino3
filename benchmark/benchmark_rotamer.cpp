/**
 * Copyright 2007-2017 IMP Inventors. All rights reserved.
 */

#include <IMP/flags.h>
#include <IMP/Model.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/hierarchy_tools.h>
#include <IMP/atom/pdb.h>



#include <IMP/rmf/atom_io.h>
#include <IMP/rmf/restraint_io.h>
#include <IMP/rmf/frames.h>
#include <IMP/rmf/particle_io.h>

#include <IMP/SetLogState.h>
#include <IMP/domino3/Factor.h>
#include <IMP/domino3/DistanceFactor.h>
#include <IMP/domino3/ExcludedVolumeFactor.h>
#include <IMP/domino3/Updater.h>
#include <RMF/FileHandle.h>
#include <IMP/domino3/XYZStates.h>
#include <IMP/rotamer/RotamerLibrary.h>
#include <IMP/core/XYZ.h>
#include <IMP/multifit/anchors_reader.h>
#include <algorithm>

namespace {
    
        
    std::string rotamer_lib = IMP::domino3::get_example_path("rotamer.lib");
    IMP::AddStringFlag asf("lib", "Input lib anchors", &rotamer_lib);
    std::string pdb_path = IMP::domino3::get_example_path("example_distance_restraints.txt");
    IMP::AddStringFlag oasf("pdb", "Input file restraints", &pdb_path);
    boost::int64_t iterations = 100;
    IMP::AddIntFlag aif("iterations", "Number of iterations",
                              &iterations);
    
    
    
    
    void run_it(IMP::Model *) {
//        IMP::rotamer::RotamerLibrary rl;
//        rl.read_library_file(rotamer_lib);
//        DistantRestraints restraints=read_restraints(restraints_path);
//        IMP::set<int> residues;
//        for( int i = 0; i < restraints.size(); i++){
//            residues.insert(restraints[i]->from);
//            residues.insert(restraints[i]->to);
//        }
//        IMP::ParticlesTemp ps;
//        for( int i = 0; i < vs_org.size(); i++){
//            IMP::Particle * p = new IMP::Particle(m);
//            ps.push_back(p);
//            IMP::core::XYZR::setup_particle(m, p->get_index(),  IMP::algebra::Sphere3D(vs_org[i], 1));
//        }
//        IMP::ParticleIndexes pis = IMP::get_indexes(ps);
//        for( int i = 0; i < pis.size(); i++){
//            std::ostringstream oss;
//            oss  << pis[i];
//            m->get_particle(pis[i])->set_name(oss.str());
//        }
//        //
//        IMP::algebra::Vector3Ds vs(pis.size());
//        for (unsigned int i = 0; i< vs.size(); ++i) {
//            vs[i] = IMP::core::XYZ(m, pis[i]).get_coordinates();
//        }
//        IMP_NEW(IMP::domino3::XYZStates, states, (m, vs));
//        IMP_NEW(IMP::domino3::StatesTable, st, (m));
//        
//        // fill statetable with marginals
//        for (unsigned int i = 0; i <= residues.size(); ++i) {
//            IMP_NEW(IMP::domino3::Marginals, marginals, (m, pis[i], vs.size()));
//            marginals->set_uniform();
//            st->add(pis[i], states, marginals);
//        }
//        
//        
//        IMP::domino3::Factors factors;
//        for( int i = 0; i < restraints.size(); i++){
//            IMP::ParticleIndexPair cur_pair(pis[restraints[i]->from], pis[restraints[i]->to]);
//            IMP_NEW(IMP::domino3::DistanceFactor, dn,(m, cur_pair, restraints[i]->dist, 1, st));
//            factors.push_back(dn);
//            dn->set_was_used(true);
//        }
//
//        std::cout << "Node size: " << factors.size() << std::endl;
//        IMP::domino3::add_neighbors(factors);
//        //     InteractionGraph ig = get_interaction_graph(rs, pst);
//        //      SubsetGraph jt = get_junction_tree(ig);
//        
//        IMP_NEW(IMP::domino3::Updater, ud, (factors, "updater"));
//        
//        std::cout << "before" << std::endl;
//        st->show_marginal();
//        ud->update(iterations);
//        IMP::domino3::update_state_table(factors,st);
//        IMP::domino3::print_graph(factors);
//        std::cout << "after" << std::endl;
//        
//        st->show_marginal();
//        double probability_to_see_best_match = 1;
//        std::vector<int> order;
//        for(int i = 0; i < residues.size(); i++){
//            std::vector<std::pair<double, int> > pos_order;
//            IMP::domino3::Marginals * m  = st->get_marginals_by_order(i);
//            for(int state = 0; state < m->get_number(); state++)
//                pos_order.push_back(std::make_pair( exp(m->get_current_marginal(state)), state));
//            sort(pos_order.begin(),pos_order.end(), std::greater<std::pair<double, int> >());
//            probability_to_see_best_match *= pos_order[0].first;
//            order.push_back(pos_order[0].second);
//        }
//        
//        for(int i = 0 ; i < order.size();i++){
//            std::cout << order[i] << " ";
//        }
//        std::cout << std::endl;
        
        //    std::cout << "Weighted RMSD: " << calc_weighted_rmsd(vs,ps,st) << std::endl;
    }
}



int main(int argc, char **argv) {
    IMP::setup_from_argv(argc, argv, "Experiment with loopy domino");
    IMP_NEW(IMP::Model, m, ());
    //    IMP::atom::Hierarchy pdb = IMP::atom:: (input, m,new IMP::atom::CAlphaPDBSelector());
    

    //  IMP::atom::Atoms atoms
    //    = IMP::atom::get_by_type(pdb, IMP::atom::ATOM_TYPE);
    
    
    run_it(m);
    
    return 0;
}

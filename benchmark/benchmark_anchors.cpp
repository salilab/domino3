/**
 * Copyright 2007-2013 IMP Inventors. All rights reserved.
 */

#include <IMP/base/flags.h>
#include <IMP/kernel/Model.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/hierarchy_tools.h>
#include <IMP/atom/pdb.h>



#include <IMP/rmf/atom_io.h>
#include <IMP/rmf/restraint_io.h>
#include <IMP/rmf/frames.h>
#include <IMP/rmf/particle_io.h>

#include <IMP/base/SetLogState.h>
#include <IMP/domino3/Factor.h>
#include <IMP/domino3/DistanceFactor.h>
#include <IMP/domino3/ExcludedVolumeFactor.h>
#include <IMP/domino3/Updater.h>
#include <RMF/FileHandle.h>
#include <IMP/domino3/XYZStates.h>
#include <IMP/core/XYZ.h>
#include <IMP/multifit/anchors_reader.h>
#include <algorithm>

namespace {
    
    
    struct DistantRestraint{
        int from;
        int to;
        double dist;
        DistantRestraint( int pfrom,int pto, float pdist):from(pfrom),to(pto),dist(pdist){}
        
    };
    typedef  std::vector<DistantRestraint * > DistantRestraints;
    DistantRestraints read_restraints(std::string path){
        std::cout << "Read chem sim. file: " << path << std::endl;
        std::ifstream infile(path.c_str());
        std::string line;
        DistantRestraints ret_array;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            int  trash, from, to;
            double dist;
            iss >> trash >> from >> to >> dist;
            
            ret_array.push_back(new DistantRestraint(from, to, dist) );
            
        }
        return ret_array;
    }
    
    
    double calc_weighted_rmsd(IMP::algebra::Vector3Ds &vs, const IMP::kernel::ParticlesTemp &ps,IMP::domino3::StatesTable *pst){
        IMP::kernel::ParticleIndexes pis = IMP::kernel::get_indexes(ps);
        double ret_sum=0.0;
        for (unsigned int i = 0; i < ps.size(); ++i) {
            IMP::domino3::Marginals * marg = pst->get_marginals(pis[i]);
            IMP::algebra::Vector3D vector_i = vs[i];
            for(int y = 0; y < marg->get_number(); y++){
                IMP::algebra::Vector3D vector_y = vs[y];
                double dist=IMP::algebra::get_distance(vector_i, vector_y);
                ret_sum+=exp(marg->get_current_marginal(y))*dist;
            }
        }
        return ret_sum;
    }
    
    
    std::string anchors = IMP::domino3::get_example_path("example_distance_anchors.txt");
    IMP::base::AddStringFlag asf("anchors", "Input file anchors", &anchors);
    std::string restraints_path = IMP::domino3::get_example_path("example_distance_restraints.txt");
    IMP::base::AddStringFlag oasf("restraints", "Input file restraints", &restraints_path);
    double dist = std::numeric_limits<double>::infinity();
    IMP::base::AddFloatFlag dasf("distance", "Maximum distance for restraints",
                                 &dist);
    boost::int64_t iterations = 100;
    IMP::base::AddIntFlag aif("iterations", "Number of iterations",
                              &iterations);
    
    
    
    
    void run_it(IMP::kernel::Model *m, const IMP::algebra::Vector3Ds &vs_org) {
        DistantRestraints restraints=read_restraints(restraints_path);
        IMP::base::set<int> residues;
        for( int i = 0; i < restraints.size(); i++){
            residues.insert(restraints[i]->from);
            residues.insert(restraints[i]->to);
        }
        IMP::ParticlesTemp ps;
        for( int i = 0; i < vs_org.size(); i++){
            IMP::kernel::Particle * p = new IMP::kernel::Particle(m);
            ps.push_back(p);
            IMP::core::XYZR::setup_particle(m, p->get_index(),  IMP::algebra::Sphere3D(vs_org[i], 1));
        }
        IMP::kernel::ParticleIndexes pis = IMP::kernel::get_indexes(ps);
        for( int i = 0; i < pis.size(); i++){
            std::ostringstream oss;
            oss  << pis[i];
            m->get_particle(pis[i])->set_name(oss.str());
        }
        //
        IMP::algebra::Vector3Ds vs(pis.size());
        for (unsigned int i = 0; i< vs.size(); ++i) {
            vs[i] = IMP::core::XYZ(m, pis[i]).get_coordinates();
        }
        IMP_NEW(IMP::domino3::XYZStates, states, (m, vs));
        IMP_NEW(IMP::domino3::StatesTable, st, (m));
        
        // fill statetable with marginals
        for (unsigned int i = 0; i <= residues.size(); ++i) {
            IMP_NEW(IMP::domino3::Marginals, marginals, (m, pis[i], vs.size()));
            marginals->set_uniform();
            st->add(pis[i], states, marginals);
        }
        
        
        IMP::domino3::Factors factors;
        for( int i = 0; i < restraints.size(); i++){
            IMP::kernel::ParticleIndexPair cur_pair(pis[restraints[i]->from], pis[restraints[i]->to]);
            IMP_NEW(IMP::domino3::DistanceFactor, dn,(m, cur_pair, restraints[i]->dist, .1, st));
            factors.push_back(dn);
            dn->set_was_used(true);
        }

        std::cout << "Node size: " << factors.size() << std::endl;
        IMP::domino3::add_neighbors(factors);
        //     InteractionGraph ig = get_interaction_graph(rs, pst);
        //      SubsetGraph jt = get_junction_tree(ig);
        
        IMP_NEW(IMP::domino3::Updater, ud, (factors, "updater"));
        
        std::cout << "before" << std::endl;
        st->print_marginal();
        ud->update(iterations);
        IMP::domino3::update_state_table(factors,st);
        IMP::domino3::print_graph(factors);
        std::cout << "after" << std::endl;
        
        st->print_marginal();
        double probability_to_see_best_match = 1;
        std::vector<int> order;
        for(int i = 0; i < residues.size(); i++){
            std::vector<std::pair<double, int> > pos_order;
            IMP::domino3::Marginals * m  = st->get_marginals_by_order(i);
            for(int state = 0; state < m->get_number(); state++)
                pos_order.push_back(std::make_pair( exp(m->get_current_marginal(state)), state));
            sort(pos_order.begin(),pos_order.end(), std::greater<std::pair<double, int> >());
            probability_to_see_best_match *= pos_order[0].first;
            order.push_back(pos_order[0].second);
        }
        
        for(int i = 0 ; i < order.size();i++){
            std::cout << order[i] << " ";
        }
        std::cout << std::endl;
        
        //    std::cout << "Weighted RMSD: " << calc_weighted_rmsd(vs,ps,st) << std::endl;
    }
}



int main(int argc, char **argv) {
    IMP::base::setup_from_argv(argc, argv, "Experiment with loopy domino");
    IMP_NEW(IMP::kernel::Model, m, ());
    //    IMP::atom::Hierarchy pdb = IMP::atom:: (input, m,new IMP::atom::CAlphaPDBSelector());
    
    IMP::multifit::AnchorsData data = IMP::multifit::read_anchors_data(anchors.c_str());
    for(int i = 0; i < data.edges_.size(); i++){
        std::cout << data.edges_[i].first << " " << data.edges_[i].second << std::endl;
    }
    std::cout << "bla" << std::endl;
    //  IMP::atom::Atoms atoms
    //    = IMP::atom::get_by_type(pdb, IMP::atom::ATOM_TYPE);
    
    
    run_it(m, data.points_);
    
    return 0;
}

/**
 * Copyright 2007-2014 IMP Inventors. All rights reserved.
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
#include <algorithm>

namespace {

    
    
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
  
    
  std::string input = IMP::domino3::get_example_path("gaussian.pdb");
  IMP::base::AddStringFlag asf("input", "Input file name", &input);
  std::string output = "out.rmf";
  IMP::base::AddStringFlag oasf("output", "Output rmf name", &output);
  double dist = std::numeric_limits<double>::infinity();
  IMP::base::AddFloatFlag dasf("distance", "Maximum distance for restraints",
                              &dist);

  boost::int64_t iterations = 1;
  IMP::base::AddIntFlag aif("iterations", "Number of iterations",
                            &iterations);

    
    
    
  void run_it(IMP::kernel::Model *m, const IMP::kernel::ParticlesTemp &ps) {
    IMP::kernel::ParticleIndexes pis = IMP::kernel::get_indexes(ps);
      for( int i = 0; i < pis.size(); i++){
          std::ostringstream oss;
          oss  << pis[i];
          m->get_particle(pis[i])->set_name(oss.str());
      }

    IMP::algebra::Vector3Ds vs(pis.size());
    for (unsigned int i = 0; i< vs.size(); ++i) {
      vs[i] = IMP::core::XYZ(m, pis[i]).get_coordinates();
    }
    IMP_NEW(IMP::domino3::XYZStates, states, (m, vs));
    IMP_NEW(IMP::domino3::StatesTable, st, (m));

    // fill statetable with marginals 
    RMF::FileHandle f = RMF::create_rmf_file(output);
    for (unsigned int i = 0; i< pis.size(); ++i) {
      IMP_NEW(IMP::domino3::Marginals, marginals, (m, pis[i], vs.size()));
      marginals->set_uniform();
      st->add(pis[i], states, marginals);
    }

    st->set_rmf(f.get_root_node());

    IMP::domino3::Factors factors;
      int error_count=0;
      double avg_dist=0;
      double counter=0;
    for (unsigned int i = 0; i < pis.size(); ++i) {
      IMP::core::XYZR di(m, pis[i]);
      for (unsigned int j = 0; j < i; j++) {

        IMP::core::XYZR dj(m, pis[j]);
        di.set_coordinates(vs[i]);
        dj.set_coordinates(vs[j]);

        double cur_dist = IMP::core::get_distance(di, dj);
          avg_dist +=cur_dist;
          counter++;
        IMP::kernel::ParticleIndexPair cur_pair(pis[i], pis[j]);
        //std::cout << pis[i] << ":" << pis[j] << " Curr dist: " << cur_dist << std::endl;
        if (cur_dist < dist) {

            IMP_NEW(IMP::domino3::DistanceFactor, dn,(m, cur_pair, cur_dist, 7, st));
            factors.push_back(dn);
            dn->set_was_used(true);
            
        }
//        else {
//          IMP_NEW(IMP::domino3::ExcludedVolumeNode, dn,
//                  (m, cur_pair, st));
//          nodes.push_back(dn);
//        }
      }
    }
    //std::cout << "AVG Dist:" << avg_dist/counter;
    //std::cout << "Node error size: " << error_count << std::endl;
    //std::cout << "Node size: " << factors.size() << std::endl;
    IMP::domino3::add_neighbors(factors);
 //     InteractionGraph ig = get_interaction_graph(rs, pst);
//      SubsetGraph jt = get_junction_tree(ig);

    IMP_NEW(IMP::domino3::Updater, ud, (factors, "updater"));

    RMF::FrameID cf = f.add_frame("frame", RMF::FRAME);
    st->add_to_frame();
    //std::cout << "before" << std::endl;
    //st->show_marginal();

      
    for (unsigned int i = 0; i < iterations; ++i) {
      ud->update(100);
      cf = f.add_frame("frame", RMF::FRAME);
      st->add_to_frame();
    }
    IMP::domino3::update_state_table(factors,st);
    //IMP::domino3::print_graph(factors);
    //std::cout << "after" << std::endl;

    //st->show_marginal();
    //std::cout << "Weighted RMSD: " << calc_weighted_rmsd(vs,ps,st) << std::endl;

    double  diagonal[] =  {1,	1,  1	,1 , 1	,1, 1, 0.999085, 0.845624, 0.846366}; 
      for( int i = 0; i < pis.size(); i++){
          IMP::domino3::Marginals * m0 = st->get_marginals_by_order(i);
          //std::cout <<  diagonal[i] << " " << exp(m0->get_current_marginal(i)) << std::endl;
          if( std::abs(diagonal[i] - exp(m0->get_current_marginal(i))) >0.1){
              std::cout << "Error in Test" << std::endl;
              exit(-1);
          }
      }
      //std::cout << "Test successful" <<std::endl;
  }
}



int main(int argc, char **argv) {
  IMP::base::setup_from_argv(argc, argv, "Experiment with loopy domino");
  IMP_NEW(IMP::kernel::Model, m, ());
    IMP::atom::Hierarchy pdb = IMP::atom::read_pdb(input, m,new IMP::atom::CAlphaPDBSelector());
  IMP::atom::Atoms atoms
    = IMP::atom::get_by_type(pdb, IMP::atom::ATOM_TYPE);


    run_it(m, IMP::kernel::ParticlesTemp(atoms.begin(), atoms.end()));

  return 0;
}

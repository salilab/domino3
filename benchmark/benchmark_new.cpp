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
#include <IMP/domino3/Node.h>
#include <IMP/domino3/DistanceNode.h>
#include <IMP/domino3/ExcludedVolumeNode.h>
#include <IMP/domino3/Updater.h>
#include <RMF/FileHandle.h>
#include <IMP/domino3/XYZStates.h>
#include <IMP/core/XYZ.h>
#include <algorithm>

namespace {

  std::string input = IMP::atom::get_example_path("1d3d-protein.pdb");
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

    IMP::domino3::Nodes nodes;
    for (unsigned int i = 0; i < pis.size(); ++i) {
      IMP::core::XYZR di(m, pis[i]);
      for (unsigned int j = 0; j < i; j++) {
//        if(rand() % 100 > 70)
//            break;
        IMP::core::XYZR dj(m, pis[j]);
        di.set_coordinates(vs[i]);
        dj.set_coordinates(vs[j]);

        double cur_dist = IMP::core::get_distance(di, dj);
        IMP::kernel::ParticleIndexPair cur_pair(pis[i], pis[j]);
        if (cur_dist < dist) {
          IMP_NEW(IMP::domino3::DistanceNode, dn,
                  (m, cur_pair, cur_dist, .1, st));
          nodes.push_back(dn);
        }
//        else {
//          IMP_NEW(IMP::domino3::ExcludedVolumeNode, dn,
//                  (m, cur_pair, st));
//          nodes.push_back(dn);
//        }
      }
    }
      std::cout << "Node size: " << nodes.size() << std::endl;
    IMP::domino3::add_neighbors(nodes);
//      InteractionGraph ig = get_interaction_graph(rs, pst);
//      SubsetGraph jt = get_junction_tree(ig);

    IMP_NEW(IMP::domino3::Updater, ud, (nodes, "updater"));

    RMF::FrameHandle cf = f.get_current_frame().add_child("frame", RMF::FRAME);
    st->add_to_frame();
    std::cout << "before" << std::endl;
    st->print_marginal();

      
    for (unsigned int i = 0; i < iterations; ++i) {
      ud->update(6);
      cf = cf.add_child("frame", RMF::FRAME);
      st->add_to_frame();
    }
    IMP::domino3::update_state_table(nodes,st);
      IMP::domino3::print_graph(nodes);
    std::cout << "after" << std::endl;

    st->print_marginal();
  }
}


int main(int argc, char **argv) {
  IMP::base::setup_from_argv(argc, argv, "Experiment with loopy domino");
  IMP_NEW(IMP::kernel::Model, m, ());
    IMP::atom::Hierarchy pdb = IMP::atom::read_pdb(input, m,new IMP::atom::CAlphaPDBSelector());
  IMP::atom::Atoms atoms
    = IMP::atom::get_by_type(pdb, IMP::atom::ATOM_TYPE);

//    for( int i = 0; i < pis.size(); i++){
//        std::ostringstream oss;
//        oss  << pis[i];
//        m->get_particle(pis[i])->set_name(oss.str());
//    }

    run_it(m, IMP::kernel::ParticlesTemp(atoms.begin(), atoms.end()));

  return 0;
}

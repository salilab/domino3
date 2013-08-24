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
#include <IMP/domino3/ErrorDistanceNode.h>
#include <IMP/domino3/SeaNode.h>

#include <IMP/domino3/ExcludedVolumeNode.h>
#include <IMP/domino3/Updater.h>
#include <RMF/FileHandle.h>
#include <IMP/domino3/IndexStates.h>
#include <IMP/core/XYZ.h>
#include <algorithm>

namespace {
  typedef std::pair<IMP::kernel::Particle *,IMP::kernel::Particle *> ParticlePair;
  typedef IMP::base::map<ParticlePair, double > ParticleScoreLookup;
  typedef IMP::base::map<std::string, IMP::kernel::Particle *>  StringToParticle;
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

    
    void read_in_pair_file(IMP::kernel::Model *m,
                           std::string path,
                           ParticleScoreLookup &scores,
                           StringToParticle &id_to_particle ){
        std::ifstream infile(path.c_str());
        std::string line;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            std::string a, b;
            double score;
            iss >> a >> b >> score;
            if(!id_to_particle[a])
                id_to_particle[a] = new IMP::kernel::Particle(m);
            if(!id_to_particle[b])
                id_to_particle[b] = new IMP::kernel::Particle(m);
            scores[std::make_pair(id_to_particle[a],id_to_particle[b])]=score;
            std::cout << a << " " << b << " " << score << std::endl;
        }

    }
    
  void run_it(IMP::kernel::Model *m) {
      StringToParticle id_to_particle;
      ParticleScoreLookup scores;
      read_in_pair_file(m,
                        std::string("/Users/aluucard/Documents/workspace/imp/imp/modules/domino3/examples/sea.txt"),
                        scores,
                        id_to_particle);
      std::cout << id_to_particle.size() << std::endl;
      IMP::kernel::Particles ps;

      for (StringToParticle::iterator i = id_to_particle.begin(); i != id_to_particle.end(); ++i){
          IMP::kernel::Particle * p = i->second;
          ps.push_back(p);
      }
      IMP_NEW(IMP::domino3::StatesTable, st, (m));
      IMP_NEW(IMP::domino3::IndexStates, states, (m, ps));

      for (StringToParticle::iterator i = id_to_particle.begin(); i != id_to_particle.end(); ++i){
          IMP::kernel::Particle * p = i->second;
          IMP_NEW(IMP::domino3::Marginals, marginals, (m, p->get_index(), id_to_particle.size()));
          marginals->make_next_zero();
          st->add(p->get_index(), states, marginals);
      }
      
      for (ParticleScoreLookup::iterator i = scores.begin(); i != scores.end(); ++i)
      {
          ParticlePair particlePair=i->first;
          double score = i->second;
          IMP::kernel::Particle * p1 = particlePair.first;
          IMP::kernel::Particle * p2 = particlePair.second;
          IMP::kernel::ParticleIndex p1_index     = p1->get_index();
          IMP::kernel::ParticleIndex p2_index     = p2->get_index();
          IMP::domino3::Marginals * m1 = st->get_marginals(p1_index);
          double converted_score=m1->convert_to_space(score);
          m1->add_to_next_marginal(p2_index.get_index(), converted_score);
          std::cout << i->second << std::endl;
      }
      for (StringToParticle::iterator i = id_to_particle.begin(); i != id_to_particle.end(); ++i){
          IMP::kernel::Particle * p = i->second;
          IMP::domino3::Marginals * m1 = st->get_marginals(p->get_index());
          m1->set_current_from_next();
      }
      st->print_marginal();
      IMP::domino3::Nodes nodes;
      for (ParticleScoreLookup::iterator i = scores.begin(); i != scores.end(); ++i)
      {
          ParticlePair particlePair=i->first;
          IMP::kernel::Particle * p1 = particlePair.first;
          IMP::kernel::Particle * p2 = particlePair.second;
          IMP::kernel::ParticleIndex p1_index     = p1->get_index();
          IMP::kernel::ParticleIndex p2_index     = p2->get_index();
          IMP::kernel::ParticleIndexPair cur_pair(p1_index, p2_index);

          IMP_NEW(IMP::domino3::SeaNode, sn,(m, cur_pair, st));
          nodes.push_back(sn);
          sn->set_was_used(true);
      }
//      IMP_NEW(IMP::domino3::SeaNode, sn,(m, cur_pair, st));
//      nodes.push_back(sn);
//      sn->set_was_used(true);
//      std::cout << "Node error size: " << error_count << std::endl;
      std::cout << "Node size: " << nodes.size() << std::endl;
    IMP::domino3::add_neighbors(nodes);
//
    IMP_NEW(IMP::domino3::Updater, ud, (nodes, "updater"));
//
//    std::cout << "before" << std::endl;
//    st->print_marginal();
//
//      
    for (unsigned int i = 0; i < iterations; ++i) {
      ud->update(100);
    }
    IMP::domino3::update_state_table(nodes,st);
    IMP::domino3::print_graph(nodes);

    st->print_marginal();
//    std::cout << "after" << std::endl;
//
//
  }
}



int main(int argc, char **argv) {
  IMP::base::setup_from_argv(argc, argv, "Experiment with loopy domino");
  IMP_NEW(IMP::kernel::Model, m, ());
    
    std::vector<IMP::kernel::Particle *> enzyms;


    run_it(m);

  return 0;
}

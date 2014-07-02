/**
 *  \file IMP/domino3/Factor.cpp
 *
 *  Copyright 2007-2014 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino3/Factor.h>
#include <boost/unordered_map.hpp>

IMPDOMINO3_BEGIN_NAMESPACE

Factor::Factor(Model *m,
           const ParticleIndexes &pis,
           StatesTable *pst,
           std::string name):
  ModelObject(m, name), pis_(pis) {
  std::sort(pis_.begin(), pis_.end());
  mine_.resize(pis_.size());
  inputs_.resize(pis_.size());
  for (unsigned int i = 0; i < pis_.size(); ++i) {
//    mine_[i] = new Marginals(m, pis_[i],
//                             pst->get_states(pis_[i])->get_number());
      
      
      IMP_NEW(IMP::domino3::Marginals, marginals, (m, pis_[i], pst->get_states(pis_[i])->get_number()));
      mine_[i] = marginals;
    mine_[i]->set_name(mine_[i]->get_name()+" "+this->get_name());
    mine_[i]->set_uniform();
    // give marginals values from state list (but its also uniform?!?
    mine_[i]->merge_probabilities_from_list(MarginalsList(1,
                             pst->get_marginals(pis_[i])));
  }
}

void Factor::update() {
  // propergate my list to neighbours. Inputs is from neighbor (MPI call in futur)
  for (unsigned int i = 0; i < pis_.size(); ++i) {
      for(int j = 0; j < inputs_[i].size(); j++) {
          inputs_[i][j]->check_current_normalized();
      }
  }
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i]->merge_probabilities_from_list(inputs_[i]);
    mine_[i]->make_next_zero();
  }
    
  do_update();
  // propergate my list to neighbours 
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i]->set_current_from_next();
    mine_[i]->check_current_normalized();
  }
}

void Factor::set_matching_inputs(Factor *n){
    ParticleIndexes pis = n->get_particle_indexes();
    for (unsigned int i = 0; i < pis.size(); ++i) {
        // just add variables that are really match
        kernel::ParticleIndexes::const_iterator it = std::find(pis_.begin(),
                                                               pis_.end(),
                                                               pis[i]);
        if (it != pis_.end()) {
            unsigned int offset = it - pis_.begin();
            // not so good for parallel approaches
            inputs_[offset].push_back(n->get_marginals()[i]);
            //      n->get_marginals()[i]->show_marginals();
        }
    }
}

void Factor::add_neighbor(Factor *n) {
  IMP_INTERNAL_CHECK(n != this,"Omg its wrong");
  neighbors_.push_back(n);
}

void add_neighbors_by_factor_edges(const FactorEdges &factor_edges){
	for (unsigned int i = 0; i < factor_edges.size(); ++i) {
		const FactorEdge * edge = &factor_edges[i];
        edge->from_->add_neighbor(edge->to_);
        edge->from_->set_matching_inputs(edge->to_);
        edge->to_->add_neighbor(edge->from_);
        edge->to_->set_matching_inputs(edge->from_);
	}
}



FactorGraph get_node_graph(const FactorsTemp &nodes) {
  FactorGraph ret(nodes.size());
  FactorGraphVertexName vm = boost::get(boost::vertex_name, ret);
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    vm[i] = nodes[i];
    // create self edge too
    for (unsigned int j = 0; j <= i; ++j) {
      kernel::ParticleIndexes intersection;
      std::set_intersection(nodes[i]->get_particle_indexes().begin(),
                            nodes[i]->get_particle_indexes().end(),
                            nodes[j]->get_particle_indexes().begin(),
                            nodes[j]->get_particle_indexes().end(),
                            std::back_inserter(intersection));
      if (!intersection.empty()) {
        boost::add_edge(i, j, intersection, ret);
      }
    }
  }
  return ret;
}


void add_neighbors(const FactorsTemp &nodes) {
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    kernel::ParticleIndexes pisi = nodes[i]->get_particle_indexes();
    for (unsigned int j = 0; j < i; ++j) {
      kernel::ParticleIndexes pisj = nodes[j]->get_particle_indexes();
      kernel::ParticleIndexes intersection;
      std::set_intersection(pisi.begin(),
                            pisi.end(),
                            pisj.begin(),
                            pisj.end(),
                            std::back_inserter(intersection));
      if (!intersection.empty()) {
        nodes[i]->add_neighbor(nodes[j]);
        nodes[i]->set_matching_inputs(nodes[j]);
        nodes[j]->add_neighbor(nodes[i]);
        nodes[j]->set_matching_inputs(nodes[i]);
      }
    }
  }
}


void update_state_table(const FactorsTemp &nodes,const StatesTable *pst) {
    boost::unordered_map<ParticleIndex,MarginalsList> map_to_merge;
    for (unsigned int i = 0; i < nodes.size(); ++i) {
        MarginalsList node_marginals =  nodes[i]->get_marginals();
        for(MarginalsList::size_type i = 0; i != node_marginals.size(); i++) {
            ParticleIndex particle_index = node_marginals[i]->get_particle_index();
            map_to_merge[particle_index].push_back(node_marginals[i]);
        }
    }
    for(boost::unordered_map<ParticleIndex,MarginalsList>::iterator iter = map_to_merge.begin();
        iter != map_to_merge.end(); ++iter)
    {
        ParticleIndex k =  iter->first;
        MarginalsList node_list =  iter->second;
        Marginals * pst_marginal = pst->get_marginals(k);
        pst_marginal->merge_probabilities_from_list(node_list);
    }
}
void print_graph(const FactorsTemp &nodes){

    for (unsigned int i = 0; i < nodes.size(); ++i) {
        kernel::ParticleIndexes particles=nodes[i]->get_particle_indexes();
        std::cout << particles << " -> ";
        FactorsTemp neighbors = nodes[i]->get_neighbors();
        for(int j = 0; j < neighbors.size(); j++){
            std::cout << neighbors[j]->get_particle_indexes() << " ";
        }
        std::cout << std::endl;
    }
}


IMPDOMINO3_END_NAMESPACE

#include <IMP/domino3/Node.h>

IMPDOMINO3_BEGIN_NAMESPACE

Node::Node(Model *m,
                 const ParticleIndexes &pis,
                 domino::ParticleStatesTable *pst,
                 std::string name):
  ModelObject(m, name), pis_(pis) {
  std::sort(pis_.begin(), pis_.end());
  mine_.resize(pis_.size());
  inputs_.resize(pis_.size());
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i] = new Marginals(m, pis_[i],
                          pst->get_number_of_states(m->get_particle(pis_[i])));
  }
}

void Node::add_input_marginal(ParticleIndex pi, Marginal *m) {
  inputs[pi].push_back(m);
}

void Node::update() {
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i]->update_current_from_list(inputs_[i]);
  }

  do_update();

  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i]->update_current_next();
  }
}

IMPDOMINO3_END_NAMESPACE

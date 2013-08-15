#include <IMP/domino3/DistanceNode.h>
#include <IMP/core/XYZR.h>

IMPDOMINO3_BEGIN_NAMESPACE

namespace {
  std::string get_distance_name(const kernel::ParticleIndexPair &pis) {
    std::ostringstream oss;
    oss << "Distance-" << pis[0] << "-" << pis[1];
    return oss.str();
  }
}

DistanceNode::DistanceNode(kernel::Model *m,
               const kernel::ParticleIndexPair &pis,
               double distance, double allowed_error,
               StatesTable *pst):
  Node(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_distance_name(pis)), distance_(distance),
  allowed_error_(allowed_error) {
  // precompile allowed states, use some sort of location structure later
  States *ps0 = pst->get_states(pis[0]), *ps1 = pst->get_states(pis[1]);
  for (unsigned int i = 0; i < ps0->get_number(); ++i) {
    ps0->load(i, pis[0]);
    for (unsigned int j = 0; j < ps1->get_number(); ++j) {
      ps1->load(j, pis[1]);
      double d = IMP::core::get_distance(IMP::core::XYZR(m, pis[0]),
                                         IMP::core::XYZR(m, pis[1]));
      if (std::abs(d- distance) < allowed_error) {
        allowed_states_.push_back(IntPair(i, j));
      }
    }
  }
  IMP_USAGE_CHECK(allowed_states_.size() > 0, "No allowed state could be found");
}

void DistanceNode::do_update() {
  Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
  for (unsigned int i = 0; i < allowed_states_.size(); ++i) {
    double cur = m0->mult_two_marginals(m0->get_current_marginal(allowed_states_[i].first), m1->get_current_marginal(allowed_states_[i].second));
    m0->add_to_next_marginal(allowed_states_[i].first, cur);
    m1->add_to_next_marginal(allowed_states_[i].second, cur);
  }
}

IMPDOMINO3_END_NAMESPACE

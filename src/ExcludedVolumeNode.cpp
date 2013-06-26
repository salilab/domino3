#include <IMP/domino3/ExcludedVolumeNode.h>

IMPDOMINO3_BEGIN_NAMESPACE

namespace {
  std::string get_ev_name(const kernel::ParticleIndexPair &pis) {
    std::ostringstream oss;
    oss << "EV-" << pis[0] << "-" << pis[1];
    return oss.str();
  }
}

ExcludedVolumeNode::ExcludedVolumeNode(Model *m,
               const kernel::ParticleIndexPair &pis,
               domino::ParticleStatesTable *pst):
  Node(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_distance_name(pis)) {
}

void ExcludedVolumeNode::do_update() {
  Marginals *m0 = get_marginals()[0], *m1 = get_margins()[1];
  for (unsigned int i = 0; i < m0->size(); ++i) {
    double cur = m0->get_marginal(i) * m1->get_margin(i);
    m0->add_to_marginal(i, 1.0 - cur);
    m1->add_to_marginal(i, 1.0 - cur);
  }
}

IMPDOMINO3_END_NAMESPACE

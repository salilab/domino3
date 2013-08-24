#include <IMP/domino3/SeaNode.h>
#include <IMP/core/XYZR.h>

IMPDOMINO3_BEGIN_NAMESPACE
namespace {
    std::string get_sea_name(const kernel::ParticleIndexPair &pis) {
        std::ostringstream oss;
        oss << "Sea-" << pis[0] << "-" << pis[1];
        return oss.str();
    }
}

SeaNode::SeaNode(kernel::Model *m,const kernel::ParticleIndexPair &pis,
               StatesTable *pst):
  Node(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_sea_name(pis)){
  // precompile allowed states, use some sort of location structure later
  this->pis = pis;
  this->pst = pst;
}


void SeaNode::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
    States *ps0 = pst->get_states(pis[0]), *ps1 = pst->get_states(pis[1]);
    for (unsigned int i = 0; i < ps0->get_number(); ++i) {
        for (unsigned int j = 0; j < ps1->get_number(); ++j) {
               double cur = m0->mult_two_marginals(m0->get_current_marginal(i),
                                                    m1->get_current_marginal(j));
               m0->add_to_next_marginal(i, cur);
               m1->add_to_next_marginal(j, cur);
        }
    }
}

IMPDOMINO3_END_NAMESPACE

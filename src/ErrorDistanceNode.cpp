#include <IMP/domino3/ErrorDistanceNode.h>
#include <IMP/core/XYZR.h>

IMPDOMINO3_BEGIN_NAMESPACE
namespace {
    std::string get_error_distance_name(const kernel::ParticleIndexPair &pis) {
        std::ostringstream oss;
        oss << "ErrorDistance-" << pis[0] << "-" << pis[1];
        return oss.str();
    }
}

ErrorDistanceNode::ErrorDistanceNode(kernel::Model *m,
               const kernel::ParticleIndexPair &pis,
               double distance, double allowed_error,
               StatesTable *pst):
  Node(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_error_distance_name(pis)), distance_(distance),
  allowed_error_(allowed_error) {
  // precompile allowed states, use some sort of location structure later
  this->pis = pis;
  this->pst = pst;
      this->distance = distance;
      this->m = m;
      this->allowed_error = allowed_error;
}

double ErrorDistanceNode::distance_to_probability(double x){
    x=x/this->allowed_error;
    return 1-std::max(exp(-pow(x,2)*M_PI),0.00000001);
}


void ErrorDistanceNode::do_update() {
    
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];

    States *ps0 = pst->get_states(pis[0]), *ps1 = pst->get_states(pis[1]);
    for (unsigned int i = 0; i < ps0->get_number(); ++i) {
        ps0->load(i, pis[0]);
        for (unsigned int j = 0; j < ps1->get_number(); ++j) {
            ps1->load(j, pis[1]);
            double d = IMP::core::get_distance(IMP::core::XYZR(m, pis[0]),
                                               IMP::core::XYZR(m, pis[1]));
            
                double cur = m0->mult_two_marginals(m0->get_current_marginal(i),
                                                    m1->get_current_marginal(j));
                cur+=log(distance_to_probability(d-this->distance));
                m0->add_to_next_marginal(i, cur);
                m1->add_to_next_marginal(j, cur);
            
        }
    }
    
}

IMPDOMINO3_END_NAMESPACE
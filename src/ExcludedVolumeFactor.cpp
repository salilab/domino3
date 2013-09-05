#include <IMP/domino3/ExcludedVolumeFactor.h>

IMPDOMINO3_BEGIN_NAMESPACE

namespace {
  std::string get_ev_name(const kernel::ParticleIndexPair &pis) {
    std::ostringstream oss;
    oss << "EV-" << pis[0] << "-" << pis[1];
    return oss.str();
  }
}

ExcludedVolumeFactor::ExcludedVolumeFactor(Model *m,
               const kernel::ParticleIndexPair &pis,
               StatesTable *pst):
  Factor(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_ev_name(pis)) {
}

void ExcludedVolumeFactor::do_update() {
  Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
    const int size_i=m0->get_number();
    const int size_j=m1->get_number();
    for (unsigned int i = 0; i < size_i; ++i) {
        const double m0_i_current_marginal = m0->get_current_marginal(i);
        for (unsigned int j = 0; j < size_j; ++j) {
            const double m1_j_current_marginal = m1->get_current_marginal(j);
            double cur = 0;
            if(i==j){
                cur = log(0.00001);
            }else{
                cur = log(1.0);
            }
            m0->add_to_next_marginal(i,cur + m1_j_current_marginal);
            m1->add_to_next_marginal(j,cur + m0_i_current_marginal);
        }
    }
}

IMPDOMINO3_END_NAMESPACE

#include <IMP/domino3/Probability2DFactor.h>
#include <IMP/core/XYZR.h>

IMPDOMINO3_BEGIN_NAMESPACE
namespace {
    std::string get_prob_name(const kernel::ParticleIndexPair &pis) {
        std::ostringstream oss;
        oss << "ProbabilityFactor-" << pis[0] << "-" << pis[1];
        return oss.str();
    }
}

Probability2DFactor::Probability2DFactor(kernel::Model *m,const kernel::ParticleIndexPair &pis,
                                     StatesTable *pst, boost::shared_array<FP> &log_probability):
  Factor(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
         pst, get_prob_name(pis)), log_probability_(log_probability),pis_(pis),pst_(pst) {
}


void Probability2DFactor::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
    States *ps0 = pst_->get_states(pis_[0]), *ps1 = pst_->get_states(pis_[1]);
    const int size_i = ps0->get_number();
    const int size_j = ps1->get_number();
    for (unsigned int i = 0; i < size_i; ++i) {
        const FP m0_i_current_marginal = m0->get_current_marginal(i);
        for (unsigned int j = 0; j < size_j; ++j) {
            const FP log_prob_i_j=this->log_probability_[i*size_j+j];
            const FP m1_j_current_marginal = m1->get_current_marginal(j);
            m0->add_to_next_marginal(i,log_prob_i_j + m1_j_current_marginal);
            m1->add_to_next_marginal(j,log_prob_i_j + m0_i_current_marginal);
        }
    }
}

IMPDOMINO3_END_NAMESPACE

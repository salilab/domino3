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
                                     StatesTable *pst, boost::shared_array<double> &probability):
  Factor(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
         pst, get_prob_name(pis)), probability_(probability),pis_(pis),pst_(pst) {
}


void Probability2DFactor::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
    States *ps0 = pst_->get_states(pis_[0]), *ps1 = pst_->get_states(pis_[1]);
    for (unsigned int i = 0; i < ps0->get_number(); ++i) {
        for (unsigned int j = 0; j < ps1->get_number(); ++j) {
               m0->add_to_next_marginal(i, log(this->probability_[i*ps1->get_number()+j])+m1->get_current_marginal(j));
               m1->add_to_next_marginal(j, log(this->probability_[i*ps1->get_number()+j])+m0->get_current_marginal(i));
        }
    }
}

IMPDOMINO3_END_NAMESPACE

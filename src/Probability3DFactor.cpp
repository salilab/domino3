#include <IMP/domino3/Probability3DFactor.h>
#include <IMP/core/XYZR.h>

IMPDOMINO3_BEGIN_NAMESPACE
namespace {
    std::string get_prob3d_name(const kernel::ParticleIndexTriplet &pis) {
        std::ostringstream oss;
        oss << "Probability3dFactor-" << pis[0] << "-" << pis[1]<< "-" << pis[2] ;
        return oss.str();
    }
}

Probability3DFactor::Probability3DFactor(kernel::Model *m,const kernel::ParticleIndexTriplet &pis,
                                         StatesTable *pst, boost::shared_array<double> &probability):
Factor(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_prob3d_name(pis)), probability_(probability),pis_(pis),pst_(pst) {
}


void Probability3DFactor::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1], *m2 = get_marginals()[2];
    States *ps0 = pst_->get_states(pis_[0]), *ps1 = pst_->get_states(pis_[1]), *ps2 = pst_->get_states(pis_[2]);
    for (unsigned int i = 0; i < ps0->get_number(); ++i) {
        for (unsigned int j = 0; j < ps1->get_number(); ++j) {
            for (unsigned int z = 0; z < ps2->get_number(); ++z) {
                int probability_index = i*ps1->get_number() * ps2->get_number()+j* ps2->get_number() + z;
                m0->add_to_next_marginal(i, log(this->probability_[probability_index])+m1->get_current_marginal(j)+m2->get_current_marginal(z));
                m1->add_to_next_marginal(j, log(this->probability_[probability_index])+m0->get_current_marginal(i)+m2->get_current_marginal(z));
                m2->add_to_next_marginal(z, log(this->probability_[probability_index])+m0->get_current_marginal(i)+m1->get_current_marginal(j));
            }
        }
    }
}

IMPDOMINO3_END_NAMESPACE

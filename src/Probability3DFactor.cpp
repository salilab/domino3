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
                                         StatesTable *pst, boost::shared_array<double> &log_probability):
Factor(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_prob3d_name(pis)), log_probability_(log_probability),pis_(pis),pst_(pst) {
    
}


void Probability3DFactor::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1], *m2 = get_marginals()[2];
    States *ps0 = pst_->get_states(pis_[0]), *ps1 = pst_->get_states(pis_[1]), *ps2 = pst_->get_states(pis_[2]);
    const int size_i=ps0->get_number();
    const int size_j=ps1->get_number();
    const int size_z=ps2->get_number();

    for (unsigned int i = 0; i < size_i; ++i) {
        const double m0_i_current_marginal = m0->get_current_marginal(i);
        for (unsigned int j = 0; j < size_j; ++j) {
            const double m1_j_current_marginal = m1->get_current_marginal(j);
            for (unsigned int z = 0; z < size_z; ++z) {
                const double m2_z_current_marginal = m2->get_current_marginal(z);
                const int probability_index = i * size_j * size_z +j * size_z + z;
                const double log_prob_i_j_z = this->log_probability_[probability_index];
                m0->add_to_next_marginal(i, log_prob_i_j_z+m1_j_current_marginal+m2_z_current_marginal);
                m1->add_to_next_marginal(j, log_prob_i_j_z+m0_i_current_marginal+m2_z_current_marginal);
                m2->add_to_next_marginal(z, log_prob_i_j_z+m0_i_current_marginal+m1_j_current_marginal);
            }
        }
    }
}

IMPDOMINO3_END_NAMESPACE

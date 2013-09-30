#include <IMP/domino3/Probability3DFactor.h>
#include <IMP/core/XYZR.h>
#include <IMP/domino3/LogMathFunctions.h>

IMPDOMINO3_BEGIN_NAMESPACE
namespace {
    std::string get_prob3d_name(const kernel::ParticleIndexTriplet &pis) {
        std::ostringstream oss;
        oss << "Probability3dFactor-" << pis[0] << "-" << pis[1]<< "-" << pis[2] ;
        return oss.str();
    }
}

Probability3DFactor::Probability3DFactor(kernel::Model *m,const kernel::ParticleIndexTriplet &pis,
                                         StatesTable *pst, Probability3D * log_probability):
Factor(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_prob3d_name(pis)), log_probability(log_probability),pis_(pis),pst_(pst) {
    
    States *ps2 = pst_->get_states(pis_[2]);

}

//
//void Probability3DFactor::do_update() {
//    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1], *m2 = get_marginals()[2];
//    States *ps0 = pst_->get_states(pis_[0]), *ps1 = pst_->get_states(pis_[1]), *ps2 = pst_->get_states(pis_[2]);
//    const int size_i=ps0->get_number();
//    const int size_j=ps1->get_number();
//    const int size_z=ps2->get_number();
//
//    for (unsigned int i = 0; i < size_i; ++i) {
//        const FP m0_i_current_marginal = m0->get_current_marginal(i);
//        for (unsigned int j = 0; j < size_j; ++j) {
//            const FP m1_j_current_marginal = m1->get_current_marginal(j);
//            for (unsigned int z = 0; z < size_z; ++z) {
//                const FP m2_z_current_marginal = m2->get_current_marginal(z);
//                const FP log_prob_i_j_z = log_probability->get(i,j,z);
//                m0->add_to_next_marginal(i, log_prob_i_j_z+m1_j_current_marginal+m2_z_current_marginal);
//                //next[j] = 4*(prob + i_prob + z prob)
//                m1->add_to_next_marginal(j, log_prob_i_j_z+m0_i_current_marginal+m2_z_current_marginal);
//                m2->add_to_next_marginal(z, log_prob_i_j_z+m0_i_current_marginal+m1_j_current_marginal);
//            }
//
//        }
//
//    }
//}

//

void Probability3DFactor::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1], *m2 = get_marginals()[2];
    States *ps0 = pst_->get_states(pis_[0]), *ps1 = pst_->get_states(pis_[1]), *ps2 = pst_->get_states(pis_[2]);
    const int size_i=ps0->get_number();
    const int size_j=ps1->get_number();
    const int size_z=ps2->get_number();
    // memory for z marginal
    FP __attribute__((aligned(16))) z_cur_m2_float[size_z+4];
    std::fill(&z_cur_m2_float[0],  &z_cur_m2_float[size_z+4],LogMathFunctions::convert_to_space(0.0f));
    FP * m2_z_current_marginal = m2->get_current_marginal();
    std::copy (m2_z_current_marginal, m2_z_current_marginal+size_z, z_cur_m2_float );
    __m128 * z_cur_m2 =(__m128 * ) &z_cur_m2_float[0];
    
    FP __attribute__((aligned(16))) z_next_m2_float[size_z+4];
    std::fill(&z_next_m2_float[0],  &z_next_m2_float[size_z+4],LogMathFunctions::convert_to_space(0.0f));
    FP * m2_z_next_marginal = m2->get_next_marginal();
    std::copy (m2_z_next_marginal, m2_z_next_marginal+size_z, z_next_m2_float );
    __m128 * z_next_m2 =(__m128 * ) &z_next_m2_float[0];
    
    FP init_value=LogMathFunctions::convert_to_space(0.0f);
    for (unsigned int i = 0; i < size_i; ++i) {
        const FP m0_i_current_marginal = m0->get_current_marginal(i);
        __m128 i_cur_m0 = _mm_set1_ps(m0_i_current_marginal);
        __m128 i_next_m0 = _mm_set1_ps(init_value);

        for (unsigned int j = 0; j < size_j; ++j) {
            const FP m1_j_current_marginal = m1->get_current_marginal(j);
            __m128 j_cur_m1 = _mm_set1_ps(m1_j_current_marginal);
            __m128 j_next_m1 = _mm_set1_ps(init_value);
            __m128 * prob_row = (__m128*) this->log_probability->get_ptr(i, j, 0);
            for (unsigned int z = 0; z < (size_z/4)+1; z++) {
                
                __m128 tmp  = LogMathFunctions::mult_sse(prob_row[z],j_cur_m1);
                tmp         = LogMathFunctions::mult_sse(tmp,z_cur_m2[z]);
                i_next_m0   = LogMathFunctions::add_sse(i_next_m0,tmp);

                tmp         = LogMathFunctions::mult_sse(prob_row[z],i_cur_m0);
                tmp         = LogMathFunctions::mult_sse(tmp,z_cur_m2[z]);
                j_next_m1   = LogMathFunctions::add_sse(j_next_m1,tmp);
                float * p = (float *)&j_next_m1;
                
                tmp         = LogMathFunctions::mult_sse(prob_row[z],i_cur_m0);
                tmp         = LogMathFunctions::mult_sse(tmp,j_cur_m1);
                z_next_m2[z]= LogMathFunctions::add_sse(z_next_m2[z],tmp);

            }
            // update j
            m1->add_to_next_marginal(j, (FP) LogMathFunctions::get_first (j_next_m1) );
            m1->add_to_next_marginal(j, (FP) LogMathFunctions::get_second(j_next_m1) );
            m1->add_to_next_marginal(j, (FP) LogMathFunctions::get_third (j_next_m1) );
            m1->add_to_next_marginal(j, (FP) LogMathFunctions::get_fourth(j_next_m1) );            
        }
        // update i
        m0->add_to_next_marginal(i, (FP) LogMathFunctions::get_first (i_next_m0 ));
        m0->add_to_next_marginal(i, (FP) LogMathFunctions::get_second(i_next_m0 ));
        m0->add_to_next_marginal(i, (FP) LogMathFunctions::get_third (i_next_m0 ));
        m0->add_to_next_marginal(i, (FP) LogMathFunctions::get_fourth(i_next_m0 ));
    }
    // update z
    std::copy (z_next_m2_float, z_next_m2_float+size_z, m2_z_next_marginal );
}


IMPDOMINO3_END_NAMESPACE

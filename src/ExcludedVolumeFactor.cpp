/**
 *  \file IMP/domino3/ExcludedVolumeFactor.cpp
 *
 *  Copyright 2007-2015 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino3/ExcludedVolumeFactor.h>
#include <IMP/domino3/LogMathFunctions.h>
IMPDOMINO3_BEGIN_NAMESPACE

namespace {
    std::string get_ev_name(const ParticleIndexPair &pis) {
        std::ostringstream oss;
        oss << "EV-" << pis[0] << "-" << pis[1];
        return oss.str();
    }
}



ExcludedVolumeFactor::ExcludedVolumeFactor(Model *m,
                                           const ParticleIndexPair &pis,
                                           StatesTable *pst):
Factor(m, ParticleIndexes(pis.begin(), pis.end()),
       pst, get_ev_name(pis)) {
}

//void ExcludedVolumeFactor::do_update() {
//  Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
//    const int size_i=m0->get_number();
//    const int size_j=m1->get_number();
//    for (unsigned int i = 0; i < size_i; ++i) {
//        const FP m0_i_current_marginal = m0->get_current_marginal(i);
//        for (unsigned int j = 0; j < size_j; ++j) {
//            const FP m1_j_current_marginal = m1->get_current_marginal(j);
//            FP cur = 0;
//            if(i==j){
//                cur = LogMathFunctions::convert_to_space(0.00001);
//            }else{
//                cur = LogMathFunctions::convert_to_space(1.0);
//            }
//            m0->add_to_next_marginal(i,cur + m1_j_current_marginal);
//            m1->add_to_next_marginal(j,cur + m0_i_current_marginal);
//        }
//    }
//}



void ExcludedVolumeFactor::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
    const int size_i=m0->get_number();
    const int size_j=m1->get_number();
    
    // init vector for j curr
    FP __attribute__((aligned(16))) j_cur_m1_float[size_j+4];
    std::fill(&j_cur_m1_float[0],  &j_cur_m1_float[size_j+4],LogMathFunctions::convert_to_space(0.0f));
    FP * m1_j_current_marginal = m1->get_current_marginal();
    std::copy (m1_j_current_marginal, m1_j_current_marginal+size_j, j_cur_m1_float );
    __m128 * j_cur_m1 =(__m128 * ) &j_cur_m1_float[0];
    // init vector for j next
    FP __attribute__((aligned(16))) j_next_m1_float[size_j+4];
    std::fill(&j_next_m1_float[0],  &j_next_m1_float[size_j+4],LogMathFunctions::convert_to_space(0.0f));
    FP * m1_j_next_marginal = m1->get_next_marginal();
    std::copy (m1_j_next_marginal, m1_j_next_marginal+size_j, j_next_m1_float );
    __m128 * j_next_m1 =(__m128 * ) &j_next_m1_float[0];
    // set prob. vector
    FP __attribute__((aligned(16))) prob_float[size_j+4];
    std::fill(&prob_float[0],  &prob_float[size_j+4],LogMathFunctions::convert_to_space(1.0f));
    __m128 * prob_vec =(__m128 * ) &prob_float[0];
    
    FP init_value=LogMathFunctions::convert_to_space(0.0f);
    for (unsigned int i = 0; i < size_i; ++i) {
        const FP m0_i_current_marginal = m0->get_current_marginal(i);
        __m128 i_cur_m0 = _mm_set1_ps(m0_i_current_marginal);
        __m128 i_next_m0 = _mm_set1_ps(init_value);
        //set for diagonal
        prob_float[i] = LogMathFunctions::convert_to_space(0.00001);
        for (unsigned int j = 0; j < (size_j/4)+1; ++j) {
            __m128 tmp  = LogMathFunctions::mult_sse(prob_vec[j],j_cur_m1[j]);
            i_next_m0   = LogMathFunctions::add_sse(i_next_m0,tmp);
            
            tmp         = LogMathFunctions::mult_sse(prob_vec[j],i_cur_m0);
            j_next_m1[j]= LogMathFunctions::add_sse(j_next_m1[j],tmp);
            
        }
        //reset probability
        prob_float[i] = LogMathFunctions::convert_to_space(1.0f);
        // update i
        m0->add_to_next_marginal(i, (FP) LogMathFunctions::get_first (i_next_m0 ));
        m0->add_to_next_marginal(i, (FP) LogMathFunctions::get_second(i_next_m0 ));
        m0->add_to_next_marginal(i, (FP) LogMathFunctions::get_third (i_next_m0 ));
        m0->add_to_next_marginal(i, (FP) LogMathFunctions::get_fourth(i_next_m0 ));
    }
    // update j
    std::copy (j_next_m1_float, j_next_m1_float+size_j, m1_j_next_marginal );
    
}


IMPDOMINO3_END_NAMESPACE

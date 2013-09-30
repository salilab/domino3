#include <IMP/domino3/DistanceFactor.h>
#include <IMP/core/XYZR.h>

IMPDOMINO3_BEGIN_NAMESPACE
namespace {
    std::string get_distance_name(const kernel::ParticleIndexPair &pis) {
        std::ostringstream oss;
        oss << "Distance-" << pis[0] << "-" << pis[1];
        return oss.str();
    }
}

DistanceFactor::DistanceFactor(kernel::Model *m,
               const kernel::ParticleIndexPair &pis,
               FP distance, FP allowed_error,
               StatesTable *pst):
  Factor(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_distance_name(pis)), distance_(distance),
  allowed_error_(allowed_error) {
  // precompile allowed states, use some sort of location structure later
  this->pis = pis;
  this->pst = pst;
  this->distance = distance;
  this->m = m;
  this->allowed_error = allowed_error;
  States *ps0 = pst->get_states(pis[0]), *ps1 = pst->get_states(pis[1]);
  this->distances = new FP[ps0->get_number()*ps1->get_number()];
      std::cout << get_distance_name(pis) << std::endl;
  for (unsigned int i = 0; i < ps0->get_number(); ++i) {
  ps0->load(i, pis[0]);
     for (unsigned int j = 0; j < ps1->get_number(); ++j) {
              ps1->load(j, pis[1]);
         distances[i*ps0->get_number()+j] = IMP::core::get_distance(IMP::core::XYZR(m, pis[0]),
                                            IMP::core::XYZR(m, pis[1]));
     }
  }
}

FP DistanceFactor::distance_to_probability(FP x){
    x=x/this->allowed_error;
    return std::max(exp(-pow(x,2)*M_PI),0.00000001);
//    return std::max(1/sqrt(2*M_PI)*exp(-pow(x,2)*1/2),0.00000001);
}


void DistanceFactor::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
    States *ps0 = pst->get_states(pis[0]), *ps1 = pst->get_states(pis[1]);
    for (unsigned int i = 0; i < ps0->get_number(); ++i) {
        for (unsigned int j = 0; j < ps1->get_number(); ++j) {
            FP d = this->distances[i*ps0->get_number()+j];
            FP prob = distance_to_probability(d-this->distance);
             m0->add_to_next_marginal(i,log(prob)+m1->get_current_marginal(j));
             m1->add_to_next_marginal(j,log(prob)+m0->get_current_marginal(i));
        }
    }
}

IMPDOMINO3_END_NAMESPACE

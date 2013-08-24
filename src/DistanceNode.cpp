#include <IMP/domino3/DistanceNode.h>
#include <IMP/core/XYZR.h>

IMPDOMINO3_BEGIN_NAMESPACE
namespace {
    std::string get_distance_name(const kernel::ParticleIndexPair &pis) {
        std::ostringstream oss;
        oss << "Distance-" << pis[0] << "-" << pis[1];
        return oss.str();
    }
}

DistanceNode::DistanceNode(kernel::Model *m,
               const kernel::ParticleIndexPair &pis,
               double distance, double allowed_error,
               StatesTable *pst):
  Node(m, kernel::ParticleIndexes(pis.begin(), pis.end()),
       pst, get_distance_name(pis)), distance_(distance),
  allowed_error_(allowed_error) {
  // precompile allowed states, use some sort of location structure later
  this->pis = pis;
  this->pst = pst;
  this->distance = distance;
  this->m = m;
  this->allowed_error = allowed_error;
  States *ps0 = pst->get_states(pis[0]), *ps1 = pst->get_states(pis[1]);
  this->distances = new double[ps0->get_number()*ps1->get_number()];
      std::cout << get_distance_name(pis) << std::endl;
  for (unsigned int i = 0; i < ps0->get_number(); ++i) {
  ps0->load(i, pis[0]);
     for (unsigned int j = 0; j < ps1->get_number(); ++j) {
              ps1->load(j, pis[1]);
         distances[i*ps0->get_number()+j] = IMP::core::get_distance(IMP::core::XYZR(m, pis[0]),
                                            IMP::core::XYZR(m, pis[1]));
         std::cout << distances[i*ps0->get_number()+j] << " ";
     }
     std::cout << std::endl;
  }
}

double DistanceNode::distance_to_probability(double x){
    x=x/this->allowed_error;
    return std::max(exp(-pow(x,2)*M_PI),0.00000001);
}


void DistanceNode::do_update() {
    Marginals *m0 = get_marginals()[0], *m1 = get_marginals()[1];
    States *ps0 = pst->get_states(pis[0]), *ps1 = pst->get_states(pis[1]);
    for (unsigned int i = 0; i < ps0->get_number(); ++i) {
        for (unsigned int j = 0; j < ps1->get_number(); ++j) {
            double d = this->distances[i*ps0->get_number()+j];
            double cur = m0->mult_two_marginals(m0->get_current_marginal(i),
                                                    m1->get_current_marginal(j));
            double prob = distance_to_probability(d-this->distance);
                cur+= log(prob);
  //          m0->show_marginals();
  //          m1->show_marginals();

            m0->add_to_next_marginal(i, cur);
            m1->add_to_next_marginal(j, cur);
            
        }
    }
}

IMPDOMINO3_END_NAMESPACE

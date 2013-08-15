/** \file domino3/Marginals.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_MARGINAL_H
#define IMPDOMINO3_MARGINAL_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/domino3/MathFunctions.h>
#include <IMP/base/Object.h>
#include <IMP/base/ConstVector.h>
#include <IMP/base/object_macros.h>

IMPDOMINO3_BEGIN_NAMESPACE

class Marginals;

IMP_OBJECTS(Marginals, MarginalsList);



/** Store the marginal for a variable. */
class IMPDOMINO3EXPORT Marginals: public base::Object {
    boost::scoped_array<double> current_, next_;
    unsigned int size_;
    double change_;
    kernel::ParticleIndex pi_;
    MathFunctions  *math;
    
public:
    Marginals(kernel::Model *m, kernel::ParticleIndex pi, unsigned int size);
    
    double get_current_marginal(unsigned int state) const {
        IMP_USAGE_CHECK(state < this->get_number(),"state in get marginal is > size");
        return current_[state];
    }
    
    double mult_two_marginals(double val1,double val2) const {
        return math->mult(val1,val2);
    }
    
    void add_to_next_marginal(unsigned int state, double value) { //change to add next
        IMP_USAGE_CHECK(state < this->get_number(),"state in add marginal is > size");
        next_[state] = math->add(next_[state],value);
    }
    
    void multp_to_next_marginal(unsigned int state, double value) {
        IMP_USAGE_CHECK(state < this->get_number(),"state in multp marginal is > size");
        next_[state] = math->mult(next_[state],value);
    }
    
    
    void set_uniform();
    
    void set_random();
    
    void show_marginals(){
        std::cout << this->get_name() << std::endl;
        for(int i = 0; i < this->get_number(); i++){
            std::cout << this->current_[i] << " ";
        }
        std::cout << std::endl;
        for(int i = 0; i < this->get_number(); i++){
            std::cout << this->next_[i] << " ";
        }
        std::cout << std::endl;

    }
    
    void calculate_joint_probability(boost::scoped_array<double> &array, const Marginals *marginals){
        for (unsigned int j = 0; j < this->get_number(); ++j) {
            array[j] = math->mult(array[j], marginals->current_[j]);
        }
    }
    
    void check_current_normalized(){
        double total = math->convert_to_space(0);
        for(int i =0;i < size_;i++){
            total = math->add(total,current_[i]);
        }
//        IMP_USAGE_CHECK(std::abs(total - 1.0)  < 0.01,
//                        "Not normalized" << total);
    }
    
    void check_next_normalized(){
        double total = std::accumulate(next_.get(), next_.get() + size_, 0.0);
//        IMP_USAGE_CHECK(std::abs(total - 1.0)  < 0.01,
//                        "Not normalized" << total);
    }
    
    void make_next_zero(){
        double zero_value = math->convert_to_space(0.0);
        std::fill(next_.get(), next_.get() + this->get_number(),zero_value); 
    }
    
    void merge_probabilities_from_list(const MarginalsListTemp &others) {
        std::copy(current_.get(), current_.get() + this->get_number(),next_.get());
        check_current_normalized();
        for (unsigned int i = 0; i < others.size(); ++i) {
            IMP_USAGE_CHECK(others[i]->get_number() == this->get_number(), "size not match");
            calculate_joint_probability(next_,others[i].get());
            math->normalize(next_.get(), this->get_number());
        }
        set_current_from_next();
    }
    /** Eventually this will be atomic. */
    void set_current_from_next() {
        math->normalize(next_.get(), this->get_number());
        check_current_normalized();
        change_ = 0;
        for (unsigned int i = 0; i < this->get_number(); ++i) {
            change_ += std::abs(math->convert_to_linear(next_[i]) - math->convert_to_linear(current_[i]));
        }
        using namespace std;
        swap(current_, next_);
    }
    
    /** Return a metric on the change (currently L0, could change) */
    double get_change() const {
        return change_;
    }
    
    kernel::ParticleIndex get_particle_index() const {
        return pi_;
    }
    
    unsigned int get_number() const {
        return size_;
    }
    
    double get_entropy() const;
    
    IMP_OBJECT_METHODS(Marginals);
};

IMPDOMINO3_END_NAMESPACE


#endif // IMPDOMINO3_MARGINAL_H

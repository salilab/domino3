/** \file IMP/domino3/Marginals.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_MARGINALS_H
#define IMPDOMINO3_MARGINALS_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/domino3/LogMathFunctions.h>
#include <IMP/base/Object.h>
#include <IMP/base/ConstVector.h>
#include <IMP/base/object_macros.h>

IMPDOMINO3_BEGIN_NAMESPACE

class Marginals;

IMP_OBJECTS(Marginals, MarginalsList);



/** Store the marginal for a variable. */
class IMPDOMINO3EXPORT Marginals: public base::Object {
    boost::scoped_array<FP> current_, next_;
    unsigned int size_;
    FP change_;
    kernel::ParticleIndex pi_;
    
public:
    Marginals(kernel::Model *m, kernel::ParticleIndex pi, unsigned int size);
    
    inline FP * get_current_marginal(){
        return current_.get();
    }
    
    inline FP get_current_marginal(unsigned int state) const {
        IMP_USAGE_CHECK(state < this->get_number(),"state in get marginal is > size");
        return current_[state];
    }
    
    inline FP * get_next_marginal(){
        return next_.get();
    }
    
    inline FP get_next_marginal(unsigned int state) const {
        IMP_USAGE_CHECK(state < this->get_number(),"state in get marginal is > size");
        return next_[state];
    }
    
    inline FP mult_two_marginals(FP val1,FP val2) const {
        return LogMathFunctions::mult(val1,val2);
    }
    
    inline void add_to_next_marginal(unsigned int state, FP value) { //change to add next
        IMP_USAGE_CHECK(state < this->get_number(),"state in add marginal is > size");
        next_[state] = LogMathFunctions::add(next_[state],value);
    }
    
    inline void multp_to_next_marginal(unsigned int state, FP value) {
        IMP_USAGE_CHECK(state < this->get_number(),"state in multp marginal is > size");
        next_[state] = LogMathFunctions::mult(next_[state],value);
    }
    
    
    void set_uniform();
    
    void set_random();
    
    void set_init_vector(boost::scoped_array<FP> &array);
    
    inline FP convert_to_space(FP x ){
        return LogMathFunctions::convert_to_space(x);
    }
    
    void show_marginals(){
        std::cout << this->get_name() << std::endl;
        for(unsigned i = 0; i < this->get_number(); i++){
            std::cout << this->current_[i] << " ";
        }
        std::cout << std::endl;
        for(unsigned i = 0; i < this->get_number(); i++){
            std::cout << this->next_[i] << " ";
        }
        std::cout << std::endl;

    }
    
    inline void compute_joint_probability(boost::scoped_array<FP> &array, const Marginals *marginals){
        for (unsigned int j = 0; j < this->get_number(); ++j) {
            array[j] = LogMathFunctions::mult(array[j], marginals->current_[j]);
        }
    }
    
    inline void compute_average_probability(boost::scoped_array<FP> &array, const Marginals *marginals){
        FP avg_term = LogMathFunctions::convert_to_space(2.0);

        for (unsigned int j = 0; j < this->get_number(); ++j) {
            array[j] = LogMathFunctions::dev(LogMathFunctions::add(array[j], marginals->current_[j]),avg_term);
        }
    }
    
    inline void check_current_normalized(){
        FP total = LogMathFunctions::convert_to_space(0);
        for(unsigned int i =0;i < size_;i++){
            total = LogMathFunctions::add(total,current_[i]);
        }
//        IMP_USAGE_CHECK(std::abs(total - 1.0)  < 0.01,
//                        "Not normalized" << total);
    }
    
    inline void check_next_normalized(){
//      FP total = std::accumulate(next_.get(), next_.get() + size_, 0.0);
//        IMP_USAGE_CHECK(std::abs(total - 1.0)  < 0.01,
//                        "Not normalized" << total);
    }
    
    inline void make_next_zero(){
        FP zero_value = LogMathFunctions::convert_to_space(0.0);
        std::fill(next_.get(), next_.get() + this->get_number(),zero_value); 
    }
    
    inline void merge_probabilities_from_list(const MarginalsListTemp &others) {
        std::copy(current_.get(), current_.get() + this->get_number(),next_.get());
        check_current_normalized();
        for (unsigned int i = 0; i < others.size(); ++i) {
            IMP_USAGE_CHECK(others[i]->get_number() == this->get_number(), "size not match");
            compute_joint_probability(next_,others[i].get());
            LogMathFunctions::normalize(next_.get(), this->get_number());
        }
        set_current_from_next();
    }
    /** Eventually this will be atomic. */
    inline void set_current_from_next() {
        LogMathFunctions::normalize(next_.get(), this->get_number());
        check_current_normalized();
        change_ = 0;
        for (unsigned int i = 0; i < this->get_number(); ++i) {
            change_ += std::abs(LogMathFunctions::convert_to_linear(next_[i]) - LogMathFunctions::convert_to_linear(current_[i]));
        }
        using namespace std;
        swap(current_, next_);
    }
    
    //! Return a metric on the change (currently L0, could change)
    FP get_change() const {
        return change_;
    }
    
    kernel::ParticleIndex get_particle_index() const {
        return pi_;
    }
    
    unsigned int get_number() const {
        return size_;
    }
    
    FP get_entropy() const;
    
    IMP_OBJECT_METHODS(Marginals);
};

IMPDOMINO3_END_NAMESPACE


#endif /* IMPDOMINO3_MARGINALS_H */

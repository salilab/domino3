#include <IMP/domino3/Updater.h>

IMPDOMINO3_BEGIN_NAMESPACE

Updater::Updater(const FactorsTemp &graph,
                 std::string name):
base::Object(name), factors_(graph) {
    for (unsigned int i = 0; i < factors_.size(); ++i) {
        factors_[i]->set_index(i);
    }
    set_change_threshold(.01);
    fill_queue();
}

float Updater::do_update() {
    base::set<kernel::ParticleIndex> changed;
    std::random_shuffle ( cur_queue_.begin(), cur_queue_.end() );
    float sum_change = 0.0f;
    for (ActiveSet::const_iterator it = cur_queue_.begin();
         it != cur_queue_.end(); ++it) {
        //    IMP_LOG_TERSE("Updating " << (*it)->get_name() << std::endl);
        (*it)->update();
        for (unsigned int i = 0; i < (*it)->get_marginals().size(); ++i) {
            if ((*it)->get_marginals()[i]->get_change() > change_threshold_) {
                changed.insert((*it)->get_particle_indexes()[i]);
            }
            sum_change += (*it)->get_marginals()[i]->get_change();
        }
    }
    // later search only relevant nodes using the graph
    for (unsigned int i = 0; i< factors_.size(); ++i) {
        for (unsigned int j = 0; j < factors_[i]->get_particle_indexes().size(); ++j) {
            if (changed.find(factors_[i]->get_particle_indexes()[j]) != changed.end()) {
                add_factor_to_next_set(factors_[i]);
            }
        }
    }
    swap_active_sets();
    next_queue_ = ActiveSet();
    return sum_change;
}

void Updater::fill_queue() {
    // lazy
    for (unsigned int i = 0; i < factors_.size(); ++i) {
        add_factor_to_next_set(factors_[i]);
    }
    swap_active_sets();
}

void Updater::update(unsigned int iterations) {
    float before_change = 0.0f;
    for (unsigned int i = 0; i < iterations; ++i) {
        std::cout << "Iteration " << i << std::endl;
        float current_change = do_update();
        std::cout << "Change: " << std::abs(current_change - before_change) << std::endl;
        //    if(std::abs(current_change - before_change) < change_threshold_){
        //       break;
        //    }
        before_change = current_change;
    }
}

IMPDOMINO3_END_NAMESPACE

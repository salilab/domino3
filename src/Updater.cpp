/**
 *  \file IMP/domino3/Updater.cpp
 *
 *  Copyright 2007-2015 IMP Inventors. All rights reserved.
 *
 */

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

void Updater::do_update() {
    boost::unordered_set<ParticleIndex> changed;
    std::random_shuffle ( cur_queue_.begin(), cur_queue_.end() );
    this->change  = 0.0f;
    this->entropy  = 0.0f;
    for (ActiveSet::const_iterator it = cur_queue_.begin();
         it != cur_queue_.end(); ++it) {
        //    IMP_LOG_TERSE("Updating " << (*it)->get_name() << std::endl);
        (*it)->update();
        for (unsigned int i = 0; i < (*it)->get_marginals().size(); ++i) {
            if ((*it)->get_marginals()[i]->get_change() > change_threshold_) {
                changed.insert((*it)->get_particle_indexes()[i]);
            }
            this->change  += (*it)->get_marginals()[i]->get_change();
            this->entropy += (*it)->get_marginals()[i]->get_entropy();
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
    
}

void Updater::fill_queue() {
    // lazy
    for (unsigned int i = 0; i < factors_.size(); ++i) {
        add_factor_to_next_set(factors_[i]);
    }
    swap_active_sets();
}

void Updater::update(unsigned int iterations) {
    float before_change  = 0.0f;
    float before_entropy = 0.0f;
    for (unsigned int i = 0; i < iterations; ++i) {
        IMP_LOG_VERBOSE("Iteration " << i << std::endl);
        do_update();
        float entropy_change = std::abs(this->entropy - before_entropy);
        float current_change = std::abs(this->change  - before_change);
        IMP_LOG_VERBOSE("Change: "  << current_change << std::endl);
        IMP_LOG_VERBOSE("Entropy: " << entropy_change << std::endl);
        if(entropy_change < 0.01 && current_change < 0.01 )
            break;
        before_change  = this->change ;
        before_entropy = this->entropy;
    }
}

IMPDOMINO3_END_NAMESPACE

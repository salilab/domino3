/**
 * Copyright 2007-2013 IMP Inventors. All rights reserved.
 */

#include <IMP/base/flags.h>
#include <IMP/kernel/Model.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/hierarchy_tools.h>
#include <IMP/atom/pdb.h>


#include <IMP/rmf/atom_io.h>
#include <IMP/rmf/restraint_io.h>
#include <IMP/rmf/frames.h>
#include <IMP/rmf/particle_io.h>

#include <IMP/base/SetLogState.h>
#include <IMP/domino3/Factor.h>
#include <IMP/domino3/Probability2DFactor.h>
#include <IMP/domino3/Probability3DFactor.h>
#include <IMP/domino3/ExcludedVolumeFactor.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <IMP/domino3/Updater.h>
#include <RMF/FileHandle.h>
#include <IMP/domino3/IndexStates.h>
#include <IMP/core/XYZ.h>
#include <algorithm>

namespace {
    std::string output = "out.rmf";
    IMP::base::AddStringFlag oasf("output", "Output rmf name", &output);
    std::string dock_file = IMP::domino3::get_example_path("33_dock.txt");
    IMP::base::AddStringFlag dock("dock", "Dock score file",&dock_file);
    std::string sea_file = IMP::domino3::get_example_path("4_true_sea.txt");
    IMP::base::AddStringFlag sea("sea", "SEA score file",&sea_file);
    std::string chem_file = IMP::domino3::get_example_path("4_33_chemsim_rxn.txt");
    IMP::base::AddStringFlag chem("chem", "Chem score file",&chem_file);
    
    boost::int64_t iterations = 100;
    IMP::base::AddIntFlag aif("iterations", "Number of iterations",
                              &iterations);
    boost::int64_t enzyme_size=3;
    IMP::base::AddIntFlag enzymes("enzymes", "Number of enzymes",
                                  &enzyme_size);
    boost::int64_t ligand_size=4;
    IMP::base::AddIntFlag ligands("ligands", "Number of ligands",
                                  &ligand_size);
    static int enzym_id_counter = 0;
    struct Enzyme {
        std::string enzyme_name;
        const unsigned id = enzym_id_counter++;
        Enzyme(std::string name) : enzyme_name(name){}
    };
    
    static int ligand_id_counter = 0;
    struct Ligand {
        std::string ligand_name;
        const unsigned id = ligand_id_counter++;
        Ligand(std::string name) : ligand_name(name){}
    };
    typedef std::pair<Enzyme *,Enzyme *> EnzymeEnzymePair;
    typedef std::pair<Enzyme *,Ligand *> EnzymeLigandPair;
    typedef std::pair<Ligand *,Ligand *> LigandLigandPair;
    typedef std::pair<Enzyme *,LigandLigandPair> EnzymeLigandLigandPair;

    typedef IMP::base::map<EnzymeEnzymePair, double > EnzymeEnzymeScoreLookup;
    typedef IMP::base::map<EnzymeLigandPair, double > EnzymeLigandScoreLookup;
    typedef IMP::base::map<EnzymeLigandLigandPair, double> EnzymeLigandLigandScoreLookup;
    typedef IMP::base::map<std::string, Enzyme *>  StringToEnzyme;
    typedef IMP::base::map<std::string, Ligand *>  StringToLigand;
    
    void read_in_sea_file(IMP::kernel::Model *m,
                           std::string path,
                           EnzymeEnzymeScoreLookup &scores,
                           StringToEnzyme &id_to_enzyme ){
        std::cout << "Read sea file: " << path << std::endl;

        std::ifstream infile(path.c_str());
        std::string line;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            std::string a, b;
            double score;
            iss >> a >> b >> score;
            if(!id_to_enzyme[a])
                id_to_enzyme[a] = new Enzyme(a);
            if(!id_to_enzyme[b])
                id_to_enzyme[b] = new Enzyme(b);
            scores[std::make_pair(id_to_enzyme[a],id_to_enzyme[b])]=score;
        }
    }
    
    
    void read_in_dock_file(IMP::kernel::Model *m,
                           std::string path,
                           EnzymeLigandScoreLookup &scores,
                           StringToEnzyme &id_to_enzyme,
                           StringToLigand &id_to_ligand){
        std::cout << "Read dock file: " << path << std::endl;
        std::ifstream infile(path.c_str());
        std::string line;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            std::string a, b;
            double score;
            iss >> a >> b >> score;
            if(!id_to_enzyme[a])
                id_to_enzyme[a] = new Enzyme(a);
            if(!id_to_ligand[b])
                id_to_ligand[b] = new Ligand(b);
            scores[std::make_pair(id_to_enzyme[a],id_to_ligand[b])]=score;
        }
    }
    
    
    void read_in_chem_file(IMP::kernel::Model *m,
                           std::string path,
                           EnzymeLigandLigandScoreLookup &scores,
                           StringToEnzyme &id_to_enzyme,
                           StringToLigand &id_to_ligand){
        std::cout << "Read chem sim. file: " << path << std::endl;
        std::ifstream infile(path.c_str());
        std::string line;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            std::string a, b, c;
            double score;
            iss >> a >> b >> c >> score;
            if(!id_to_enzyme[a])
                id_to_enzyme[a] = new Enzyme(a);
            if(!id_to_ligand[b])
                id_to_ligand[b] = new Ligand(b);
            if(!id_to_ligand[c])
                id_to_ligand[c] = new Ligand(c);
            LigandLigandPair ligand_ligand_pair = std::make_pair(id_to_ligand[b],id_to_ligand[c]);
            scores[std::make_pair(id_to_enzyme[a],ligand_ligand_pair)]=score;
        }
    }
    
    
    void create_linear_graph(IMP::domino3::Factors &factors,
                             IMP::domino3::StatesTable * st,
                             IMP::kernel::Model *m,
                             int enzyme_size,int ligand_size,
                             boost::shared_array<double> &sea_probability_vector,
                             boost::shared_array<double> &dock_probability_vector,
                             boost::shared_array<double> &chem_probability_vector){
        IMP::domino3::FactorEdgesTemp edges;
        IMP::kernel::ParticleIndexes pis = st->get_particle_indexes();
        int ENZYME_SIZE=enzyme_size;
        int LIGAND_SIZE=ligand_size;
        // create SEA factors 
        for(int i = 0; i < ENZYME_SIZE-1; i++) {
            IMP::kernel::ParticleIndexPair cur_pair1(pis[i], pis[i+1]);
            IMP_NEW(IMP::domino3::Probability2DFactor, pf1,(m, cur_pair1, st, sea_probability_vector));
            factors.push_back(pf1);
            pf1->set_was_used(true);
        }
        // create dock factors
        for(int i = 0; i < ENZYME_SIZE; i++) {
            IMP::kernel::ParticleIndexPair cur_pair3(pis[i], pis[ENZYME_SIZE+i]);
            IMP_NEW(IMP::domino3::Probability2DFactor, pf3,(m, cur_pair3, st, dock_probability_vector));
            factors.push_back(pf3);
            pf3->set_was_used(true);
        }
        // create chem. sim. factors
        for(int i = 0; i < ENZYME_SIZE; i++) {
            IMP::kernel::ParticleIndexTriplet cur_triplet1(pis[i], pis[ENZYME_SIZE+i], pis[ENZYME_SIZE+i+1]);
            IMP_NEW(IMP::domino3::Probability3DFactor, pf3d1,(m, cur_triplet1, st, chem_probability_vector));
            factors.push_back(pf3d1);
            pf3d1->set_was_used(true);
        }
        // create exclude double enzyme
        for(int x = 0 ; x< ENZYME_SIZE; x++){
            for(int y = 0; y < x ; y++){
                IMP::kernel::ParticleIndexPair excl_pair(pis[x], pis[y]);
                IMP_NEW(IMP::domino3::ExcludedVolumeFactor, excl_factor,(m, excl_pair, st));
                factors.push_back(excl_factor);
                excl_factor->set_was_used(true);
            }
        }
        // create exclude double ligand
        for(int x = ENZYME_SIZE ; x< (ENZYME_SIZE+LIGAND_SIZE); x++){
            for(int y = ENZYME_SIZE; y < x ; y++){
                IMP::kernel::ParticleIndexPair excl_pair(pis[x], pis[y]);
                IMP_NEW(IMP::domino3::ExcludedVolumeFactor, excl_factor,(m, excl_pair, st));
                factors.push_back(excl_factor);
                excl_factor->set_was_used(true);
            }
        }        
        std::cout << "Node size: " << factors.size() << std::endl;
        IMP::domino3::add_neighbors(factors);
    }
    
    void normalize(boost::shared_array<double> &probability_vector,
                   int size_x,
                   int size_y,int size_z,
                   bool normalize_square){
        std::cout << "Normalize probability vector of size: " << " x: " << size_x << " y: " << size_y << " z: " << size_z << std::endl;

        double total = 0;
        for(int x = 0; x < size_x; x++){
            double * it = probability_vector.get()+(x*size_y*size_z);
            if(normalize_square)
                total = std::accumulate(it, it + (size_z*size_y), 0.0);
            for(int y = 0; y < size_y; y++){
                double * row = it + (y*size_y);
                if(!normalize_square)
                    total = std::accumulate(row, row + size_z, 0.0);
                for(int z = 0; z < size_z; z++){
                    *(row+z)=log(*(row+z)/total);
                }
            }
        }
    }
    void run_it(IMP::kernel::Model *m) {
        StringToEnzyme id_to_enzyme;
        EnzymeEnzymeScoreLookup enzyme_enzyme_scores;
        StringToLigand id_to_ligand;
        EnzymeLigandScoreLookup enzyme_ligand_scores;
        EnzymeLigandLigandScoreLookup enzyme_ligand_ligand_scores;
        IMP_USAGE_CHECK(boost::filesystem::exists(sea_file) , "SEA File dosen't exist");
        IMP_USAGE_CHECK(boost::filesystem::exists(dock_file), "Dock File dosen't exist");
        IMP_USAGE_CHECK(boost::filesystem::exists(chem_file), "Chem File dosen't exist");

        read_in_sea_file(m,
                          sea_file,
                          enzyme_enzyme_scores,
                          id_to_enzyme);
        int enzymes_in_sea=id_to_enzyme.size();
        read_in_dock_file(m,
                         dock_file,
                         enzyme_ligand_scores,
                         id_to_enzyme,
                         id_to_ligand);
        int enzymes_in_dock=id_to_enzyme.size();
        int ligands_in_dock=id_to_ligand.size();
        read_in_chem_file(m,
                          chem_file,
                          enzyme_ligand_ligand_scores,
                          id_to_enzyme,
                          id_to_ligand);
        int enzymes_in_chem=id_to_enzyme.size();
        int ligands_in_chem=id_to_ligand.size();
        IMP_USAGE_CHECK(enzymes_in_chem == enzymes_in_sea,"Not the same amount of enzymes in SEA and Chem.");
        IMP_USAGE_CHECK(enzymes_in_dock == enzymes_in_sea,"Not the same amount of enzymes in SEA and Dock.");
        IMP_USAGE_CHECK(ligands_in_chem == ligands_in_dock,"Not the same amount of ligands in Dock and Chem.");


        std::cout << id_to_enzyme.size() << std::endl;
        std::vector<int> enzyme_ids;
        std::vector<int> ligand_ids;

        
        // read structure to probability_vector and normalize
        boost::shared_array<double> sea_probability_vector(new double[id_to_enzyme.size()*id_to_enzyme.size()]);
        std::fill (sea_probability_vector.get(),sea_probability_vector.get()+(id_to_enzyme.size()*id_to_enzyme.size()),0);
        for (EnzymeEnzymeScoreLookup::iterator i = enzyme_enzyme_scores.begin(); i != enzyme_enzyme_scores.end(); ++i)
        {
            EnzymeEnzymePair enzymePair=i->first;
            double score = i->second;
            Enzyme * e1 = enzymePair.first;
            Enzyme * e2 = enzymePair.second;
            if(e1->id == e2->id){
                sea_probability_vector[e1->id*id_to_enzyme.size()+e2->id] = score/10;
            }else{
                sea_probability_vector[e1->id*id_to_enzyme.size()+e2->id] = score;
            }
        }
        normalize(sea_probability_vector,1,id_to_enzyme.size(),id_to_enzyme.size(),true);

        boost::shared_array<double> dock_probability_vector(new double[id_to_enzyme.size()*id_to_ligand.size()]);
        std::fill (dock_probability_vector.get(),dock_probability_vector.get()+(id_to_enzyme.size()*id_to_ligand.size()),0);
        for (EnzymeLigandScoreLookup::iterator i = enzyme_ligand_scores.begin(); i != enzyme_ligand_scores.end(); ++i)
        {
            EnzymeLigandPair enzymePair=i->first;
            double score = i->second;
            Enzyme * e1 = enzymePair.first;
            Ligand * l1 = enzymePair.second;
            dock_probability_vector[e1->id*id_to_ligand.size()+l1->id] = score;
//            dock_probability_vector[e1->id*id_to_ligand.size()+l1->id] =-3.49650756146648;
        }
//        normalize(dock_probability_vector,1,id_to_enzyme.size(),id_to_ligand.size(),false);

        const unsigned int chem_size = id_to_enzyme.size()*id_to_ligand.size()*id_to_ligand.size();
        boost::shared_array<double> chem_probability_vector(new double[chem_size]);
        std::fill (chem_probability_vector.get(),chem_probability_vector.get()+chem_size,0);
        for (EnzymeLigandLigandScoreLookup::iterator i = enzyme_ligand_ligand_scores.begin(); i != enzyme_ligand_ligand_scores.end(); ++i)
        {
            EnzymeLigandLigandPair enzymePair=i->first;
            double score = i->second;
            Enzyme * e1 = enzymePair.first;
            LigandLigandPair ligand_pair = enzymePair.second;
            Ligand * l1 = ligand_pair.first;
            Ligand * l2 = ligand_pair.second;
            unsigned int x = e1->id;
            unsigned int y = l1->id;
            unsigned int z = l2->id;

            unsigned int index = x*id_to_ligand.size()*id_to_ligand.size() + y*id_to_ligand.size() + z;
            if(y == z)
                chem_probability_vector[index] = score/100;
            else
                chem_probability_vector[index] = score;
        }
        normalize(chem_probability_vector,id_to_enzyme.size(),id_to_ligand.size(),id_to_ligand.size(),true);


        for (StringToEnzyme::iterator i = id_to_enzyme.begin(); i != id_to_enzyme.end(); ++i){
            Enzyme * p = i->second;
            enzyme_ids.push_back(p->id);
        }
        for (StringToLigand::iterator i = id_to_ligand.begin(); i != id_to_ligand.end(); ++i){
            Ligand * p = i->second;
            ligand_ids.push_back(p->id);
        }
        IMP_NEW(IMP::domino3::StatesTable, st, (m));
        IMP_NEW(IMP::domino3::IndexStates, enzyme_states, (m, enzyme_ids));
        IMP_NEW(IMP::domino3::IndexStates, ligand_states, (m, ligand_ids));
        // Enzymes
        for (int i = 0; i < enzyme_size; i++){
            IMP::kernel::Particle * p = new IMP::kernel::Particle(m);
            IMP_NEW(IMP::domino3::Marginals, marginals, (m, p->get_index(), id_to_enzyme.size()));
            marginals->set_uniform();
            st->add(p->get_index(), enzyme_states, marginals);
        }
        // Ligands
        for (int i = 0; i < ligand_size; i++){
            IMP::kernel::Particle * p = new IMP::kernel::Particle(m);
            IMP_NEW(IMP::domino3::Marginals, marginals, (m, p->get_index(), id_to_ligand.size()));
            marginals->set_uniform();
            st->add(p->get_index(), ligand_states, marginals);
        }
        st->print_marginal();

        IMP::domino3::Factors factors;
        create_linear_graph(factors,st,m,
                            enzyme_size,ligand_size,
                            sea_probability_vector,
                            dock_probability_vector,
                            chem_probability_vector);
        
        //
        IMP_NEW(IMP::domino3::Updater, ud, (factors, "updater"));
        ud->update(iterations);
        IMP::domino3::update_state_table(factors,st);
        IMP::domino3::print_graph(factors);
        //
        st->print_marginal();
    }
}



int main(int argc, char **argv) {
    IMP::base::setup_from_argv(argc, argv, "Experiment with loopy domino");
    IMP_NEW(IMP::kernel::Model, m, ());
    run_it(m);
    return 0;
}

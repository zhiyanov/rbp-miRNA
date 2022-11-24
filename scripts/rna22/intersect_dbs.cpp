#include <string>
#include <set>
#include <vector>
#include <fstream>
#include <iostream>


std::set<std::string>* fill_set(std::string db_path) {
    std::fstream fstr;
    fstr.open(
        db_path,
        std::ios::in
    );

    std::set<std::string>* db_set = new std::set<std::string>();
    std::string str;
    int iteration_number = 0;
    while (std::getline(fstr, str)) {
        if (iteration_number % 100000 == 0) {
            std::cout << iteration_number << "\n";
        }
        
        iteration_number += 1;
        
        db_set->insert(str);
    }

    return db_set;
}

int main() {
    auto first_db_set = fill_set(
        "/home/dude/huge/dude/rbp-miRNA/rna22/rna22_energy_seed_pvalue.tsv");
    auto second_db_set = fill_set(
        "/home/dude/huge/dude/rbp-miRNA/rna22/rna22_heteroduplex_seed_pvalue.tsv");
    
    int intersection_size = 0;
    for (std::string str: *first_db_set) { 
        if (second_db_set->find(str) != second_db_set->end()) {
            intersection_size += 1;
        }
    }
    
    std::cout << first_db_set->size() << " ";
    std::cout << intersection_size << " ";
    std::cout << second_db_set->size() << "\n";

    delete first_db_set;
    delete second_db_set;
    return 0;
}

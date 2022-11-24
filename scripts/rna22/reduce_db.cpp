#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <set>
#include <string>
#include <functional>

typedef std::tuple<std::string, std::string,
        std::string, std::string> t4;
typedef std::function<bool(std::string)> cond;


int reduce_rna22_by_condition(
        std::string path_to_rna22,
        std::string path_to_output,
        cond condition) {
    
    std::fstream rna22_db, rna22_reduced_db;
    rna22_db.open(
        path_to_rna22,
        std::ios::in
    );
    rna22_reduced_db.open(
        path_to_output,
        std::ios::out
    );

    if (!rna22_db.is_open()) {
        return -1;
    }
    
    if (!rna22_reduced_db.is_open()) {
        return -2;
    }
    
    std::cout << "Reduction was started\n";

    std::string peak;
    int iteration_number = 0, success_number = 0;
    while (std::getline(rna22_db, peak)) {
        if (iteration_number % 100000 == 0) {
            std::cout << iteration_number << " " << success_number << "\n";
        }
        
        iteration_number += 1;

        if (!condition(peak)) {
            continue;
        }
        
        rna22_reduced_db << peak << "\n";
        success_number += 1;
    }

    rna22_db.close();
    rna22_reduced_db.close();

    return 0;
}

std::set<t4> create_gene_transcript_chr_strand_set(
        std::string path_to_gtf_clip) {
    
    std::fstream gtf_clip_db;
    gtf_clip_db.open(
        path_to_gtf_clip,
        std::ios::in
    );
    
    std::set<t4> gtf_clip_set;

    std::string peak;
    bool header_flag = true;
    while (std::getline(gtf_clip_db, peak)) {
        if (header_flag) {
            header_flag = false;
            continue;
        }

        std::istringstream iss(peak);
        std::vector<std::string> tokens;
        std::string token;

        while (std::getline(iss, token, '\t')) {
            tokens.push_back(token);
        }

        t4 sp{
            tokens[9], tokens[10],
            tokens[2], tokens[3]
        };
        
        gtf_clip_set.insert(sp); 
    }
    
    gtf_clip_db.close();
    
    return gtf_clip_set;
}

int reduce_rna22_by_gene_transcript_chr_strand() {
    std::set<t4> gtf_clip_set = create_gene_transcript_chr_strand_set(
        "/home/dude/huge/dude/rbp-miRNA/clip/clip_db_annotated.tsv"
    );
    
    std::cout << "Reduction by gene, transcript, chr, strand\n";

    auto gtf_set_condition = [gtf_clip_set](std::string peak) {
        std::istringstream _iss(peak);
        std::vector<std::string> _tokens;
        std::string _token;

        while (std::getline(_iss, _token, '\t')) {
            _tokens.push_back(_token);
        }
    
        std::istringstream iss(_tokens[1]);
        std::vector<std::string> tokens;
        std::string token;

        while (std::getline(iss, token, '_')) {
            tokens.push_back(token);
        }

        if (tokens[3] == "1") {
            tokens[3] = "+";
        } else {
            tokens[3] = "-";
        }

        t4 sp{
            tokens[0], tokens[1],
            tokens[2], tokens[3]
        };
        
        if (gtf_clip_set.find(sp) == gtf_clip_set.end()) {
            return false;
        }
        
        return true;
    };

    return reduce_rna22_by_condition(
        "/home/dude/huge/bulk/RNA22/Homo_Sapiens_mRNA_ENSEMBL96_mirbasev22.txt",
        "/home/dude/huge/dude/rbp-miRNA/rna22/rna22_gene_transcript_chr_strand.tsv",
        gtf_set_condition
    );
}

int reduce_rna22_by_energy_seed_pvalue() {
    float energy_bound = -16.0;
    int seed_matching_bound = 6;
    float pvalue_bound = 0.01;
    
    std::cout << "Reduction by energy, seed, pvalue\n";

    auto condition = [
        energy_bound,
        seed_matching_bound,
        pvalue_bound
    ](std::string peak) {

        std::istringstream _iss(peak);
        std::vector<std::string> _tokens;
        std::string _token;

        while (std::getline(_iss, _token, '\t')) {
            _tokens.push_back(_token);
        }

        float energy = std::stof(_tokens[5]);
        float pvalue = std::stof(_tokens[15]);
        std::string miRNA_align = _tokens[9];

        int miRNA_align_len = miRNA_align.length();
        for (int i = 0; i < miRNA_align_len / 2; ++i) {
            char tmp = miRNA_align[i];
            miRNA_align[i] = miRNA_align[miRNA_align_len - i - 1]; 
            miRNA_align[miRNA_align_len - i - 1] = tmp;
        }
        
        int seed_mathcing_number = 0;
        for (int i = 1; i < 7; ++i) {
            if (miRNA_align[i] == ')') {
                seed_mathcing_number += 1;
            }
        }

        if (pvalue > pvalue_bound) {
            return false;
        }

        if (energy > energy_bound) {
            return false;
        }

        if (seed_mathcing_number < seed_matching_bound) {
            return false;
        }
        
        return true;
    };
    
    return reduce_rna22_by_condition(
        "/home/dude/huge/bulk/RNA22/Homo_Sapiens_mRNA_ENSEMBL96_mirbasev22.txt",
        "/home/dude/huge/dude/rbp-miRNA/rna22/rna22_energy_seed_pvalue.tsv",
        condition
    ); 
}

int reduce_rna22_by_heteroduplex_seed_pvalue() {
    int heteroduplex_matching_bound = 12;
    int seed_matching_bound = 5;
    float pvalue_bound = 0.01;

    std::cout << "Reduction by heteroduplex, seed, pvalue\n";

    auto condition = [
        heteroduplex_matching_bound,
        seed_matching_bound,
        pvalue_bound
    ](std::string peak) {

        std::istringstream _iss(peak);
        std::vector<std::string> _tokens;
        std::string _token;

        while (std::getline(_iss, _token, '\t')) {
            _tokens.push_back(_token);
        }

        float pvalue = std::stof(_tokens[15]);
        std::string mRNA_align = _tokens[8];
        std::string miRNA_align = _tokens[9];


        int miRNA_align_len = miRNA_align.length();
        for (int i = 0; i < miRNA_align_len / 2; ++i) {
            char tmp = miRNA_align[i];
            miRNA_align[i] = miRNA_align[miRNA_align_len - i - 1]; 
            miRNA_align[miRNA_align_len - i - 1] = tmp;
        }
        
        int seed_mathcing_number = 0;
        for (int i = 1; i < 7; ++i) {
            if ((miRNA_align[i] == ')') && ((mRNA_align[i] == '('))) {
                seed_mathcing_number += 1;
            }
        }
           
        int heteroduplex_matching_number = 0;
        for (int i = 0; i < miRNA_align_len; ++i) {
            if (miRNA_align[i] == ')') {
                heteroduplex_matching_number += 1;
            }
        } 

        if (pvalue > pvalue_bound) {
            return false;
        }

        if (seed_mathcing_number < seed_matching_bound) {
            return false;
        }

        if (heteroduplex_matching_number < heteroduplex_matching_bound) {
            return false;
        }
        
        return true;
    };
    
    return reduce_rna22_by_condition(
        "/home/dude/huge/bulk/RNA22/Homo_Sapiens_mRNA_ENSEMBL96_mirbasev22.txt",
        "/home/dude/huge/dude/rbp-miRNA/rna22/rna22_heteroduplex_seed_pvalue.tsv",
        condition
    ); 
}

int main() {
    reduce_rna22_by_gene_transcript_chr_strand();
    // reduce_rna22_by_energy_seed_pvalue();
    // reduce_rna22_by_heteroduplex_seed_pvalue();
    return 0;
}

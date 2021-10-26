

//////////////// Generating mRNA
void generating_mRNA(){
    // for (unsigned long long i = 0; i < seq_len; i++)
    //     if (seq[i] == 'T')
    //         seq[i] = 'U';
}


//////////////// Detecting genes
void detecting_genes(){
    // struct gene_map_s {
    //     //
    //     unsigned long long genes_counter;
    //     //Gene start position (AUG)
    //     unsigned long long gene_start[MAX_GENES];
    //     //Gene stop position (UAA, UAG, UGA)
    //     unsigned long long gene_end[MAX_GENES];
    // };
    // gene_map_s gene_map;


// Algorithm:
    // ----------
    //     start_pos = lookup AUG in seq;
    // stop_pos = lookup(UAA or UAG or UGA) in seq + start_pos;
    // if (start_pos && stop_pos)   {
    //     gene_map.gene_start[gene_map.genes_counter) = start_pos;
    //     gene_map.gene_stop[gene_map.genes_counter) = stop_pos;
}

//////////////// Generating an amino acid chain (protein) 
void generating_amino_acid_chain(){

}

//////////////// Detecting probable mutation zones
void detecting_mutations(){

}

//////////////// Calculating the matching score of two sequences
void calculating_matching_score(){

}

//////////////// Hamming calculation
void hamming(int b1, int b2){
    int x = b1 ^ b2;
    return __builtin_popcount(x);
}
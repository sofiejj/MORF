%% MAP bases to numbers (alphabet index) - counterpart: map_num2base.m.

function protein_num_seq = map_base2num(prot_base_seq)

    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    
    [~, protein_num_seq] = ismember(prot_base_seq, Alphabet);

end

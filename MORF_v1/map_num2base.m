%% Map numbers to base letters using alphabet key - counterpart: 
% map_base2num.m.

function prot_base_seq = map_num2base(protein_num_seq)

    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    
    % Index directly into the Alphabet string:
    prot_base_seq = Alphabet(protein_num_seq);

end
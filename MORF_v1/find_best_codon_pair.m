function [bestA, bestB, bestMatches] = find_best_codon_pair(synCodonsNumA, synCodonsNumB)
%FIND_BEST_CODON_PAIR
% synCodonsNumA : 3 x NA numeric codons (alphabet key)
% synCodonsNumB : 3 x NB numeric codons (alphabet key)
%
% Returns:
%   bestA       : 3x1 numeric codon for seq A
%   bestB       : 3x1 numeric codon for seq B
%   bestMatches : number of matching bases (0..3)

    [~, NA] = size(synCodonsNumA);
    [~, NB] = size(synCodonsNumB);

    % Initialise with something valid
    bestMatches = -Inf;
    bestA       = synCodonsNumA(:,1);
    bestB       = synCodonsNumB(:,1);

    for i = 1:NA
        for j = 1:NB
            cA = synCodonsNumA(:,i);
            cB = synCodonsNumB(:,j);

            % Count how many positions match (0..3)
            matches = sum(cA == cB);

            if matches > bestMatches
                bestMatches = matches;
                bestA       = cA;
                bestB       = cB;
            end
        end
    end
end
%% Sliding Hamming Similarity: 
% determines percentage overlap of sequences A and B over a full
% sliding window.
% INPUT: full base seqs for proteins A and B.
% OUTPUT: best overlap of A and B in terms of percentage base match
% following full sliding window analysis (fraction as well as index window
% provided).
% The algorithm ALWAYS defines offset using: B(1) relative to A(1)
% N = number of Co-Encoding Regions (CERs) desired for output.
% L = minCERlength (in terms of number of nucleotides, INTEGER).
% bestOffset and bestFrac are 1xN numeric vectors.

function [bestOffset, bestFrac, A, B] = slidingHamming(seqA, seqB, N, L)
    
    A = map_base2num(seqA); % numeric mapping.
    B = map_base2num(seqB);

    A = A(:)'; 
    B = B(:)';

    nA = numel(A);
    nB = numel(B);

    if nB > nA % once l_min/L/minCERlength is integrated, this is no longer relevant.
        error('This version assumes seqB is not longer than seqA.');
    end
    % offsets = location of seqB START codon relative to seqA START.
    offsets = L - nB : nA - L; % SB: this is all wrong! I want L to be minimum OVERLAP length.
    %offsets   = 0 : (nA - nB); % THIS keeps the sliding seq B entirely within seq A.
    matchFrac = zeros(size(offsets));

    for k = 1:numel(offsets) % IDENTIFYING INDICES, iA and iB, WITHIN CERs given by offsets.
        d = offsets(k);
        % length iA and iB should always match.
        if d < 0 % if seqB STARTS BEFORE seqA

            iA = 1 : (d + nB);
            iB = (-d + 1) : nB;

        elseif d > (nA - nB) % if seqB ENDS AFTER seqA. SB 21/06/2026: this currently isn't quite right.

            iA = (1 + d) : nA; % The length of these should not go below L.
            iB = 1 : (nA - d);
        
        else % if seqB WITHIN seqA

            iA = (1 + d) : (d + nB);
            iB = 1 : nB;

        end


        % full overlap: B(1..nB) ↔ A(1+d .. d+nB)
        % iA = (1 + d) : (d + nB);
        % iB = 1:nB;
        len_CER(1,k) = length(iA); % Check that these do not go below L.
        len_CER(2,k) = length(iB); % track length iB for DEBUG.
        matchFrac(k) = sum(A(iA) == B(iB)) / len_CER(1,k); % checks for nucleotide matches/mismatches.
    end

    % What if N > length(matchFrac)?
    if N <= length(matchFrac)
        [bestFrac, idx] = maxk(matchFrac,N); % great that matchFrac is a list.
        bestOffset = offsets(idx);
    else 
        error('Number of CERs requested exceeds total number of CERs identified with the given parameters.\n');
    end
end


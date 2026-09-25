% MORF proof-of-concept implementation.
% INDEXING NOTE: if first element in A is index 1 for forward reading, then
% first element of A when reading backwards is -1.
% nt = nucleotides.
% Whichever seq starts first *matters*. This is considered the ref seq.
close all;
clear all; clc;

% MORF colour palette
purple   = [0.42, 0.16, 0.55];   % Original / initial
pink     = [0.91, 0.25, 0.55];   % Optimised
plum     = [0.25, 0.08, 0.32];   % Emphasis / best result
lavender = [0.88, 0.82, 0.91];   % Mismatch / secondary background
charcoal = [0.25, 0.25, 0.25];   % Text / neutral reference lines

%% Sliding Hamming Similarity

seqA = seq_chooser('GFP'); % GFP string
seqB = seq_chooser('dsRed'); % dsRed string

num_CERs = 10; % USER SPECIFIED (integer).
minCERlength = 339; % l_min: minimum allowed length of CER in terms of number of nucleotides. INTEGER. 339 is half length of current seqB (dsRed).

[bestOffset, bestFrac, num_seqA, num_seqB] = slidingHamming(seqA, seqB, num_CERs, minCERlength); % num_seqA, num_seqB = 'numeric' full protein sequences.
% length of respective sequences (in nt):
len_seqA = numel(num_seqA);
len_seqB = numel(num_seqB);
%% Multi-CER Optimisation
% Loop over CER identifiers (i.e. 1 : num_CERs):
%To keep per CER iteration: cer, num_CERs, bestOffset, bestFrac, num_seqA,
%num_seqB, len_seqA, len_seqB
for cer = 1 : num_CERs

   offset = bestOffset(cer);
    frac = bestFrac(cer);
    if offset < 0 % SB: seqB starts before seqA, seqA considered shifted.
        fprintf('SeqB starts before seqA.\n');
        ref_seq = num_seqB; 
        ref = 'B';
        RF_B = 0; % considered ref.
        A_CER_start = 1;
        B_CER_start = 1 - offset;
        A_CER_end = len_seqB + offset;
        B_CER_end = len_seqB;

        % Calculate rf of A:
        if rem(B_CER_start, 3)==1
            RF_A = 0;
            fprintf('SeqA in reading frame 0.\n');
        elseif rem(B_CER_start, 3)==2
            RF_A = 1;
            fprintf('SeqA in reading frame +1.\n');
        elseif rem(B_CER_start, 3)==0
            RF_A = 2;
            fprintf('SeqA in reading frame +2.\n');
        else
            error('Error in sequence B division by three. Please debug.\n');
        end

    elseif offset > (len_seqA-len_seqB)
        fprintf('SeqB ends outwith seqA.\n');
        ref_seq = num_seqA; 
        ref = 'A';
        RF_A = 0;
        A_CER_start = 1 + offset;
        B_CER_start = 1;
        A_CER_end = len_seqA;
        B_CER_end = len_seqA - offset; % End index of the overlapping region in sequence B.

        % Calculate rf of B:
        if rem(A_CER_start, 3)==1
            RF_B = 0;
            fprintf('SeqB in reading frame 0.\n');
        elseif rem(A_CER_start, 3)==2
            RF_B = 1;
            fprintf('SeqB in reading frame +1.\n');
        elseif rem(A_CER_start, 3)==0
            RF_B = 2;
            fprintf('SeqB in reading frame +2.\n');
        else
            error('Error in sequence A division by three. Please debug.\n');
        end

    else % CER entirely within seqA:
        fprintf('SeqB entirely within seqA.\n');
        ref_seq = num_seqA; 
        ref = 'A';
        RF_A = 0;
        A_CER_start = 1 + offset;
        B_CER_start = 1;
        A_CER_end = A_CER_start + len_seqB - 1; 
        B_CER_end = len_seqB;

        % Calculate rf of B:
        if rem(A_CER_start, 3)==1
            RF_B = 0; % first nts of B cods aligned with first nts of A cods.
            fprintf('SeqB in reading frame 0.\n');
        elseif rem(A_CER_start, 3)==2
            RF_B = 1; % first nts of B cods aligned with middle nts of A cods.
            fprintf('SeqB in reading frame +1.\n');
        elseif rem(A_CER_start, 3)==0
            RF_B = 2; % first nts of B cods aligned with last nts of A cods.
            fprintf('SeqB in reading frame +2.\n');
        else
            error('Error in sequence A division by three. Please debug.\n');
        end
    end

    % Check lengths - should always match:
    A_CER = length(A_CER_start : A_CER_end); % in terms of nts.
    B_CER = length(B_CER_start : B_CER_end);

    bestPercent = frac * 100;

    %% Basic Plots
    fig_start_idx = 4*cer-3;

    figure(fig_start_idx); % PLOT FULL SEQUENCES
    plot(1:len_seqA,num_seqA);hold on;
    plot(1+offset:offset+len_seqB,num_seqB);%plot(A_CER_start+offset:A_CER_start+offset+ln_seqB-1, num_seqB);
    legend('GFP (M62653.1)','DsRed'); % OFFSET in no. nucleotides where offset = +-1 is A(+-2), A(1)=A(-1).
    title(sprintf('Visual of CER %.f Location, OFFSET %.f, (%.2f%% Matching Bases) for GFP and DsRed (Sliding Hamming Similarity)',cer,offset,bestPercent)); % should start idx relative to seqA be provided?
    xlabel('Base Index');ylabel('Numeric Base Value From Alphabet Key');

    CER_window_A = num_seqA(A_CER_start:A_CER_end); 
    CER_window_B = num_seqB(B_CER_start:B_CER_end); % CER window in nts.
    
    if offset < 0 % 14/09: ref = 'B' and we need to check RF_A to place padding correctly.
        padding = zeros(1, RF_A);
        extend = size(padding,2);
        if extend == 1
            CER_window_A_extension = [CER_window_A num_seqA(A_CER_end+1)];
            CER_window_B_extension = [num_seqB(B_CER_start-1) CER_window_B];
        elseif extend == 2
            CER_window_A_extension = [CER_window_A num_seqA(A_CER_end+1) num_seqA(A_CER_end+2)];
            CER_window_B_extension = [num_seqB(B_CER_start-2) num_seqB(B_CER_start-1) CER_window_B];
        else % no extension (RF==0)
            CER_window_A_extension = CER_window_A;
            CER_window_B_extension = CER_window_B;
        end
        CER_window_A_padded = [padding CER_window_A];
        CER_window_B_padded = [CER_window_B padding];
    else %offset > 0 so ref='A' and need to check RF_B for padding:
        padding = zeros(1, RF_B);
        extend = size(padding,2);
        if extend == 1
            CER_window_A_extension = [num_seqA(A_CER_end-1) CER_window_A];
            CER_window_B_extension = [CER_window_B num_seqB(B_CER_start+1)];
        elseif extend == 2
            CER_window_A_extension = [num_seqA(A_CER_end-2) num_seqA(A_CER_end-1) CER_window_A];
            CER_window_B_extension = [CER_window_B num_seqB(B_CER_start+1) num_seqB(B_CER_start+2) ];
        else  % no extension
            CER_window_A_extension = CER_window_A;
            CER_window_B_extension = CER_window_B;
        end
        CER_window_A_padded = [CER_window_A padding];
        CER_window_B_padded = [padding CER_window_B];
    end

    optimisation_codA_setOne = reshape(CER_window_A_padded,3,[]); 
    optimisation_codB_setOne = reshape(CER_window_B_extension,3,[]);
    optimisation_codA_setTwo = reshape(CER_window_A_extension,3,[]);
    optimisation_codB_setTwo = reshape(CER_window_B_padded,3,[]);
    % CER_window_A_triplet = reshape(CER_window_A_padded,3,[]); % SB 27/07: codon fold inc. PADDING.
    % CER_window_B_triplet = reshape(CER_window_B_padded,3,[]);
    optimisation_array = [optimisation_codA_setOne; optimisation_codB_setOne; optimisation_codA_setTwo; optimisation_codB_setTwo]; % 'setOne' for optimising seqB, 'setTwo' for optimising seqA.
    % CER_cod_windows = [CER_window_A_triplet; CER_window_B_triplet]; % the CER folded into TRIPLETS, inc. padding. 
    
    figure(fig_start_idx+1); % PLOT CLOSEUP. 
    plot(A_CER_start:A_CER_end,CER_window_A);
    hold on;
    plot(A_CER_start:A_CER_end, CER_window_B);
    title(sprintf('Closeup of Co-Encoding Region %.f, OFFSET %.f',cer,offset));xlabel('Base Index');ylabel('Numeric Base Value from Alphabet Key');legend('GFP (M62653.1)','dsRed');
    
    %% CALCULATE nucleotide matches and mismatches:
    nt_match_array = ( CER_window_A == CER_window_B ); % output: logical array, 1 = NUCLEOTIDE MATCH location, 0 = NUCLEOTIDE MISMATCH location. SB 29/07: CORRECT.
    
    nt_match_array = nt_match_array(:);  % ensure column vector
    N = numel(nt_match_array); 
    
    figure(fig_start_idx+2); % PLOT MATCH ARRAY as DNA 'barcode'. %SB 21/06/2026: add offset.
    for k = 1:N
        if nt_match_array(k)
            rectangle('Position', [k 0 1 1], 'FaceColor', purple, 'EdgeColor', 'none'); % green match
        else
            rectangle('Position', [k 0 1 1], 'FaceColor', lavender, 'EdgeColor', 'none'); % grey mismatch
        end
    end
    xlim([0 N]);
    ylim([0 1]);
    axis off;
    title(sprintf('CER %.f Initial (%.2f%% matching)', cer, mean(nt_match_array)*100));

%% CODON match array (inc. padding & extensions, see lines 127-139):

    triplet_match_array(1:3,:) = ( optimisation_array(1:3,:) == optimisation_array(4:6,:) ); % SB 05/08: for optimising seqB.
    triplet_match_array(4:6,:) = ( optimisation_array(7:9,:) == optimisation_array(10:12,:) ); % SB 05/08: for optimising seqA.
    
    optimisationLen = size(triplet_match_array,2);
    number_of_seqs = size(triplet_match_array,1) / 3;

    %% Identify codon pairs requiring substitution:
    codonTable = build_numeric_codon_table(); % table of subs.

    for n = 1 : number_of_seqs 
        top_row = 6*n - 5; % These numbers match with a row count of 12 in optimisation_array.
        top_match_array = 3*n - 2;
        bottom_row = 6*n - 3;
        bottom_match_array = 3*n;
        codon_array_A_window = optimisation_array(top_row:bottom_row,:);
        codon_array_B_window = optimisation_array(top_row+3:bottom_row+3,:);
    
        for i = 1 : optimisationLen
        
    
            if sum(triplet_match_array(top_match_array:bottom_match_array,i))~=3 % when sum~=3, codon MISMATCH
                fprintf('Codon MISMATCH at optimisation index: %.f\n', i); % it is technically a NUCLEOTIDE mismatch - only a codon mismatch if no RF shift.
                codA = codon_array_A_window(:,i); 
                codB = codon_array_B_window(:,i);
                
                if n==2 % GET codon A 'synonyms':
                    idxA = codonTable.codonIndexNum(codonTable.makeKey(codA));
                    synIdxA       = codonTable.synonymIdx{idxA};
                    synCodonsNumA = codonTable.codonsNumeric(:, synIdxA);

                    synCodonsNumB = codB; % SB 05/08: loop over CER A, hold B fixed.
                else % n==1
                    % GET codon B 'synonyms':
                    idxB = codonTable.codonIndexNum(codonTable.makeKey(codB));
                    synIdxB       = codonTable.synonymIdx{idxB};
                    synCodonsNumB = codonTable.codonsNumeric(:, synIdxB);

                    synCodonsNumA = codA; % SB 05/08: loop over CER B, hold A fixed.
                end
        
                % Possible outcomes: - perfect sub, better sub, no sub.
                [bestA, bestB, bestMatches] = find_best_codon_pair(synCodonsNumA, synCodonsNumB);
        
                % Compare with original match quality:
                origMatches = sum(codA == codB);
            
                if bestMatches > origMatches
                    % Apply improved codons
                    optimisation_array(top_row:bottom_row,i) = bestA; % SB 05/08: when looping over B, A should remain unchanged, and vice versa.
                    optimisation_array(top_row+3:bottom_row+3,i) = bestB; %14/09: for n=1, A fixed, B subbed (rows 4:6 in opt array).

                    fprintf('New codon pair at optimisation index %.f\n',i); 
                else
                    % No improvement – log for later reporting
                    % noSolutionIdx(end+1) = i;
                    fprintf('No better codon pair found at optimisation index %.f\n',i);
                end
    
            else % codon MATCH
                fprintf('Codon MATCH at optimisation index: %.f\n', i);
                fprintf('Skipping to next codon pair...\n');
            end

        end
        % 14/09: need to carry through 'optimised B' for A analysis.
        % flatten opt_B into nts (1st iteration):
        if n == 1
            B_extend_opt_nts = reshape(optimisation_array(4:6,:), 1, []); 
            if offset < 0 % ref = 'B'
                B_padded_new_nts = [B_extend_opt_nts(1+extend:end) zeros(1,extend)]; % carry updates forward for A optimisation.
            else % ref = 'A'
                B_padded_new_nts = [zeros(1,extend) B_extend_opt_nts(1:end-extend)];
            end
    
            optimisation_array(10:12,:) = reshape(B_padded_new_nts, 3, []);  % B updates (1st iteration) carry forward to A updates (2nd iteration).
        end

        triplet_match_array(1:3,:) = ( optimisation_array(1:3,:) == optimisation_array(4:6,:) ); 
        triplet_match_array(4:6,:) = ( optimisation_array(7:9,:) == optimisation_array(10:12,:) ); 

    end

    % final altered seqs in: opt_array(7:9) and opt_array(10:12)
    % -extension and -padding, respectively. Most up to date triplet_match_array is
    % triplet_match_array(4:6,:).
    
    %% Calculate new match percentage:
    new_opt_window_A = reshape( optimisation_array(7:9,:), 1, [] ); 
    new_opt_window_B = reshape( optimisation_array(10:12,:), 1, [] );
    if offset < 0 % ref = 'B'
        new_CER_window_A(cer,:) = new_opt_window_A(1:end-extend); 
        new_CER_window_B(cer,:) = new_opt_window_B(1:end-extend);
    else % ref = 'A'
        new_CER_window_A(cer,:) = new_opt_window_A(1+extend:end);
        new_CER_window_B(cer,:) = new_opt_window_B(1+extend:end);
    end

    % Insert optimised CERs into full (numeric) sequences:
    opt_num_seqA(cer,:) = [num_seqA(1:A_CER_start-1) new_CER_window_A(cer,:) num_seqA(A_CER_end+1:end)]; %16/09: struct here?
    opt_num_seqB(cer,:) = [num_seqB(1:B_CER_start-1) new_CER_window_B(cer,:) num_seqB(B_CER_end+1:end)];
    new_matches_nt = sum(new_CER_window_A(cer,:) == new_CER_window_B(cer,:));   % scalar
    total   = numel(new_CER_window_A(cer,:));                  % = 3 * L % SB : is this not simply cerLen? Not exactly: it is [cerLen * 3], i.e. length of nt seq.
    
    pct(cer) = 100 * (new_matches_nt / total); % improved percentage matches for CER cer.
    
    fprintf('New base matches percentage after substitutions: %.2f%%\n',pct(cer));
    
    %% Produce new match BARCODE:
    match_array_new = ( new_CER_window_A(cer,:) == new_CER_window_B(cer,:) );
    % match_array_bases = reshape(match_array_new,1,[]); % Reshape from codon triplets to nucleotide sequence for barcode.
    % N = numel(nt_match_array); 
    
    figure(fig_start_idx+3); % PLOT MATCH ARRAY as DNA 'barcode'. 
    for k = 1:N
        if match_array_new(k)
            rectangle('Position', [k 0 1 1], 'FaceColor', pink, 'EdgeColor', 'none'); % green match
        else
            rectangle('Position', [k 0 1 1], 'FaceColor', lavender, 'EdgeColor', 'none'); % grey mismatch
        end
    end
    xlim([0 N]);
    ylim([0 1]);
    axis off;
    title(sprintf('CER %.f Post-Optimisation (%.2f%% matching)', cer, pct(cer)));

    clearvars -except cer pct num_CERs bestOffset bestFrac num_seqA num_seqB len_seqA len_seqB seqA seqB opt_num_seqA opt_num_seqB purple pink plum lavender charcoal% variable clear between CER iterations to start over.

end

% Identify best match score:
[bestPct,cer_idx] = max(pct);
fprintf('Best CER = %.f, with match score %.2f%%, increase from match score %.2f%%.\n',cer_idx,bestPct,bestFrac(cer_idx)*100);

% Map back to bases (from numeric):
opt_seqA = map_num2base(opt_num_seqA);
opt_seqB = map_num2base(opt_num_seqB);

% Check for AA preservation:
aa_dsRed_original  = translate_nt(seqB);
aa_GFP_original = translate_nt(seqA);
aa_dsRed_optimised = translate_nt(opt_seqB);
aa_GFP_optimised = translate_nt(opt_seqA);

dsRed_preserved = all(aa_dsRed_optimised == aa_dsRed_original, 2);
GFP_preserved = all(aa_GFP_optimised == aa_GFP_original, 2);

% Calculate means:
mean_init_pct = mean(bestFrac*100);
mean_opt_pct = mean(pct);

% Generate results fig:
figure(100); 
plot(1:num_CERs,bestFrac*100,'o','MarkerSize',12,'MarkerFaceColor',purple,'MarkerEdgeColor','k','LineWidth',2);hold on;
plot(1:num_CERs,pct,'o','MarkerSize',12, 'MarkerFaceColor',pink,'MarkerEdgeColor','k','LineWidth',2);hold on;
xlim([1,num_CERs]);
yline(bestPct,'--',sprintf('BEST: %.2f%%',bestPct),'LineWidth',3,'Color',plum);
yline(mean_init_pct,'--',sprintf('Mean Initial Match Score %.2f%%',mean_init_pct),'Color',charcoal)
yline(mean_opt_pct,'--',sprintf('Mean Post-Optimisatin Match Score %.2f%%',mean_opt_pct),'Color',charcoal)
xline(cer_idx,'--',sprintf('BEST CER %d',cer_idx),'LineWidth',3,'Color',plum);
grid on;
title('SeqA = GFP, seqB = dsRed, l_{min} = length(seqB)'); % need more info about the proteins.
xlabel('Co-Encoding Region');ylabel('Nucleotide Match Score (%)');
legend('Initial Match Score','Post-Optimisation Match Score','Best Final Score');

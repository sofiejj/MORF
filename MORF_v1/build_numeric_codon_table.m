function codonTable = build_numeric_codon_table()
%BUILD_CODON_TABLE_ALPHABETKEY
% Construct a codon substitution table using the global Alphabet key
% (A=1, B=2, ..., Z=26), e.g. 'ATC' -> [1; 20; 3].
%
% Fields:
%   .Alphabet       - 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
%   .codonsNumeric  - 3xN numeric codons (A=1,...,Z=26)
%   .codonIndexNum  - Map: '1_20_3' -> column index
%   .synonymIdx     - synonymIdx{i} = indices of synonymous codons
%   .aa             - amino acid per codon (optional)
%   .makeKey        - helper: triplet (3x1) -> key string

    %----------------------------------------
    % Alphabet and encoders (match map_base2num)
    %----------------------------------------
    codonTable.Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    encodeBase = @(b) find(codonTable.Alphabet == b, 1, 'first');
    decodeBase = @(n) codonTable.Alphabet(n);

    %----------------------------------------
    % Letter codon map (AA -> codon strings)
    %----------------------------------------
    codonMap = containers.Map;

    codonMap('A') = {'GCT','GCC','GCA','GCG'};
    codonMap('C') = {'TGT','TGC'};
    codonMap('D') = {'GAT','GAC'};
    codonMap('E') = {'GAA','GAG'};
    codonMap('F') = {'TTT','TTC'};
    codonMap('G') = {'GGT','GGC','GGA','GGG'};
    codonMap('H') = {'CAT','CAC'};
    codonMap('I') = {'ATT','ATC','ATA'};
    codonMap('K') = {'AAA','AAG'};
    codonMap('L') = {'TTA','TTG','CTT','CTC','CTA','CTG'};
    codonMap('M') = {'ATG'};%{'ATG','GTG','TTG'};
    codonMap('N') = {'AAT','AAC'};
    codonMap('P') = {'CCT','CCC','CCA','CCG'};
    codonMap('Q') = {'CAA','CAG'};
    codonMap('R') = {'CGT','CGC','CGA','CGG','AGA','AGG'};
    codonMap('S') = {'TCT','TCC','TCA','TCG','AGT','AGC'};
    codonMap('T') = {'ACT','ACC','ACA','ACG'};
    codonMap('V') = {'GTT','GTC','GTA','GTG'};
    codonMap('W') = {'TGG'};
    codonMap('Y') = {'TAT','TAC'};
    codonMap('*') = {'TAA','TAG','TGA'};

    %----------------------------------------
    % Flatten into lists of codons + AA labels
    %----------------------------------------
    AAs     = keys(codonMap);
    codList = {};
    aaList  = {};

    for i = 1:numel(AAs)
        aa   = AAs{i};
        cods = codonMap(aa);
        codList = [codList, cods];                 %#ok<AGROW>
        aaList  = [aaList, repmat({aa},1,numel(cods))]; %#ok<AGROW>
    end

    N = numel(codList);
    codonsNumeric = zeros(3, N);

    %----------------------------------------
    % Convert each codon string to numeric [1..26]
    %----------------------------------------
    for k = 1:N
        trip = codList{k};            % e.g. 'ATC'
        codonsNumeric(:,k) = arrayfun(encodeBase, trip(:));
        % 'A'->1, 'T'->20, 'C'->3  => [1;20;3]
    end

    %----------------------------------------
    % Build numeric-key map for lookup
    %
    % Use a safe key format with separators, e.g. '1_20_3'
    %----------------------------------------
    codonIndexNum = containers.Map('KeyType','char','ValueType','int32');

    makeKey = @(trip) sprintf('%d_%d_%d', trip(1), trip(2), trip(3));

    for k = 1:N
        trip = codonsNumeric(:,k);
        key  = makeKey(trip);
        codonIndexNum(key) = k;
    end

    %----------------------------------------
    % Build synonym index list (same AA = synonyms)
    %----------------------------------------
    synonymIdx = cell(1, N);

    for i = 1:numel(AAs)
        aa  = AAs{i};
        idx = find(strcmp(aaList, aa));   % all codons for this AA
        for j = 1:numel(idx)
            synonymIdx{idx(j)} = idx;
        end
    end

    %----------------------------------------
    % Pack into struct
    %----------------------------------------
    codonTable.codonsNumeric = codonsNumeric;
    codonTable.codonIndexNum = codonIndexNum;
    codonTable.synonymIdx    = synonymIdx;
    codonTable.aa            = aaList;

    % Helpers (if you ever need letters again)
    codonTable.encodeBase = encodeBase;
    codonTable.decodeBase = decodeBase;
    codonTable.makeKey    = makeKey;

end
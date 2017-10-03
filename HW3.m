%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

match=2;
mismatch=-1;
ofdiag=ones(4)-eye(4);
s=match*eye(4)+mismatch*ofdiag;
seq1='GTAATCC';
seq2='GTATCCG';
[score,align,start]=swalign(seq1,seq2,'Alphabet','nt','ScoringMatrix',s,'Gapopen',1,'Showscore',true);
align

%%
%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

ERK1=getgenbank('NM_002746');
ERK2=getgenbank('NM_002745');
[score1, align, start1]=swalign(ERK1.Sequence,ERK2.Sequence,'Alphabet','nt','Showscore',true);
showalignment(align);
%%Bingyan Wu: 1053/1506 (70%) alignment

%%
% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

ERK1_AA=getgenpept(ERK1.CDS.protein_id);
ERK2_AA=getgenpept(ERK2.CDS.protein_id);

[pscore, palign, pstart]=swalign(ERK1_AA.Sequence,ERK2_AA.Sequence);
showalignment(palign);
%%Bingyan Wu: 305/346 (88%) alignment

%%
% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

mouseERK1=getgenbank('X64605.1');
mouseERK1_AA=getgenpept(mouseERK1.CDS.protein_id);
mouseERK2=getgenbank('D10939.1');
mouseERK2_AA=getgenpept(mouseERK2.CDS.protein_id);

[score1, align1, start1]=swalign(ERK1.Sequence, mouseERK1.Sequence);
showalignment(align1);
[score2, align2, start2]=swalign(ERK2.Sequence, mouseERK2.Sequence);
showalignment(align2);


[score_AA1, align_AA1, start_AA1]=swalign(ERK1_AA.Sequence, mouseERK1_AA.Sequence);
showalignment(align_AA1);
[score_AA2, align_AA2, start_AA2]=swalign(ERK2_AA.Sequence, mouseERK2_AA.Sequence);
showalignment(align_AA2);

%%Bingyan Wu: for ERK1, the gene sequence is 90% similar, ERK2 is 75%
%%similar; the amino acid for ERK1 is 97% similar and that for ERK2 is 99%
%%similar

%%
%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

topNhits('NM_002746',5)

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

ClosestMatch('NM_002746')

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

ClosestMatch('NM_005228') % Homo sapiens EGFR
% the closest human gene: 'NM_005228'    'Homo sapiens epidermal growth factor receptor (E?'
% the closest non-human gene: 'XM_519102'    'PREDICTED: Pan troglodytes epidermal growth fact?'

ClosestMatch('L37053') %Rh blood group D antigen [ Gorilla gorilla (western gorilla) ]
% the closest human gene: 'NM_016124'    'Homo sapiens Rh blood group D antigen (RHD), tra?'
% the closest non-human gene: 'NM_001279597 XM_004025193'    'Gorilla gorilla Rh blood group D antigen (RHD), ?'

%% it makes sense that the human ortholog of Rh group D antigen was returned and the closest non-human gene was found in Gorilla



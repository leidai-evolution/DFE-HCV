%% pre-processing of HCV NS5A sequence alignment
%updated: 06/23/2016, LD
%sequence data: HCV LANL database
%genotype: all or 2; NS5A region; aligned
%manually examined by MEGA7

%% loading
fastaFile='./data_LANL/HCV_genotypeall_NS5A.fasta';
[label,S]=fastaread(fastaFile);
n_seq = length(label);

%% find reference: H77 genome
label_ref='NC_004102';
index_ref=find(~cellfun(@isempty,strfind(label,label_ref)));

%idenfity gaps in reference
bp_start=69;
bp_end=328;
% gap1=203-bp_start+1;
% gap2=314-bp_start+1;
gap_pos=strfind(S{index_ref}(bp_start:bp_end),'-');
seq_pos=setdiff(1:bp_end-bp_start+1,gap_pos); %remove gaps in reference

%%
seq_temp=S{index_ref}(bp_start:bp_end);
seq_nt_ref=seq_temp(seq_pos)
seq_aa_ref=nt2aa(seq_nt_ref,'ACGTOnly','false')

%% exclude sequences that cause gaps in the reference sequence
tic;
count=0;
for i=1:n_seq
    seq_temp = S{i}(bp_start:bp_end);
    if seq_temp(gap_pos(1))~='-' || seq_temp(gap_pos(2))~='-'  %note: two sequences identified
        index_include(i)=0;
    else
        index_include(i)=1;
        count=count+1;
        seq_nt{count}=seq_temp(seq_pos);
        seq_nt{count}=strrep(seq_nt{count},'-','N');  % setting all gaps equal to unknown base, N 
                                                %  valid since ref seq is continuous and so these gaps are either deletions or unknowns
    end
end
toc;

%% translate to amino acid sequence
%sequences w/ gaps or ambiguous nt may lead to 'X' after translation to aa
tic;
for i=1:length(seq_nt)
    seq_aa{i} = nt2aa(seq_nt{i},'ACGTOnly','false');      % translating accepting non ACGT bases
    seq_int{i} = aa2int(seq_aa{i});                       % coding aa as int
    seq_int{i}(seq_int{i} > 21) = 21;                     % setting all ambiguous residues to 21
end
toc;

%% transform to matrix format: easier for downstream analysis (entropy)
seq_int_mat=zeros(length(seq_int),86);
for i=1:length(seq_int)
    seq_int_mat(i,:)=cell2mat(seq_int(i));
end
%exclude sequences w/ ambiguous residues (aa int=21)
count=0;
for i=1:size(seq_int_mat,1)
    if isempty(find(seq_int_mat(i,:)==21))
        count=count+1;
        seq_int_mat_select(count,:)=seq_int_mat(i,:);
    end
end

%% calculate entropy
N=20;
entropy=entropy_msa(seq_int_mat_select,N);

%%
figure(1);
site=18:103;
stem(site,entropy);
xlabel('Site');
ylabel('Entropy (bits)');
set(gca,'xlim',[18 103],'ylim',[0 max(entropy)]);
% title(strcat(seq_label,' lineage'));

%% save processed data
%downstream analysis: compare entropy to fitness cost
% save('./data_LANL/genotypeall_NS5AD1.mat','seq_int_mat_select','entropy');



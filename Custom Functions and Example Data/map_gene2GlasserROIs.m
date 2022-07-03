function [ROI_expression_mat] = map_gene2GlasserROIs(expression_csv, samples_csv, ontology_csv, probes_csv)
% map_gene2GlasserROIs: This function assigns microarray expression values
% in each sample from an AHBA subject's 'expression.csv' file to the
% nearest Glasser parcellation ROIs, using Euclidean distance measured in
% MNI coordinate space.

%To run this function, the following dependencies must be added to the Matlab search path:
%    - CONN toolbox: available for download at <https://web.conn-toolbox.org/>
%    - Midnight Scan Club Codebase: available for download at <https://github.com/MidnightScanClub/MSCcodebase>
%    - Richiardi_Data_File_S2.csv file available for download in supplement of
%      <Richiardi et al, 'Correlated gene expression supports synchronous activity in brain networks', 
%      Science, 2015. DOI: 10.1126/science.1255905> 
%    - MNI template surfaces: 
%        - lh.midthickness.32k_fs_LR.surf.gii file attached in Github repository
%        - rh.midthickness.32k_fs_LR.surf.gii file attached in Github repository
%    - Glasser_SubCort.dtseries.nii file attached in Github repository

% Inputs should be filepaths to expression.csv, samples.csv, ontology.csv,
% and probes_csv files in a downloaded AHBA subject data folder. 

%Eg: [subj9861] = samples2GlasserROIs('normalized_microarray_donor9861/MicroarrayExpression.csv','normalized_microarray_donor9861/SampleAnnot.csv','normalized_microarray_donor9861/Ontology.csv','normalized_microarray_donor9861/Probes.csv');

%% Upload necessary files

%Use Glasser cifti files to define MNI coordinates of 91,282 voxels in the
%Glasser parcellation:
P = ft_read_cifti_mod('Glasser_SubCort.dtseries.nii');

coordsL = gifti('lh.midthickness.32k_fs_LR.surf.gii');
coordsL = coordsL.vertices; % XYZ coordinates for LH
coordsR = gifti('rh.midthickness.32k_fs_LR.surf.gii');
coordsR = coordsR.vertices; % XYZ coordinates for RH
coords_surf=[coordsL;coordsR]; 
surf_indices_incifti = P.brainstructure > 0;
brain_structures_incifti = P.brainstructure(surf_indices_incifti);
surf_indices_incifti_cx = surf_indices_incifti(1:size(coords_surf,1));
coords_cifti_surf = coords_surf(surf_indices_incifti_cx,:);

clear surf_indices_incifti surf_indices_incifti_cx coordsL coordsR coords_surf

%combine into one coordinates array matching length of 'P.data' for cifti
%files
all_coords_mni = coords_cifti_surf;

clear coords_cifti_surf

%Load gene expression files
expression_array = table2array(readtable(expression_csv,'ReadVariableNames',0,'Delimiter',','));
samples = readtable(samples_csv,'ReadVariableNames',1,'Delimiter',',');
probes = readtable(probes_csv,'ReadVariableNames',1,'Delimiter',',');
column_identifiers = readtable(ontology_csv,'ReadVariableNames',1,'Delimiter',',');

%reannotate probes
reannot=readtable('Richiardi_Data_File_S2.csv');
for i = 1:size(reannot,1)
    probes.entrez_id(strcmp(probes.probe_name, reannot.probe_id{i})) = reannot.Entrez_id(i);
end

%for each gene, take the mean of all probes measuring expression of that
%gene:

gene_eID_list = unique(probes.entrez_id); 
gene_eID_list = gene_eID_list(~isnan(gene_eID_list));

meaned_expression_array = zeros(length(gene_eID_list),size(expression_array,2));

for i = 1:length(gene_eID_list)
    meaned_expression_array(i,:) = mean(expression_array(probes.entrez_id == gene_eID_list(i),:),1);
end

%Replace first column with corresponding entrez gene ID's for later indexing:
meaned_expression_array(:,1) = gene_eID_list;

%% Assign the most proximate Glasser ROI to each sample site using MNI coordinates
%Locate the nearest voxel to each sample in MNI space
%and index P.data to locate the Glasser ROI containing that voxel

%Do this separately for each hemisphere, cort/subcort, and cerebellum
%Preallocate structures for each anatomical designation:

CortexL.coords = all_coords_mni(brain_structures_incifti == 1,:);
CortexR.coords = all_coords_mni(brain_structures_incifti == 2,:);

CortexL.ROIs = P.data(brain_structures_incifti == 1);
CortexR.ROIs = P.data(brain_structures_incifti == 2);

CortexL.samples = zeros(size(CortexL.coords,1), size(samples,1));
CortexR.samples = zeros(size(CortexR.coords,1), size(samples,1));

CortexL.distances = zeros(size(CortexL.coords,1), size(samples,1));
CortexR.distances = zeros(size(CortexR.coords,1), size(samples,1));

%specify maximum distance threshold to define sphere around ABA microarray 
%sample sites, and only accept cifti coordinates within that distance
%threshold:

d_max = 2;
d_max_subcx = 2;

for i = 1:size(samples,1)
    %populate the CortexL.samples array with sample #'s for which that
    %cifti voxel lies within the d_max distance threshold sphere
    %surrounding that sample site. 
    
    %also populate the CortexL.distances array with corresponding Euclidean
    %distance for each sample #. This distance will later be used to weight
    %the mean expression value for Glasser voxels with >1 assigned sample
    %site. 
        
    %left hemisphere (ontology # 4008, L)
    if contains(column_identifiers.structure_id_path{column_identifiers.id == samples.structure_id(i)}, '4008') & column_identifiers.hemisphere{column_identifiers.id == samples.structure_id(i)} == 'L'
        %explanation of below: create a logical vector of
        %length(CortexL.coords) using pdist2, which identifies rows (voxels)
        %in CortexL.samples within max_d mm of the sample site. Index
        %CortexL.samples by that logical vector, and set the value = to i,
        %the sample #. Do the same to populate the CortexL.distances array,
        %except set the value = to the calculated 'pdist2' + 0.001. (add
        %0.001 so that a true distance of 0 can be differentiated from
        %the zeros assigned to empty cells)
        CortexL.samples(pdist2(CortexL.coords(:,:),[samples.mni_x(i) samples.mni_y(i) samples.mni_z(i)]) < d_max,i) = i;
        CortexL.distances(pdist2(CortexL.coords(:,:),[samples.mni_x(i) samples.mni_y(i) samples.mni_z(i)]) < d_max,i) = pdist2(CortexL.coords(pdist2(CortexL.coords(:,:),[samples.mni_x(i) samples.mni_y(i) samples.mni_z(i)]) < d_max,:),[samples.mni_x(i) samples.mni_y(i) samples.mni_z(i)]) + 0.001;
    
    %do the same for right hemisphere (ontology # 4008, R)
    elseif contains(column_identifiers.structure_id_path{column_identifiers.id == samples.structure_id(i)}, '4008') & column_identifiers.hemisphere{column_identifiers.id == samples.structure_id(i)} == 'R'
        CortexR.samples(pdist2(CortexR.coords(:,:),[samples.mni_x(i) samples.mni_y(i) samples.mni_z(i)]) < d_max,i) = i;
        CortexR.distances(pdist2(CortexR.coords(:,:),[samples.mni_x(i) samples.mni_y(i) samples.mni_z(i)]) < d_max,i) = pdist2(CortexR.coords(pdist2(CortexR.coords(:,:),[samples.mni_x(i) samples.mni_y(i) samples.mni_z(i)]) < d_max,:),[samples.mni_x(i) samples.mni_y(i) samples.mni_z(i)]) + 0.001;

    end
        
end

%preallocate ROI by expression matrices (only for ROIs with assigned sample
%sites, so as to accomodate Matlab matrix size limit. To do this, create a 
%'non-zero index' for each anatomical subdivision, which will be used later 
%to assign appropriate Glasser ROI #'s):

CortexL.nonzero_idx = find(sum(CortexL.samples,2)>0);
CortexL.expression = zeros(size(CortexL.nonzero_idx,1),size(meaned_expression_array,1));

CortexR.nonzero_idx = find(sum(CortexR.samples,2)>0);
CortexR.expression = zeros(size(CortexR.nonzero_idx,1),size(meaned_expression_array,1));

for i=1:size(CortexL.expression,1)
    
    %Populate the L cortex expression array with the mean
    %of expression values indexed by sample# from the imported 'expression array'
    CortexL.expression(i,:) = mean(meaned_expression_array(:,1+nonzeros(CortexL.samples(CortexL.nonzero_idx(i),:))),2)';
    
end

%Do the same for R cortex, subcortex, cerebellum: 

for i=1:size(CortexR.expression,1)
    CortexR.expression(i,:) = mean(meaned_expression_array(:,1+nonzeros(CortexR.samples(CortexR.nonzero_idx(i),:))),2)';
end


%add in a first column assigning Glasser ROIs (1 to 379) to expression
%values (in rows) in each expression table:

CortexL.expression = horzcat(CortexL.ROIs(CortexL.nonzero_idx), CortexL.expression);
CortexR.expression = horzcat(CortexR.ROIs(CortexR.nonzero_idx), CortexR.expression);

%concatenate into one array containing expression values with assigned 
%Glasser ROIs for the whole brain: 
Whole_Brain_expression = vertcat(CortexL.expression, CortexR.expression);


%% Populate a ROI-by-gene expression matrix with expression data 
%for each ROI, indexed from Sample_ROIs_mat

%preallocate
ROI_expression_mat = zeros(360,size(meaned_expression_array,1));

%populate this ROI_expression_mat with gene expression values for 
%corresponding ROIs, indexed from column 1 of Whole_Brain_expression:

for i = 1:360
    ROI_expression_mat(i,:) = mean(Whole_Brain_expression(Whole_Brain_expression(:,1)==i,2:end), 1);
end

ROI_expression_mat = vertcat(transpose(gene_eID_list), ROI_expression_mat);

%flip R and L hemisphere to put L hemisphere first
ROI_expression_mat = ROI_expression_mat([1 182:361 2:181],:);

    
    
   


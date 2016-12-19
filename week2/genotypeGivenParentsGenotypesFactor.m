function genotypeFactor = genotypeGivenParentsGenotypesFactor(numAlleles, genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo)
% This function computes a factor representing the CPD for the genotype of
% a child given the parents' genotypes.

% THE VARIABLE TO THE LEFT OF THE CONDITIONING BAR MUST BE THE FIRST
% VARIABLE IN THE .var FIELD FOR GRADING PURPOSES

% When writing this function, make sure to consider all possible genotypes 
% from both parents and all possible genotypes for the child.

% Input:
%   numAlleles: int that is the number of alleles
%   genotypeVarChild: Variable number corresponding to the variable for the
%   child's genotype (goes in the .var part of the factor)
%   genotypeVarParentOne: Variable number corresponding to the variable for
%   the first parent's genotype (goes in the .var part of the factor)
%   genotypeVarParentTwo: Variable number corresponding to the variable for
%   the second parent's genotype (goes in the .var part of the factor)
%
% Output:
%   genotypeFactor: Factor in which val is probability of the child having 
%   each genotype (note that this is the FULL CPD with no evidence 
%   observed)

% The number of genotypes is (number of alleles choose 2) + number of 
% alleles -- need to add number of alleles at the end to account for homozygotes

genotypeFactor = struct('var', [], 'card', [], 'val', []);

% Each allele has an ID.  Each genotype also has an ID.  We need allele and
% genotype IDs so that we know what genotype and alleles correspond to each
% probability in the .val part of the factor.  For example, the first entry
% in .val corresponds to the probability of having the genotype with
% genotype ID 1, which consists of having two copies of the allele with
% allele ID 1, given that both parents also have the genotype with genotype
% ID 1.  There is a mapping from a pair of allele IDs to genotype IDs and 
% from genotype IDs to a pair of allele IDs below; we compute this mapping 
% using generateAlleleGenotypeMappers(numAlleles). (A genotype consists of 
% 2 alleles.)

[allelesToGenotypes, genotypesToAlleles] = generateAlleleGenotypeMappers(numAlleles);

% One or both of these matrices might be useful.
%
%   1.  allelesToGenotypes: n x n matrix that maps pairs of allele IDs to 
%   genotype IDs, where n is the number of alleles -- if 
%   allelesToGenotypes(i, j) = k, then the genotype with ID k comprises of 
%   the alleles with IDs i and j
%
%   2.  genotypesToAlleles: m x 2 matrix of allele IDs, where m is the 
%   number of genotypes -- if genotypesToAlleles(k, :) = [i, j], then the 
%   genotype with ID k is comprised of the allele with ID i and the allele 
%   with ID j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INSERT YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Fill in genotypeFactor.var.  This should be a 1-D row vector.
genotypeFactor.var = [genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo]
% Fill in genotypeFactor.card.  This should be a 1-D row vector.
numberOfGene = numAlleles*(numAlleles-1)/2 + numAlleles
genotypeFactor.card = [numberOfGene,numberOfGene,numberOfGene]
genotypeFactor.val = zeros(1, prod(genotypeFactor.card));
% Replace the zeros in genotypeFactor.val with the correct values.
assignments = IndexToAssignment(1:prod(genotypeFactor.card), genotypeFactor.card);
% pick columns for factors in A
indexA = AssignmentToIndex(assignments(:, 2), numberOfGene);
indexB = AssignmentToIndex(assignments(:, 3), numberOfGene);
indexC = AssignmentToIndex(assignments(:, 1), numberOfGene);

for (i = 1:length(genotypeFactor.val))
    child_allele = genotypesToAlleles(indexC(i),:);
    par1_allele = genotypesToAlleles(indexA(i),:);
    par2_allele = genotypesToAlleles(indexB(i),:);
    % two allele of child are same
    if (child_allele(1) == child_allele(2))
        al = child_allele(1)
        result = 1
        if (~ismember(al, par1_allele)|~ismember(al, par2_allele))
            genotypeFactor.val(i) = 0;
        else
            if(par1_allele(1) ~= par1_allele(2))
                result = 0.5;
            end
            if (par2_allele(1) ~= par2_allele(2))
                result = result* 0.5;
            end
            genotypeFactor.val(i) = result
        end
        % two allele of child different and both of them contained in each
        % parent
    elseif (length(intersect(par1_allele, child_allele)) == 2 & length(intersect(par2_allele, child_allele)) == 2 )
        genotypeFactor.val(i) = 0.5;
    else % one allele must be get from one and only one of the parent
        c1 = child_allele(1)
        c2 = child_allele(2)
        % c1 match par1, c2 match par2
        if (ismember(c1, par1_allele) && ismember(c2, par2_allele))
           result = 1
           if (par1_allele(1)~=c1 | par1_allele(2)~=c1)
               result = 0.5
           end
           if (par2_allele(1)~=c2 | par2_allele(2)~=c2) 
               result = result*0.5
           end
           genotypeFactor.val(i) = result;
        elseif (ismember(c1, par2_allele) && ismember(c2, par1_allele))
            % c1 match par2, c2 match par1
            result = 1
           if (par2_allele(1)~=c1 | par2_allele(2)~=c1)
               result = 0.5
           end
           if (par1_allele(1)~=c2 | par1_allele(2)~=c2)
               result = result*0.5
           end
           genotypeFactor.val(i) = result;
        else
           genotypeFactor.val(i) = 0;
        end
    end           
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
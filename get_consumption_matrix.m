%% Constructs random consumption matrix
%
% Syntax:
%           [C] = get_consumption_matrix(n_row, n_col, type, varargin)
% Example:
%           [C] = get_consumption_matrix(100, 250, 'universe', 20, 'equal')
%           [C] = get_consumption_matrix(100, 250, 'nested', 20, 0.1, 0.9)
%
%   Where:
%      C               output matrix
%      numB and numD   dimensions of output matrix (numB x numD)
%      type            'universe'
%
%   Varargin:
%       nsubs           number of substrates per bacteria
%       dist_type       'equal', 'exponential' or 'non-normalized'  
%                       distribution of substrate preferences per bacteria
%
% Construct a random matrix, filled with postitive numbers between 0 and 1.
% Every row and column must have at least one entry > 0.
% Every row and column must have at least one entry == 0.

function [C] = get_consumption_matrix(n_row, n_col, type, varargin)

if strcmp(type, 'universe')
    % Fixed substrate preference distribution (uniform distribution or
    % exponentially declining preference).
    % In each row, the specified number of substrates with the chosen 
    % distribution is assigned to randomly chosen substrates.
    % The sum over rows is 1. The sum over columns differs.
    
    numB = n_row;
    numD = n_col;
    
    C = zeros(numB, numD);
    
    % number of substrates one bacterium can take up
    if nargin==5
        nsubs     = varargin{1};
        dist_type = varargin{2};
    else
        error('Number of input arguments wrong.')
    end
    assert(numD> nsubs,'C:numD_low',...
        'the number of DOM must > the number of substances to take up')
    
    % distribution of the substrate preferences
    switch dist_type
        case {'equal', 'uniform'}
            subs_dist = repmat(1/nsubs, 1,nsubs);
        case 'exponential'
            subs_dist = exp(-(1:nsubs));
            subs_dist = subs_dist/sum(subs_dist);
        case {'non-normalized'}
            subs_dist = ones(1, nsubs);
    end
    assert(numB*nsubs>=numD, 'C:nsubs_low', sprintf(...
        'not enough taken-up substrates to cover all DOMs (>=%d)',...
        ceil(numD/numB)))
    assert((numD-nsubs)*numB>=numD, 'C:nsubs_high',...
        'too many substrates taken up to have one zero per column')
    
    for i = 1:size(C,1)
        % randomly assign the distribution across substrates
        inds      = randperm(numD, nsubs);
        C(i,inds) = subs_dist;
    end
    
    % if a column contains only zeros or no zero at all, shuffle the order 
    % of entries in randomly chosen rows
    err_0    = sum(C>0)==0;
    err_numB = sum(C>0)==numB;
    while any(err_0) || any(err_numB)
        if any(err_0) 
            % choose a random row and move one of its positive (>0) entries 
            % to the respective column of this row
            indr      = randperm(numB, 1);
            indc_pos  = find(C(indr,:)>0);
            indc      = indc_pos(randperm(length(indc_pos),1));
            C(indr, find(err_0, 1, 'first')) = C(indr, indc);
            C(indr, indc) = 0;
        elseif any(err_numB)
            % choose a random row and move a randomly chosen zero to the
            % respective column
            indr = randperm(numB, 1);
            ind0 = find(C(indr,:)==0);
            % fill entry of the column into the zero spot
            C(indr, ind0(randperm(length(ind0),1))) =...
                C(indr, find(err_numB, 1, 'first'));
            % and fill zero into the column
            C(indr, find(err_numB, 1, 'first')) = 0;
        end
        err_0    = sum(C>0)==0;
        err_numB = sum(C>0)==numB;
    end
    %assert(all(sum(C,2)-1<10^-12), 'Rowsums should be 1.')
    
%     r = sum(C);
%     mean(r)
%     sum(subs_dist)*(numB/numD)
%     % The mean columnsum is sum(subs_dist)*(numB/numD)!
%     % Or for the relative case: numB/numD!
% 
%     r = sum(C);
%     mean(subs_dist)
%     [min(r(r>0)) max(r)]
%     numB

    % all rows and columns must have at least one entry >0
    assert(all(sum(C>0)>=1) & all(sum(C>0,2)>=1), ...
        'Row(s) or column(s) with only zeros detected.')

    % all rows and columns must have at least one entry == 0
    assert(all(sum(C==0)>0) & all(sum(C==0,2)>0), ...
        'Row(s) or column(s) with no zero detected.')

else
    error('undefined input type for consumption matrix')

end




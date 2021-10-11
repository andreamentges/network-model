%% Constructs excretion matrix
%
% Syntax:
%           [E] = get_excretion_matrix(C, excr_param)
%
%   Where
%      C    consumption matrix
%      E    excretion matrix, E(C>0) = 0
%      
% If excr_param < 1 : interpret as excretion density and generate matrix
% with specified percentage of nonzero entries.
% If excr_param>=1 : interpret as number of excreted compounds per
% bacteria, i.e. number of nonzero entries per row.
%
% Constructs excretion matrix E. Which must:
% - be zero where C>0
% - contain at least one entry >0 in every row and column
% - have rows summing up to 1

function [E] = get_excretion_matrix(C, excr_param)

% numB = 100;
% numD = 100;
% cons_width = 0.3
% excr_width = 0.5
% C = get_consumption_matrix(numB, numD, 'random', cons_width)

if excr_param<1
    
    % excretion density
    excr_width = excr_param;

    % calculate proportion of zeros in C, this proportion will be zero in E
    p = sum(sum(C==0))/(size(C,1)*size(C,2));
    if p<excr_width
        excr_width = p;
    %     warning(['consumption width is too high to implement specified ',...
    %     sprintf('excretion width. Reduced to %1.2f!', p)])
        error('E:nsubs_high', 'consumption width is too high')
    end

    retry = 1;
    count = 0;
    while retry == 1
        % realize correct percentage of zeros in C==0 places
        v_bin = rand(1,sum(sum(C==0)))<excr_width/p;

        % replaces ones by decimals
        v = v_bin .* rand(size(v_bin));

        % fill the random numbers in the "free slots" (where C==0)
        E       = NaN(size(C));
        E(C>0)  = 0;
        E(C==0) = v;

        retry = any(sum(E>0,1)==0) || any(sum(E>0,2)==0);
        count = count+1;
        assert(count<1000, 'E:retry',...
            'Generation of E not successfull after 1000 tries.',... 
            'Reduce consumption width.')
    end

    % rescale for columns to sum up to 1
    E = E./repmat(sum(E,2), 1, size(E,2)); % skewed distribution of entries
    
elseif excr_param>=1 & length(excr_param)==1
    
    numB = size(C,1);
    numD = size(C,2);
    
    % number of excreted compounds per bacteria
    nexcr = excr_param;
    
    % check whether specified parameters are feasible
    n_zeros_C = size(C,2)-sum(C>0,2);
    assert(all(n_zeros_C >= nexcr), ...
        'Not enough zeros in at least one row of C to implement specified number of non-zeros in E.')
    
    % potential entries of E are at zero entries in C
    E = zeros(size(C));
    potential = C==0;
     
    for i = 1:size(C,1)
        
        % number of potential entries in current row
        npot = sum(potential(i,:));
        
        % randomly assign entries across potential entries 
        inds      = randperm(npot);
        inds      = inds(1:nexcr);
        potinds   = find(potential(i,:)>0);
        E(i,potinds(inds)) = ones(1,nexcr);
    end
    
    % if a column contains only zeros or no zero at all, shuffle the order 
    % of entries in randomly chosen rows
    err_0    = sum(E>0)==0;
    err_numB = sum(E>0)==size(E,1);
    while any(err_0) || any(err_numB)
        if any(err_0) 
            % choose a random row and move one of its positive (>0) entries 
            % to the respective column of this row
            %(edit: chose only rows which have C==0 in the respective column)
            errcol    = find(err_0, 1, 'first');
            potential_rows = find(C(:,errcol)==0);
            temp = randperm(length(potential_rows));
            indr      = potential_rows(temp(1));
            indc_pos  = find(E(indr,:)>0);
            temp = randperm(length(indc_pos));
            indc = indc_pos(temp(1)); 
            E(indr, errcol) = E(indr, indc);
            E(indr, indc) = 0;
        elseif any(err_numB)
            % choose a random row and move a randomly chosen zero to the
            % respective column
            indr = randperm(numB);
            indr = indr(1);
            ind0 = find(E(indr,:)==0);
            % fill entry of the column into the zero spot
            temp = randperm(length(ind0),1);
            E(indr, ind0(temp(1))) =...
                E(indr, find(err_numB, 1, 'first'));
            % and fill zero into the column
            E(indr, find(err_numB, 1, 'first')) = 0;
        end
        err_0    = sum(E>0)==0;
        err_numB = sum(E>0)==size(E,1);
    end
    
    % make rows sum up to 1 (to conserve mass, taken up amount of DOM per 
    % bacteria must equal the amount of excreted DOM)
    E = E./repmat(sum(E,2), 1,  size(E,2));
       
    
end

%% security checks of E

% all rows and columns must have at least one entry >0
assert(all(sum(E>0)>=1) & all(sum(E>0,2)>=1), ...
    'Row(s) or column(s) with only zeros detected.')

% all rows and columns must have at least one entry ==0
assert(all(sum(E==0)>=1) & all(sum(E==0,2)>=1), ...
    'Row(s) or column(s) with no zeros detected.')

% all row sums must equal 1 (with accuracy 10^-8)
assert(all(abs(sum(E,2)-1)<10^-8), ...
    'Column(s) dont sum up to 1.')

% E cannot be >0 where C is >0 (bacteria cannot excrete what they consume)
assert(all(all(E(C>0)==0)), 'E is >0 where C is >0.')

    




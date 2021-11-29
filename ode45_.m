% Wrapper for ode45() returning non-interpolated values, only.
%
% Syntax
%
%    [t,y]     = ode45(odefun,tspan,y0)
%    [t,y]     = ode45(odefun,tspan,y0,options)
%    sol       = ode45(___)
%    [t,y,sol] = ode45(___)
%
% Description
%
%    Parameters and return values as described in ode45-documentation. 
%
%    Differences to ode45():
%
%       1) Return values t and y according to the steps chosen
%          by the solver, i.e. _no_ interpolation using deval().
%
%       2) Argument tspan with more than 2 elements does _not_ lead to
%          interpolation using deval(), but computes solution piecewise
%          using steps chosen by the solver, only.
%
% SEE ALSO: ode45

% Christoph Feenders, 2016-05-10 - 2016-06-07

function varargout = ode45_(odefun, tspan, y0, options, varargin)

    if ~exist('options', 'var')
        options = {};
    end

    sol = ode45(odefun, tspan([1 2]), y0, options, varargin{:});
%     sol = ode15s(odefun, tspan([1 2]), y0, options, varargin{:});
    
    if length(tspan) == 2
        
        t   = sol.x;
        y   = sol.y;
        
    else
    
        t   = sol.x([1 end]);
        y   = sol.y(:, [1 end]);
        for k = 3:length(tspan)
            if nargout == 2     % do not extend sol structure to save memory
                sol = ode45(odefun, tspan([k-1 k]), y(:,end), options, varargin{:});
            else                % extend sol structure
                sol = odextend(sol, [], tspan(k));
            end
            t(k)    = sol.x(end);
            y(:,k)  = sol.y(:, end);
        end
        
    end
    
    if nargout == 1
        varargout{1} = sol;
    elseif nargout >= 2
        varargout{1} = t';
        varargout{2} = y';
        if nargout == 3
            varargout{3} = sol;
        end
    end
    
end
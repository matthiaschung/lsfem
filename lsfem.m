function varargout = lsfem(fun, tspan, y0, param)
%
% function varargout = lsfem(fcn, tspan, y0, param)
%
% Authors:
%   (c) Matthias Chung (mcchung@vt.edu)
%       Justin Krueger (kruegej2@vt.edu)
%       Honghu Liu     (hhliu@vt.edu)
%
% Date: September 2021 (ver 1.0)
%
% MATLAB Version: 9.10.0.1710957 (R2021a) Update 4
%
% Description:
%   This function numerically solves the initial value problem
%              y' = fun(t,y) subject to y(t0) = y0
%
%   using least squares finite element methods by minimizing
%     min_x ||s'(t;x) - f(t,s(t;x))||^2 + ||s(t0;x) - y0||^2
%
% Input arguments:
%   fun             - right hand side of ODE
%   tspan           - time span [t0, tend] or vector of knots
%   y0              - initial conditions y(t0)
%   param           - optional parameter
%     splineDegree  - spline degree 1, 2, or 3 (default 3)
%     nknots        - number of equidistant finite elements (default 50)
%      npoints      - number of Gauss-Legendre points per element (default 4)
%      options      - optimset for Matlab's lsqnonlin function
%
% Output arguments:
%   t  - knots
%   y  - solution at knots structure
%   odestruc - structure of lsfem solution
%
% Example:
%   >> sol = lsfem(@(t,y) -y,[0 1], 1)
%
% References: M. Chung, J. Krueger, H. Liu, Least-squares finite element
% methods for ordinary differential equations" ArXiv:, 2021.
%

% set default parameters
splineDegree = 3; nknots = 50; npoints = 4; options = optimset('Display', 'off');

% rewrite default options if needed
if nargin == nargin(mfilename)
  for i = 1:size(param,1), eval([param{i,1},'= param{i,2};']); end
end

if nargout < 2
  out.solver = 'lsfem';
  out.extdata.odefun   = fun;
  out.extdata.options  = [];
  out.extdata.varargin = {};
end

nf = length(y0);                                      % spline dimension

if length(tspan) == 2                                 % generate nodes
  knots = linspace(tspan(1), tspan(end), nknots);     % uniform discretization
else                                                  % nodes already given
  knots = tspan; nknots = length(knots);
end

[t, w] = gaussLegendre(npoints, knots);               % generate Gauss-Legendre quadrature points and weights
s = femSpline(knots, [], t, splineDegree);            % construct spline structure

x = bsxfun(@times, y0, ones(nf, nknots));  x = x(:);  % initial guess for spline parameters replicate initial conditions
rFcn = @(x) residualFcn(x, fun, y0, nf, s, w);        % define residual function

x = lsqnonlin(rFcn, x, [], [], options);              % solve least squares problem

if nargout > 1                                        % outputs
  varargout{1} = knots;
  varargout{2} = reshape(x, nf, s.nknots);
else
  out.x = knots;
  out.y = reshape(x, nf, s.nknots);
  out.eval = @(t) femSpline(s, reshape(x, nf, s.nknots), t);
  out.spline = s;
  varargout{1} = out;
end

end

function r = residualFcn(x, fcn, y0, nf, s, w)

[S, Sd] = femSpline(s, reshape(x, nf, s.nknots), []);              % evaluate spline
r = bsxfun(@times, w, Sd - fcn(s.t, S)); r = [r(:); x(1:nf) - y0]; % compute residual

end

function [x, w] = gaussLegendre(n, knots)
%
% function [x, w] = gaussLegendre(n, interval)
%
% Authors:
%   (c) Matthias Chung (mcchung@vt.edu)
%       Justin Krueger (kruegej2@vt.edu)
%       Honghu Liu     (hhliu@vt.edu)
%
% Date: September 2021 (ver 1.0)
%
% MATLAB Version: 9.10.0.1710957 (R2021a) Update 4
%
% Description: 
%   This functions returns the evaluation nodes and corresponding weights
%   for a Gauss-Legendre quadrature rule with n <= 8. If necessary, the
%   function can return the translated nodes and weights for the intervals
%   passed in as a second input
%
% Input arguments:
%   n        - number of nodes for quadrature
%   interval - interval(s) on which to translate the nodes and weights
%
% Output arguments:
%   x       - quadrature nodes
%   w       - quadrature weights

% Example:
%   [x, w] = gaussLegendre =(3, [0 1 3 4])
%
% References:
%   [1] Numerical values from "Handbook of Mathematical Functions",
%   Abramowitz and Stegun, eds., 1965 Dover (reprint), Table 25.4, p. 916
%

% setup x and w
n = fix(n);         %  Make sure number of nodes is an integer
x = zeros(1, n);    %  Preallocate x vector
w = x;              %  Preallocate w vector

% table of node and weight values given n
switch n
  case 1
    x = 0;   
    w = 2;

  case 2
    x(1) = -0.577350269189626;     x(2) = -x(1);
    w(1) =  1;                     w(2) =  w(1);

  case 3
    x(1) = -0.774596669241483;     x(3) = -x(1);   
    x(2) =  0;      
    w(1) =  0.555555555555556;     w(3) =  w(1);
    w(2) =  0.888888888888889;    

  case 4
    x(1) = -0.861136311594053;     x(4) = -x(1);
    x(2) = -0.339981043584856;     x(3) = -x(2);
    w(1) =  0.347854845137454;     w(4) =  w(1);
    w(2) =  0.652145154862546;     w(3) =  w(2);

  case 5
    x(1) = -0.906179845938664;     x(5) = -x(1);
    x(2) = -0.538469310105683;     x(4) = -x(2);
    x(3) =  0;
    w(1) =  0.236926885056189;     w(5) =  w(1);
    w(2) =  0.478628670499366;     w(4) =  w(2);
    w(3) =  0.568888888888889;

  case 6
    x(1) = -0.932469514203152;     x(6) = -x(1);
    x(2) = -0.661209386466265;     x(5) = -x(2);
    x(3) = -0.238619186083197;     x(4) = -x(3);
    w(1) =  0.171324492379170;     w(6) =  w(1);
    w(2) =  0.360761573048139;     w(5) =  w(2);
    w(3) =  0.467913934572691;     w(4) =  w(3);

  case 7
    x(1) = -0.949107912342759;     x(7) = -x(1);
    x(2) = -0.741531185599394;     x(6) = -x(2);
    x(3) = -0.405845151377397;     x(5) = -x(3);
    x(4) =  0;
    w(1) =  0.129484966168870;     w(7) =  w(1);
    w(2) =  0.279705391489277;     w(6) =  w(2);
    w(3) =  0.381830050505119;     w(5) =  w(3);
    w(4) =  0.417959183673469;

  case 8
    x(1) = -0.960289856497536;     x(8) = -x(1);
    x(2) = -0.796666477413627;     x(7) = -x(2);
    x(3) = -0.525532409916329;     x(6) = -x(3);
    x(4) = -0.183434642495650;     x(5) = -x(4);
    w(1) =  0.101228536290376;     w(8) =  w(1);
    w(2) =  0.222381034453374;     w(7) =  w(2);
    w(3) =  0.313706645877887;     w(6) =  w(3);
    w(4) =  0.362683783378362;     w(5) =  w(4);

  otherwise
    error('Gauss quadrature with %d nodes not supported', n);

end

if nargin > 1  % translate nodes and weights if requested
    x = 0.5*(kron(diff(knots), x) + kron(knots(1:end-1)+knots(2:end), ones(1,n)));
    w = 0.5*kron(diff(knots), w); w = sqrt(w); 
end

end

function varargout = femSpline(varargin)
%
% function varargout = femSpline(varargin)
%
% Author:
%   (c) Matthias Chung (mcchung@vt.edu)
%       Justin Krueger (kruegej2@vt.edu)
%       Honghu Liu     (hhliu@vt.edu)
%
% Date: September 2021 (ver 1.0)
%
% MATLAB Version: 9.10.0.1710957 (R2021a) Update 4
%
% Description:
%   This function creates a structure or evaluates of linear, quadratic,
%   and, cubic spline that contains the spline's partial derivatives with
%   respect to its parameters and with respect to time and its parameters.
%   This information can then be used to evaluate the spline and the time
%   derivative of the spline as a function of the spline parameters.
%   Not-a-knot boundary conditions are used.
%
% Input arguments:
%   varargin -  input of spline structure  or evaluation points
%       knots   knot locations  (knot locations)
%       degree  choice of spline basis linear (1), quadratic (2), qubic (3)
%       t       evaluation points of spline
%       q       spline coefficients number of columns must match length of knots
%               can have multiple rows for approximation
%               dimesion of 1d function
%
% A spline structure s may have the following fields
% s =
%   struct with fields:
%            degree: degree of spline (1,2,3)
%             knots: knots of spline knots_1 < ... < knots_m
%            nknots: number of knots
%       coefficient: piecewise polynomial coefficients of s(t)
%     coefficientdt: piecewise polynomial coefficients of s'(t)
%                 t: evaluation points op spline s
%              dsdq: derivative of s with respect to spline coefficent at t
%            dsdqdt: derivative of s'(t) with respect to spline coefficent at t
%
% Output arguments:
%   varargout - returns either a spline structure or evaluates a given spline
%
% Example:
%
%  >> [sval, dsval, s] = femSpline(1:5, [sin(1:5); cos(1:5)], 1:0.5:5, 1)
%
% sval =
%     0.8415    0.8754    0.9093    0.5252    0.1411   -0.3078   -0.7568   -0.8579   -0.9589
%     0.5403    0.0621   -0.4161   -0.7031   -0.9900   -0.8218   -0.6536   -0.1850    0.2837
% dsval =
%     0.0678    0.0678   -0.7682   -0.7682   -0.8979   -0.8979   -0.2021   -0.2021   -0.2021
%    -0.9564   -0.9564   -0.5738   -0.5738    0.3363    0.3363    0.9373    0.9373    0.9373
% s =
%   struct with fields:
%
%            degree: 1
%             knots: [1 2 3 4 5]
%            nknots: 5
%       coefficient: [5×8 double]
%     coefficientdt: [5×4 double]
%                 t: [1 1.5000 2 2.5000 3 3.5000 4 4.5000 5]
%              dsdq: [5×9 double]
%            dsdqdt: [5×9 double]
%
% evaluates the spline s and its derivative dsdt with knots 1:5 and
% coefficients sin(1:5) for the first state and cos(1:5) for the second
% state at points 1:0.5:5, polynomail of degree 1 (default is 3). The third optional
% output parameter always returns the current spline structure.
%
% One is also able to just create a spline structure or overwrite parts of
% the spline structure without evaluating the spline. For instance
%
%  >> s = femSpline(1:10, 2)
%
% s =
%   struct with fields:
%            degree: 2
%             knots: [1 2 3 4 5 6 7 8 9 10]
%            nknots: 10
%       coefficient: [10×27 double]
%     coefficientdt: [10×18 double]
%
% generates a quadratic spline structure with knots 1,2,...,10 and
% piecewise polynomial coefficients (coefficient) of s and piecewise
% polynomial coefficients (coefficientdt) of dsdt. If the coefficent is
% kept empty [] the spline structure is generated:
%
% >> s = femSpline(1:5, [], 1:0.5:5, 3)
% s = 
%   struct with fields:
% 
%            degree: 3
%             knots: [1 2 3 4 5]
%            nknots: 5
%       coefficient: [5×16 double]
%     coefficientdt: [5×12 double]
%                 t: [1 1.5000 2 2.5000 3 3.5000 4 4.5000 5]
%              dsdq: [5×9 double]
%            dsdqdt: [5×9 double]
%
% the functionality also allows to update the spline structure. Assume
% the evaluation points t have not been determined or needed to be updated
% in a given structure s
%
% >> s = femSpline(s,[],1:0.1:5)
% s = 
%   struct with fields:
% 
%            degree: 3
%             knots: [1 2 3 4 5]
%            nknots: 5
%       coefficient: [5×16 double]
%     coefficientdt: [5×12 double]
%                 t: [1×41 double]
%              dsdq: [5×41 double]
%            dsdqdt: [5×41 double]
%
% Providing a non-empty second input argument always evaluates the spline. 
% Note if no evaluation points t exist in the spline structure an error
% occurs:
%
% >> s = femSpline(1:10, 2); femSpline(s,sin(1:10))
%
% Index exceeds the number of array elements (2).
% Error in femSpline (line 163)
%   elseif isempty(varargin{3}) 
%
% However, if a discretization exist the current spline can be
% evaluated with the coefficients as follows
%
% >> s = femSpline(1:5, [], 1:0.5:5); sval = femSpline(s,cos(1:5),[])
% sval =
%    0.5403    0.0603   -0.4161   -0.7969   -0.9900   -0.9293   -0.6536   -0.2277    0.2837
%
% In a similar fashion the discritization can be updated in the spline
% structure as follows
%
% >> s = femSpline(1:5, [], 1:0.5:5); s = femSpline(s,[],1:0.1:5)
% s = 
%   struct with fields:
% 
%            degree: 3
%             knots: [1 2 3 4 5]
%            nknots: 5
%       coefficient: [5×16 double]
%     coefficientdt: [5×12 double]
%                 t: [1×41 double]
%              dsdq: [5×41 double]
%            dsdqdt: [5×41 double]
%

if isstruct(varargin{1})
  if isempty(varargin{2})
    varargout{1} = constructPPMatrix(varargin{1}, varargin{3});
  elseif isempty(varargin{3})
    if nargout == 1
      varargout{1} = splineEval(varargin{1}, varargin{2});
    else
      [varargout{1}, varargout{2}] = splineEval(varargin{1}, varargin{2});
    end
  else
    s = constructPPMatrix(varargin{1}, varargin{3});
    if nargout == 1
      varargout{1} = splineEval(s, varargin{2});
    else
      [varargout{1}, varargout{2}] = splineEval(s, varargin{2});
    end
  end
else
  if nargin == 1
    varargout{1} = makeSplineCoefficients(varargin{1}, 3);
  elseif nargin == 2
    varargout{1} = makeSplineCoefficients(varargin{1}, varargin{2});
  elseif nargin == 3
    s = makeSplineCoefficients(varargin{1}, 3);
    s = constructPPMatrix(s, varargin{3});
    if isempty(varargin{2})
      varargout{1} = s;
    else
      if nargout == 1
        varargout{1} = splineEval(s, varargin{2});
      else
        [varargout{1}, varargout{2}] = splineEval(s, varargin{2});
      end
    end
  elseif nargin == 4
    s = makeSplineCoefficients(varargin{1}, varargin{4});
    s = constructPPMatrix(s, varargin{3});
    if isempty(varargin{2})
      varargout{1} = s;
    else
      if nargout == 1
        varargout{1} = splineEval(s, varargin{2});
      else
        [varargout{1}, varargout{2}] = splineEval(s, varargin{2});
      end
    end
  end
end
if nargout > 2
  varargout{3} = s;
end

end

function s = makeSplineCoefficients(knots, degree)

s.degree = degree; s.knots = knots; s.nknots = length(knots);
h = diff(s.knots)'; n = s.nknots-1; I = eye(n, n+1);

switch degree
  case 1 % linear spline
    % initialize constants for calculations
    Ja = -diag([1./h; 1].*ones(n+1,1)) + diag(1./h.*ones(n,1),1); Ja = Ja(1:end-1,:);
    s.coefficient   = [Ja; I]';
    s.coefficientdt = Ja';
    
  case 2 % quadratic spline
    J1 = [ -1/h(1), 1/h(1)+1/h(2), -1/h(2), zeros(1,n-2); [diag(h), zeros(n,1)] + [zeros(n,1),diag(h)]];
    J2 = diag([0, 2*ones(1,n)]) - diag(2*ones(1,n),-1);
    Jm = J1\J2;
    
    % initialize constants for calculations and calculate coefficients
    I0 = [zeros(n,1), eye(n)];
    Ja = bsxfun(@times, 1./(2*h), ones(1,n+1)).*(I0-I)*Jm;
    Jb = I*Jm;
    s.coefficient   = [Ja; Jb; I]';
    s.coefficientdt = [2*Ja; Jb]';
    
  case 3 % cubic spline
    lambda = h(2:n)./(h(1:n-1)+h(2:n));
    mu     = -6./(h(1:n-1).*h(2:n));
    
    % derivatives of the moments with respect to the knots    
    J1 = [-1/h(1), 1/h(1)+1/h(2), -1/h(2), zeros(1,n-2); ...
          [diag(1-lambda), zeros(n-1,2)] + [zeros(n-1,1), diag(2*ones(n-1,1)), zeros(n-1,1)] + [zeros(n-1,2), diag(lambda)]; ...
          zeros(1,n-2), -1/h(n-1), 1/h(n-1)+1/h(n), -1/h(n)];
    
    J2 = [zeros(1,n+1); ...
          [diag(-mu.*lambda), zeros(n-1,2)] + [zeros(n-1,1), diag(mu), zeros(n-1,1)] + [zeros(n-1,2), diag(-mu.*(1-lambda))]; ...
          zeros(1,n+1)];       
    Jm = J1\J2;
    
    % initialize constants for calculations and calculate coefficients
    H  = bsxfun(@times, h, ones(1,n+1));
    I0 = [zeros(n,1), eye(n)];
    Ja = 1./(6*H).*(I0-I)*Jm;
    Jb = 1/2*I*Jm;
    Jc = 1./H.*(I0-I) - H/6.*(2*I+I0)*Jm;
    s.coefficient   = [Ja; Jb; Jc; I]';
    s.coefficientdt = [3*Ja; 2*Jb; Jc]';
    
end

end

function s = constructPPMatrix(s,t) % calculate dsdq and dsdqdt
s.t = t;
s.dsdx   = ppval(mkpp(s.knots, s.coefficient,   s.nknots), t);
s.dsdxdt = ppval(mkpp(s.knots, s.coefficientdt, s.nknots), t);

end

function varargout = splineEval(s, x)
% q is required to be a dim x nknots matrix where
%  x = [x1@t1, x2@t1, ...., xn@t1, x1@t2, x2@t2, ...., xn@t2, ...]
% or for multiple states
% or  x = [x1@t1,  x1@t2, x1@t3, ...., x1@tm;
%          x2@t1,  x2@t2, x2@t3, ...., x2@tm
%           ...
%          xn@t1,  xn@t2, xn@t3, ...., xn@tm]

varargout{1} = x*s.dsdx;
if nargout > 1, varargout{2} = x*s.dsdxdt; end

end

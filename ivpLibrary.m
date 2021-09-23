function ivp = ivpLibrary(id,c)
%
% function P = ivpLibrary(id, c)
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
%   This function provides a library of initial value test problems of the
%   form:  y'(t) = f(t,y(t)) with y(t0) = y0.
%
% Input arguments:
%   choice - id number of initial value problem
%   c      - optional additional parameters for some initial value problems
%
% Output arguments:
%   ivp            - structure of initial value problem with the following fields
%     id           - id of IVP
%     name         - name of IVP problem
%     odefun       - handle of right-hand side (f(t,y))
%     dodefun      - Jacobian of right hand side of ODE (with respect to y)
%     y0           - initial condition of ivp
%     ivp.tspan    - time interval of interest
%     ivp.solution - analytic solution of ivp if exist 
%     ivp.linear   - boolean if ODE is linear or not
%     ivp.odestr   - 'string of right hand side of ODE
%     ivp.note     - additional information on the initial value problem
%     ivp.dim      - dimension of ode
%     ivp.linear   - structure of linear ODE Ey' = Ay + Bu(t) with
%           .E     - matrix E
%           .A     - matrix A
%           .B     - matrix B
%           .u     - function handle control u(t)
%
% Example:
%   >> ivpLibrary 
%
%   -------------------------------------------------
%   | id |      ode            |       name         |
%   -------------------------------------------------
%   |  1 | y' =              y |  exponential growth|
%   |  2 | y' =    y - 2e^{-t} |              Hairer|
%   |  3 | y' =   c1 y(1-y/c2) |   logistic equation|
%   |  4 | y' =      -y/(c1-y) |            kinetics|
%   -------------------------------------------------
%   ID of test problem? 
%
%   This provides a list of the initial value problems in this library and
%   request a selecton of a test problem
% 
%   >> ivpLibrary(4,0.01)
% ans = 
%   struct with fields:
% 
%           id: 4
%         name: 'kinetics'
%       odefun: @(t,y)-y./(c(1)+y)
%      dodefun: @(t,y)-c(1)./(c(1)+y).^2
%           y0: 1
%        tspan: [0 10]
%     solution: 'unknown'
%       linear: 0
%            c: 0.0100
%       odestr: '-y/(c1-y)'
%         note: 'kinetics problem'
%          dim: 1
%      odefuns: {@(t,y)-y./(c(1)+y)  @(t,y)-c(1)./(c(1)+y).^2}
%
% References:
%   M Chung, J Kruger, and H Liu. Least-Squares Finite Element Method for 
%   Ordinary Differential Equations,...
%

if nargin < 1
  fprintf('\n-------------------------------------------------\n')
  fprintf('| id |      ode            |       name         |\n')
  fprintf('-------------------------------------------------\n')
  for id = 1:4
    ivp = choose(id);
    fprintf('| %2d | y'' = %14s |%20s|\n',ivp.id,ivp.odestr,ivp.name)
  end
  fprintf('-------------------------------------------------\n') 
  id = input('ID of test problem? ');
  ivp = choose(id);
else
  if nargin < 2
    ivp = choose(id);
  else
    ivp = choose(id,c);
  end
end

end

function ivp = choose(id,c)

ivp.id = id;

switch id
  
  case 1
    ivp.name     = 'exponential growth';
    ivp.odefun   = @(t,y) y;
    ivp.dodefun  = @(t,y) eye(length(t));  % eye(nt, nt);
    ivp.y0       = 1;
    ivp.tspan    = [0, 2];
    ivp.solution = @(t) exp(t);
    ivp.odestr   = 'y';
    ivp.note     = 'standard exponential growth y'' = y wrt. y(0) = 1';
    ivp.linear.E =  1;
    ivp.linear.A =  1;
    ivp.linear.B =  0;
    ivp.linear.u = @(t) 0;
    
  case 2
    ivp.name     = 'Hairer';
    ivp.odefun   = @(t,y) y - 2*exp(-t);
    ivp.dodefun  = @(t,y) eye(length(t));
    ivp.y0       = 1;
    ivp.tspan    = [0, 30];
    ivp.solution = @(t) exp(-t);
    ivp.linear.E = 1;
    ivp.linear.A = 1;
    ivp.linear.B = 1;
    ivp.linear.u = @(t) - 2*exp(-t);
    ivp.odestr   = 'y - 2e^{-t}';
    ivp.note     = 'finite difference methods fail for this example';
    
  case 3
    ivp.name     = 'logistic equation';
    if nargin < 2, c = [1; 1]; end
    ivp.odefun   = @(t,y) c(1)*y.*(1-y/c(2));
    ivp.dodefun  = @(t,y) c(1)*(1-2*y/c(2));
    ivp.y0       = 0.1;
    ivp.tspan    = [0, 10];
    ivp.solution = @(t) c(2)*ivp.y0*exp(c(1)*t) ./ (c(2)+ivp.y0*(exp(c(1)*t)-1));
    ivp.linear   = 0;
    ivp.c        = c;
    ivp.odestr   = 'c1 y(1-y/c2)';
    ivp.note     = 'standard logistic equation';
    
  case 4
    if nargin < 2, c = 0.1; end
    ivp.name     = 'kinetics';
    ivp.odefun   = @(t,y) -y./(c(1)+y);
    ivp.dodefun  = @(t,y) -c(1)./(c(1)+y).^2;
    ivp.y0       = 1;
    ivp.tspan    = [0, 10];
    ivp.solution = 'unknown';
    ivp.linear   = 0;
    ivp.c        = c;
    ivp.odestr   = '-y/(c1-y)';
    ivp.note     = 'kinetics problem';
  
  otherwise
    error('Not a valid id.')
    
end

if isstruct(ivp.linear)
  ivp.odefun   = @(t,y) ivp.linear.A*y + ivp.linear.B*ivp.linear.u(t);
  ivp.linear.def = ' Ey''(t) = Ay(t) + Bu(t) ';
  ivp.dodefun   = @(t,y) ivp.linear.A*y + ivp.linear.B*ivp.linear.u(t);
end

ivp.dim = length(ivp.y0);
ivp.odefuns{1} = ivp.odefun; ivp.odefuns{2} = ivp.dodefun;

end

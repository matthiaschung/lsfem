% driver_lsfem.m
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
% References: M. Chung, J. Krueger, H. Liu, Least-squares finite element
% methods for ordinary differential equations" ArXiv:, 2021.

clear, close all

% standard problems
id = 4;
ivp = ivpLibrary(id);
%   | id |      ode            |       name         |
%   |  1 | y' =              y |  exponential growth|
%   |  2 | y' =    y - 2e^{-t} |              Hairer|
%   |  3 | y' =   c1 y(1-y/c2) |   logistic equation|
%   |  4 | y' =      -y/(c1-y) |            kinetics|

% solve ode
fdmSol = ode45(ivp.odefun,ivp.tspan,ivp.y0); % use ode45
femSol = lsfem(ivp.odefun,ivp.tspan,ivp.y0); % use lsfem

% evaluate ODE solution
t = linspace(ivp.tspan(1),ivp.tspan(end),500);
yfdm = deval(fdmSol,t);
yfem = femSol.eval(t);

% plot solutions
subplot(2,1,1)
plot(t,yfdm), hold on
plot(t,yfem)
leg = {fdmSol.solver,femSol.solver};
if ~strcmp(ivp.solution,'unknown')
    ytrue = ivp.solution(t); % compute analytic solution
    leg{3} = 'analytic';
    plot(t,ytrue)
end
legend(leg)
ylabel('$y^h$')
title(ivp.name)

% plot error
subplot(2,1,2)
if ~strcmp(ivp.solution,'unknown')
    semilogy(t,abs(yfdm - ytrue)),hold on
    semilogy(t,abs(yfem - ytrue))
    ylabel('abs. error [log]')
else
    semilogy(t,abs(yfem-yfdm))
    ylabel('$|y_{\rm fem} - y_{\rm fdm}|$ [log]', 'Interpreter','Latex')
    legend('$y_{\rm fem} - y_{\rm fdm}$', 'Interpreter','Latex')
end
xlabel('$t$')


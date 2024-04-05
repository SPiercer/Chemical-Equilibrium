clc; clear; close all
% Faculty of Engineering
% Cairo University
% 3rd Year Mechanical Power Engineering
% MEP-3080 Fundamentals of Combustion
% *************************************************************************
% Required to find the Adiabatic Flame Temperature(ADF)and the chemical
% combostion for the followig reaction:
% *************************************************************************
% (phi)CxHy + (x+y/4)(O2+3.76N2) >>> (a)Co2 + (b)Co + (d)H2O + (e)H2 +(f)N2
% *************************************************************************
% Notes:
% *******
% The Adiabatic Flame Temperature is at Q=0.
% Q = The Amount of Heat.
% phi = Equivalence ratio.
% The reaction is rich (No excess Oxygen).
% Water-gas shift equation: CO2 + H2 >>> CO + H2O.

% declaring global variables
global x y phi z N O f;

% declaring user inputs

% The number of Carbon atoms in the fuel
x = safeInput('Enter x (The number of carbon atoms in the fuel): ');
% The number of Hydrogen atoms in the fuel
y = safeInput('Enter y (The number of hydrogen atoms in the fuel): ');

% The Equivalence ratio
phi = safeInput('Enter phi (The equivalence ratio Φ): ');;

% The number of moles of air
z = x + (y / 4);
delayedDisp(['The number of moles of air is: ', num2str(z)]);

% The number of Nitrogen moles in reactants (air)
N = 3.76*z;
delayedDisp(['The number of Nitrogen moles in reactants (air) is: ', num2str(N)]);

% The number of Nitrogen moles in reactants (air)
O = 1*z;
delayedDisp(['The number of Oxygen moles in reactants (air) is: ', num2str(O)]);

% N Balance
f = 1/2*(O * 3.76 * 2);
delayedDisp(['The number of moles of N2 in products is: ', num2str(f)]);
longBreak();

% C Balance
delayedDisp('We balance the carbon moles by using this equation');
delayedDisp('Φ * x = a + b');
delayedDisp('hence: b = Φ*x - a, (Number of moles of CO)');
longBreak();

% O Balance
delayedDisp('We balance the oxygen moles by using this equation');
delayedDisp('2z = 2a + b + d');
delayedDisp('hence: d = 2z - 2a - b, (Number of moles of H2O)');
longBreak();

% H Balance
delayedDisp('We balance the hydrogen moles by using this equation');
delayedDisp('Φ * y = 2d + 2e');
delayedDisp('hence: e = 1/2 (Φ * y - 2d), (Number of moles of H2)');
longBreak();

% Water-gas shift equation
delayedDisp("Now we'll calculate the equilbrium constant (Kp) using this equation:");
delayedDisp('Kp = (X.CO^v.co * X.H2O^v.h2o / X.CO2^v.co2 * X.H2^v.h2) * (P.total / P.atm)^v.co+v.h2o-v.co2-v.h2');
pause(1);
delayedDisp('where X is the number of moles of the element')
delayedDisp('and P.total and P.atm are the pressure of products and atmospheric (reactants) respectivly')
longBreak();

delayedDisp('Subsituting the number of moles accordingly and knowing that v = 1 for all elements');
longBreak();

delayedDisp('We get:')
delayedDisp('Kp = b * d / a * e')
longBreak();

% Calculating the chemical composition and amount of heat by assuming a temperture of the products.
delayedDisp('since it is a rich flame, then the temperature (T) is larger than 1800 k')
longBreak();

T1 = safeInput('Please enter the first assumed temperature of the products (T1): ');
[a1,b1,d1,e1] = chemicalComposition(T1);

T2 = safeInput('Please enter the second assumed temperature of the products (T2): ');
[a2,b2,d2,e2] = chemicalComposition(T2);

function [a,b,d,e] = chemicalComposition(T)
% getting global variables
global x y phi z f;

lnKp = lnKpConst(T);

Kp = exp(lnKp);
delayedDisp(['The equilibrium constant (Kp) is: ', num2str(Kp)]);
longBreak();

% now we can calculate the chemical composition
delayedDisp('Subsituting the values of Kp, a, b, d, e in the equation Kp = b * d / a * e');
delayedDisp('Since all the values are a function of a, we can calculate a');
delayedDisp('Kp = (Φ*x - a) * (2z - 2a - (Φ*y - 2d)) / a * (1/2 (Φ*y - 2(2z - 2a - (Φ*x - a))))');
delayedDisp('Solving the equation will give us the value of a');
longBreak();

% solve for a
syms a;
eqn = Kp == ((phi * x - a) * (2 * z - 2 * a - (phi * x - a)) / (a * (1 / 2) * (phi * y - 2 * (2 * z - 2 * a - (phi * x - a)))));
S = solve(eqn, a,'Real',true);

a = double(S(S>0));
delayedDisp(['The value of a is: ', num2str(a)]);

% calculating the rest of the values
b = phi*x - a;
d = 2*z - 2*a - b;
e = 1/2 * (phi*y - 2*d);
delayedDisp(['The value of b is: ', num2str(b)]);
delayedDisp(['The value of d is: ', num2str(d)]);
delayedDisp(['The value of e is: ', num2str(e)]);

delayedDisp(['recalling the value of f is: ', num2str(f)]);
longBreak();
end

% helper functions

% A mini equilibrium constants table
function [lnKp] = lnKpConst(T)
switch T
    case 2400
        lnKp = 1.847;
    case 2000
        lnKp = 1.51;
    otherwise
        lnKp = safeInput(['Please enter the lnKp according to (', num2str(T), 'K) in your table of natural log of equilibrium constants: ']);
end
end

% safe input function to check for invalid inputs
function [x] = safeInput(prompt)
x = input(prompt);
while isempty(x)
    disp('Invalid input, please enter a valid number');
    x = input(prompt);
end
end

%  delayed display function to display the input after a delay
function delayedDisp(input)
% pause(0.5);
disp(input);

end

function longBreak()
disp(" ")
% pause(1)
end
clear

SPLINE = true;
Nz = 10000;
n = 1;
zex = 5.264; % cm
tend = 100;
Ib = 35; % A
R0 = 2.077; % cm
Rb = 0.827; % cm
g = 1.2;
ukv = 100; % [kV]
gamma = 1.0 + ukv/511.0;
betta = sqrt(1.0d0 - 1.0/(gamma*gamma));
betta2 = 1.0d0 - 1.0/(gamma*gamma);
betta_z = betta/sqrt(g*g + 1.0d0);
betta_z2 = betta2/(g*g + 1.0d0);
betta_perp2 = betta2 - betta_z2;
gamma0 = sqrt(1 - betta_perp2 - betta_z);

c = 29979245800; % [cm/s]
e = 4.803e-10; % [ед. СГС]
m = 9.1093837015e-28; % [g]

% c = 299792458; % [m/s]
% e = 1.6021766208e-19; % [Кл]
% m = 9.1093837015e-31; % [кг]

% Вычислене nu (мода ТЕ28.12)
syms x;
dydx = matlabFunction(diff(besselj(28, x), x));
% z=0:0.01:100;
% plot(z,dydx(z))
approx = 74;
nu = fzero(dydx, approx);

w_op = c * nu / R0;
ZetaEx = betta_perp2/2.0/betta_z*w_op*zex/c;


I0 = 16 / (17045.81697831) * (Ib * betta_z * besselj(28-1,nu*Rb/R0)^2) / ...
    (gamma0 * betta^6 * (nu^2 - 28^2) * besselj(28, nu)^2);


TauEnd = betta_perp2^2*w_op*tend/8/betta_z2;

% delta(z)
dz = ZetaEx/(Nz-1);
Rr_file = load('d:\Alex\Documents\Work\novozhilova\22-09-23 Продольная структура и профиль рез-ра\RR2812.dat');
ZAxis = 0:dz:ZetaEx;
Rr = zeros(Nz,1);
[b, c, d] = spline(length(Rr_file(:,1)), Rr_file(:,1), Rr_file(:,2));

x = ZAxis/ZetaEx*(Rr_file(end,1) - Rr_file(1,1)) + Rr_file(1,1);
if (SPLINE == true)
    for i=1:length(x)
        %    y(i) = uval(x(i),Rr(:,1),Rr(:,2));
        Rr(i) = seval(x(i), length(Rr_file(:,1)), Rr_file(:,1), Rr_file(:,2), b, c, d);
    end
else
    for i=1:length(x)
           Rr(i) = uval(x(i),Rr_file(:,1),Rr_file(:,2));        
    end
end

plot(ZAxis,Rr)
dz
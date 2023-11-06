function [OUTFre, OUTFim, OUTJre, OUTJim, OUTZAxis, OUTTAxis] = ...
    gyrotron(Ne, Nz, Lz, Tend, Delta, I0, R0, Rb, g, ukv, dz, dt, tol, INTT, INTZ, SPLINE) %#codegen

if nargin < 8
    fprintf('USAGE: orotron Ne Lz Tend Delta I dz dt\n')
end

convert = true;
if Nz == 0
    Nz = fix(Lz/dz) + 1;
    Nt = fix(Tend/dt) + 1;
    convert = false;
end

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

if convert == true
    ZetaEx = betta_perp2/2.0/betta_z*w_op*Lz/c;
    I0 = 16 / (17000) * (I0 * betta_z * besselj(28-1,nu*Rb/R0)^2) / ...
        (gamma0 * betta^6 * (nu^2 - 28^2) * besselj(28, nu)^2);
    
    TauEnd = betta_perp2^2*w_op*Tend/8/betta_z2;
    
    % delta(z)
    dz = ZetaEx/(Nz-1);
    dt = dz*dz;
    Nt = fix(TauEnd/dt/100000000) + 1;
    Rr_file = load('d:\Alex\Documents\Work\novozhilova\22-09-23 Продольная структура и профиль рез-ра\RR2812.dat');
    ZAxis = 0:dz:ZetaEx;
    Rr = zeros(Nz,1);
    [ba, ca, da] = spline(length(Rr_file(:,1)), Rr_file(:,1), Rr_file(:,2));
    
    x = ZAxis/ZetaEx*(Rr_file(end,1) - Rr_file(1,1)) + Rr_file(1,1);
    if (SPLINE == true)
        for i=1:length(x)
            %    y(i) = uval(x(i),Rr(:,1),Rr(:,2));
            Rr(i) = seval(x(i), length(Rr_file(:,1)), Rr_file(:,1), Rr_file(:,2), ba, ca, da);
        end
    else
        for i=1:length(x)
            Rr(i) = uval(x(i),Rr_file(:,1),Rr_file(:,2));
        end
        Rr(:) = Rr(:)/10; % в миллиметрах
    end  
    
    wc = c*nu*ones(Nz,1)./Rr;    
    kpar2 = 8*betta_z2/betta_perp2^2*(1 - wc/w_op);        
else
    ZetaEx = Lz;
    TauEnd = Tend;
    kpar2 = zeros(Nz,1);
end

ZAxis = zeros(Nz, 1);
TAxis = zeros(Nt, 1);
InitialField = zeros(Nz,1);

for i=1:Nz
    ZAxis(i) = (i-1) * dz;
end

for i=1:Nt
    TAxis(i) = (i-1) * dt;
end

ZBEG = 0;
% ZEND = .5;
ZEND = ZetaEx;
IND1 = (ZAxis > ZBEG & ZAxis < ZEND);
% InitialField(IND1,1) = 0.001*sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;
InitialField(IND1,1) = 0.01*sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;
% InitialField = ones(length(ZAxis),1) + 1i*ones(length(ZAxis),1);

% infield=[real(InitialField) imag(InitialField)];
%save('init_field.in','infield','-ascii')

if INTT < 1
    error('Too small Tend');
end

if INTZ < 1
    error('Too small Zend');
end

if INTT > 1 && INTZ > 1
    OUTNt = fix(Nt/INTT) + 1;
    OUTNz = fix(Nz/INTZ) + 1;
elseif INTT == 1 && INTZ > 1
    OUTNt = Nt;
    OUTNz = fix(Nz/INTZ) + 1;
elseif INTT > 1 && INTZ == 1
    OUTNt = fix(Nt/INTT) + 1;
    OUTNz = Nz;
else
    OUTNt = Nt;
    OUTNz = Nz;
end

%OUTJ = zeros(OUTNz,OUTNt);
%OUTF = zeros(OUTNz,OUTNt);
OUTZAxis = zeros(OUTNz,1);
OUTTAxis = zeros(OUTNt,1);

OUTZAxis(1) = 0;
for i = 2:OUTNz
    OUTZAxis(i) = (i-1)*INTZ*dz;
end

OUTTAxis(1) = 0;
for i = 2:OUTNt
    OUTTAxis(i) = (i-1)*INTT*dt;
end

fileID = fopen('input.txt','w');
fprintf(fileID,'Nz = %f\n', Nz);
fprintf(fileID,'Nt = %f\n', Nt);
fprintf(fileID,'Ne = %f\n', Ne);
fprintf(fileID,'ZetaEx = %f\n', ZetaEx);
fprintf(fileID,'TauEnd = %f\n', TauEnd);
fprintf(fileID,'Delta = %f\n', Delta);
fprintf(fileID,'I0 = %f\n', I0);
fprintf(fileID,'R0 = %f\n', R0);
fprintf(fileID,'Rb = %f\n', Rb);
fprintf(fileID,'g = %f\n', g);
fprintf(fileID,'ukv = %f\n', ukv);
fprintf(fileID,'dz = %f\n', dz);
fprintf(fileID,'dt = %f\n', dt);
fprintf(fileID,'tol = %g\n', tol);
fclose(fileID);

[OUTF, OUTJ] = gyroscr(Nz, Nt, Ne, ZAxis, TAxis, Delta, I0, dt, dz, tol, kpar2, INTT, INTZ, OUTNz, OUTNt, InitialField);

Folder = 'results/';

hash = datestr(now,30);
FolderName=sprintf('%s%s', Folder, hash);
mkdir(FolderName);

% Вывод в .mat файл
%     fileResults = sprintf('%s/%s', FolderName, 'results.mat');
%     save(fileResults,"RES","-v7.3");

% Вывод в .dat файл
OUTFre = real(OUTF);
OUTFim = imag(OUTF);
OUTJre = real(OUTJ);
OUTJim = imag(OUTJ);

OUTBvsZ = zeros(size(OUTF,1), 2*size(OUTF,2));
OUTJvsZ = zeros(size(OUTJ,1), 2*size(OUTJ,2));
OUTBvsZ(:,1:2:end-1) = real(OUTF);
OUTBvsZ(:,2:2:end) = imag(OUTF);
OUTJvsZ(:,1:2:end-1) = real(OUTJ);
OUTJvsZ(:,2:2:end) = imag(OUTJ);

OUTF = OUTF.';
OUTJ = OUTJ.';
OUTBvsT = zeros(size(OUTF,1), 2*size(OUTF,2));
OUTJvsT = zeros(size(OUTJ,1), 2*size(OUTJ,2));
OUTBvsT(:,1:2:end-1) = real(OUTF);
OUTBvsT(:,2:2:end) = imag(OUTF);
OUTJvsT(:,1:2:end-1) = real(OUTJ);
OUTJvsT(:,2:2:end) = imag(OUTJ);

OUTBvsZ = [OUTZAxis OUTBvsZ];
OUTJvsZ = [OUTZAxis OUTJvsZ];

OUTBvsT = [OUTTAxis OUTBvsT];
OUTJvsT = [OUTTAxis OUTJvsT];

fileParameters = sprintf('%s/%s', FolderName, 'parameters.txt');
fileResultsBvsZ = sprintf('%s/%s', FolderName, 'resultsBvsZ.dat');
fileResultsJvsZ = sprintf('%s/%s', FolderName, 'resultsJvsZ.dat');
fileResultsBvsT = sprintf('%s/%s', FolderName, 'resultsBvsT.dat');
fileResultsJvsT = sprintf('%s/%s', FolderName, 'resultsJvsT.dat');
fileT = sprintf('%s/%s', FolderName, 'Time.dat');
fileZ = sprintf('%s/%s', FolderName, 'Z.dat');

strParam = 'Ne = %d\nTend = %f\nLz = %f\nDelta = %f\nI = %f\ndz = %f\ndt = %f\n';
fileID = fopen(fileParameters ,'wt');
fprintf(fileID, strParam, Ne, Tend, Lz, Delta, Ic, dz, dt);
fclose(fileID);
save(fileResultsBvsZ, 'OUTBvsZ', '-ascii');
save(fileResultsJvsZ, 'OUTJvsZ', '-ascii');
save(fileResultsBvsT, 'OUTBvsT', '-ascii');
save(fileResultsJvsT, 'OUTJvsT', '-ascii');
save(fileT, 'TAxis', '-ascii');
save(fileZ, 'ZAxis', '-ascii');
end

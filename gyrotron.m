function [OUTFre, OUTFim, OUTJre, OUTJim, OUTZAxis, OUTTAxis] = gyrotron(NNe, LLz, TTend, DDelta, Icc, ddz, ddt, tol) %#codegen

if nargin < 8
    fprintf('USAGE: orotron Ne Lz Tend Delta I dz dt\n')
end

fprintf('\nNe = %f\n', NNe);
fprintf('Lz = %f\n', LLz);
fprintf('Tend = %f\n', TTend);
fprintf('Delta = %f\n', DDelta);
fprintf('Ic = %f\n', Icc);
fprintf('dz = %f\n', ddz);
fprintf('dt = %f\n', ddt);
fprintf('tol = %g\n', tol);

% Ne = str2double(NNe);
% Lz = str2double(LLz);
% Tend = str2double(TTend);
% Delta = str2double(DDelta);
% Ic = str2double(Icc);
% dz = str2double(ddz);
% dt = str2double(ddt);

Ne = NNe;
Lz = LLz;
Tend = TTend;
Delta = DDelta;
Ic = Icc;
dz = ddz;
dt = ddt;

Nz = fix(Lz/dz) + 1;
Nt = fix(Tend/dt) + 1;

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
ZEND = 0.5;
IND1 = (ZAxis > ZBEG & ZAxis < ZEND);
% InitialField(IND1,1) = 0.001*sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;
InitialField(IND1,1) = sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;
% InitialField = ones(length(ZAxis),1) + 1i*ones(length(ZAxis),1);

infield=[real(InitialField) imag(InitialField)];
%save('init_field.in','infield','-ascii')

INTT = fix(Nt/500);
INTZ = fix(Nz/500);
% INTT = 1;
% INTZ = 1;
if INTT < 1
    error('Too small Tend');    
end

if INTZ < 1
    error('Too small Zend');    
end

if INTT > 1 && INTZ > 1
    OUTNt = fix((Nt-1)/INTT) + 1;
    OUTNz = fix((Nz-1)/INTZ) + 1;
elseif INTT == 1 && INTZ > 1
    OUTNt = Nt;
    OUTNz = fix((Nz-1)/INTZ) + 1;
elseif INTT > 1 && INTZ == 1
    OUTNt = fix((Nt-1)/INTT) + 1;
    OUTNz = Nz;
else
    OUTNt = Nt;
    OUTNz = Nz;
end    

%OUTJ = zeros(OUTNz,OUTNt);
%OUTF = zeros(OUTNz,OUTNt);
OUTZAxis = zeros(OUTNz,1);
OUTTAxis = zeros(OUTNt,1);

for i = 1:OUTNz
    OUTZAxis(i) = (i-1)*INTZ*dz;
end

for i = 1:OUTNt
    OUTTAxis(i) = (i-1)*INTT*dt;
end

[OUTF, OUTJ] = gyroscr(Nz, Nt, Ne, ZAxis, TAxis, Delta, Ic, dt, dz, tol, INTT, INTZ, OUTNz, OUTNt, InitialField);

% plot(RES.ZAxis(:,1), abs(RES.OUTF(:,end)))

% if (isempty(TAxis) || isempty(OUTF))
%     msgbox('First press the button "Solve"!')
% else
%     Folder = 'results/';
%     
%     hash = datestr(now,30);
%     FolderName=sprintf('%s%s', Folder, hash);
%     mkdir(FolderName);
    
    % Âûגמה ג .mat פאיכ
%     fileResults = sprintf('%s/%s', FolderName, 'results.mat');
%     save(fileResults,"RES","-v7.3");
    
    % Âûגמה ג .dat פאיכ    
    OUTFre = real(OUTF);
    OUTFim = imag(OUTF);
    OUTJre = real(OUTJ);
    OUTJim = imag(OUTJ);
    
%     OUTBvsZ = zeros(size(OUTF,1), 2*size(OUTF,2));
%     OUTJvsZ = zeros(size(OUTJ,1), 2*size(OUTJ,2));
%     OUTBvsZ(:,1:2:end-1) = real(OUTF);
%     OUTBvsZ(:,2:2:end) = imag(OUTF);
%     OUTJvsZ(:,1:2:end-1) = real(OUTJ);
%     OUTJvsZ(:,2:2:end) = imag(OUTJ);
%     
%     OUTF = OUTF.';
%     OUTJ = OUTJ.';
%     OUTBvsT = zeros(size(OUTF,1), 2*size(OUTF,2));
%     OUTJvsT = zeros(size(OUTJ,1), 2*size(OUTJ,2));
%     OUTBvsT(:,1:2:end-1) = real(OUTF);
%     OUTBvsT(:,2:2:end) = imag(OUTF);
%     OUTJvsT(:,1:2:end-1) = real(OUTJ);
%     OUTJvsT(:,2:2:end) = imag(OUTJ);      
%     
%     OUTBvsZ = [OUTZAxis OUTBvsZ];
%     OUTJvsZ = [OUTZAxis OUTJvsZ];
%     
%     OUTBvsT = [OUTTAxis OUTBvsT];
%     OUTJvsT = [OUTTAxis OUTJvsT];
%     
%     fileParameters = sprintf('%s/%s', FolderName, 'parameters.txt');    
%     fileResultsBvsZ = sprintf('%s/%s', FolderName, 'resultsBvsZ.dat');
%     fileResultsJvsZ = sprintf('%s/%s', FolderName, 'resultsJvsZ.dat');
%     fileResultsBvsT = sprintf('%s/%s', FolderName, 'resultsBvsT.dat');
%     fileResultsJvsT = sprintf('%s/%s', FolderName, 'resultsJvsT.dat');
%     fileT = sprintf('%s/%s', FolderName, 'Time.dat');
%     fileZ = sprintf('%s/%s', FolderName, 'Z.dat');
%     
%     strParam = 'Ne = %d\nTend = %f\nLz = %f\nDelta = %f\nI = %f\ndz = %f\ndt = %f\n';
%     fileID = fopen(fileParameters ,'wt');
%     fprintf(fileID, strParam, Ne, Tend, Lz, Delta, Ic, dz, dt);
%     fclose(fileID);    
%     save(fileResultsBvsZ, 'OUTBvsZ', '-ascii');
%     save(fileResultsJvsZ, 'OUTJvsZ', '-ascii');
%     save(fileResultsBvsT, 'OUTBvsT', '-ascii');
%     save(fileResultsJvsT, 'OUTJvsT', '-ascii');
%     save(fileT, 'TAxis', '-ascii');
%     save(fileZ, 'ZAxis', '-ascii');
% end

return

end

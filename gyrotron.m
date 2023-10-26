function gyrotron(NNe, LLz, TTend, DDelta, II, ddz, ddt, tol) %#codegen

if nargin < 7
    fprintf('USAGE: orotron Ne Lz Tend Delta I dz dt\n')
end

Ne = NNe;
Lz = LLz;
Tend = TTend;
Delta = DDelta;
I = II;
dz = ddz;
dt = ddt;
DeltaZ = dz;
DeltaT = dt;

ZAxis = (0:DeltaZ:Lz)';
TAxis = (0:DeltaT:Tend)';

InitialField = zeros(length(ZAxis),1);
ZBEG = 0;
ZEND = 0.5;
IND1 = (ZAxis > ZBEG & ZAxis < ZEND);
InitialField(IND1,1) = 0.001*sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;
% InitialField(IND1,1) = sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;

infield=[real(InitialField) imag(InitialField)];
save('init_field.in','infield','-ascii')

INTERVALT = fix(size(TAxis,1)/500);
INTERVALZ = fix(size(ZAxis,1)/500);
% INTERVALT = 1;
% INTERVALZ = 1;
if INTERVALT < 1
    error('Too small Tend');    
end

if INTERVALZ < 1
    error('Too small Zend');    
end

IN.Ne = Ne;
IN.ZAxis = ZAxis;
IN.TAxis = TAxis;
IN.Delta = Delta;
IN.I = I;
IN.INTERVALT = INTERVALT;
IN.INTERVALZ = INTERVALZ;
IN.InitialField = InitialField;
IN.tol = tol

OUT = gyroscr(IN);

RES.Ne = IN.Ne;
RES.Tend = Tend;
RES.Lz = Lz;
RES.Delta = IN.Delta;
RES.I = IN.I;
RES.dz = dz;
RES.dt = dt;
RES.INTERVALT = IN.INTERVALT;
RES.INTERVALZ = IN.INTERVALZ;
RES.ZAxis = OUT.ZAxis;
RES.TAxis = OUT.TAxis;
RES.OUTB = OUT.OUTB;
RES.OUTJ = OUT.OUTJ;
% OUTFieldRe = real(OUT.OUTB);
% OUTFieldIm = imag(OUT.OUTB);
% OUTJRe = real(OUT.OUTJ);
% OUTJIm = imag(OUT.OUTJ);

% plot(RES.ZAxis(:,1), abs(RES.OUTB(:,end)))

if (isempty(RES.TAxis) || isempty(RES.OUTB))
    msgbox('First press the button "Solve"!')
else
    Folder = 'results/';
    
    hash = datestr(now,30);
    FolderName=sprintf('%s%s', Folder, hash);
    mkdir(FolderName);
    
    % Âûגמה ג .mat פאיכ
%     fileResults = sprintf('%s/%s', FolderName, 'results.mat');
%     save(fileResults,"RES","-v7.3");
    
    % Âûגמה ג .dat פאיכ    
    OUTBvsZ = zeros(size(OUT.OUTB,1), 2*size(OUT.OUTB,2));
    OUTJvsZ = zeros(size(OUT.OUTJ,1), 2*size(OUT.OUTJ,2));
    OUTBvsZ(:,1:2:end-1) = real(OUT.OUTB);
    OUTBvsZ(:,2:2:end) = imag(OUT.OUTB);
    OUTJvsZ(:,1:2:end-1) = real(OUT.OUTJ);
    OUTJvsZ(:,2:2:end) = imag(OUT.OUTJ);
    
    OUT.OUTB = OUT.OUTB.';
    OUT.OUTJ = OUT.OUTJ.';
    OUTBvsT = zeros(size(OUT.OUTB,1), 2*size(OUT.OUTB,2));
    OUTJvsT = zeros(size(OUT.OUTJ,1), 2*size(OUT.OUTJ,2));
    OUTBvsT(:,1:2:end-1) = real(OUT.OUTB);
    OUTBvsT(:,2:2:end) = imag(OUT.OUTB);
    OUTJvsT(:,1:2:end-1) = real(OUT.OUTJ);
    OUTJvsT(:,2:2:end) = imag(OUT.OUTJ);      
    
    OUTBvsZ = [OUT.ZAxis OUTBvsZ];
    OUTJvsZ = [OUT.ZAxis OUTJvsZ];
    
    OUTBvsT = [OUT.TAxis OUTBvsT];
    OUTJvsT = [OUT.TAxis OUTJvsT];
    
    fileParameters = sprintf('%s/%s', FolderName, 'parameters.txt');    
    fileResultsBvsZ = sprintf('%s/%s', FolderName, 'resultsBvsZ.dat');
    fileResultsJvsZ = sprintf('%s/%s', FolderName, 'resultsJvsZ.dat');
    fileResultsBvsT = sprintf('%s/%s', FolderName, 'resultsBvsT.dat');
    fileResultsJvsT = sprintf('%s/%s', FolderName, 'resultsJvsT.dat');
    fileT = sprintf('%s/%s', FolderName, 'Time.dat');
    fileZ = sprintf('%s/%s', FolderName, 'Z.dat');
    
    strParam = 'Ne = %d\nTend = %f\nLz = %f\nDelta = %f\nI = %f\ndz = %f\ndt = %f\n';
    fileID = fopen(fileParameters ,'wt');
    fprintf(fileID, strParam, RES.Ne, RES.Tend, RES.Lz, RES.Delta, RES.I, RES.dz, RES.dt);
    fclose(fileID);    
    save(fileResultsBvsZ, 'OUTBvsZ', '-ascii');
    save(fileResultsJvsZ, 'OUTJvsZ', '-ascii');
    save(fileResultsBvsT, 'OUTBvsT', '-ascii');
    save(fileResultsJvsT, 'OUTJvsT', '-ascii');
    save(fileT, 'TAxis', '-ascii');
    save(fileZ, 'ZAxis', '-ascii');
end

return

end

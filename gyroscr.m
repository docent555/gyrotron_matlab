function OUT = gyroscr(IN) %#codegen

fprintf('\nTime')

ZAxis = zeros(length(IN.ZAxis),1);
TAxis = zeros(length(IN.TAxis),1);
InitialField = complex(zeros(length(IN.InitialField),1));

ZAxis(:,1) = IN.ZAxis;
TAxis(:,1) = IN.TAxis;
INTERVALT = IN.INTERVALT;
INTERVALZ = IN.INTERVALZ;
Delta = IN.Delta;
I = IN.I;
Ne = IN.Ne;
InitialField(:,1) = IN.InitialField(:,1);
tol = IN.tol;

if INTERVALZ > 1
    IZ = 0:INTERVALZ:length(ZAxis);
    IZ(1) = 1;
    SIZEZ = length(IZ);
else
    IZ = 1:INTERVALZ:length(ZAxis);
end

if INTERVALT > 1 && INTERVALZ > 1
    TAxisNew = zeros(fix(length(TAxis(:,1))/INTERVALT)+1,1);
    OUT.OUTJ = complex(zeros(SIZEZ,fix(length(TAxis(:,1))/INTERVALT)+1));
    OUT.OUTB = complex(zeros(SIZEZ,fix(length(TAxis(:,1))/INTERVALT)+1));
elseif INTERVALT ==1 && INTERVALZ > 1
    TAxisNew = zeros(fix(length(TAxis(:,1))/INTERVALT),1);
    OUT.OUTJ = complex(SIZEZ,fix(length(TAxis(:,1))/INTERVALT));
    OUT.OUTB = complex(SIZEZ,fix(length(TAxis(:,1))/INTERVALT));
elseif INTERVALT > 1 && INTERVALZ == 1
    TAxisNew = zeros(fix(length(TAxis(:,1))/INTERVALT)+1,1);
    OUT.OUTJ = complex(zeros(size(InitialField,1),fix(length(TAxis(:,1))/INTERVALT)+1));
    OUT.OUTB = complex(zeros(size(InitialField,1),fix(length(TAxis(:,1))/INTERVALT)+1));
else % (INTERVALT == 1) && (INTERVALZ == 1)
    TAxisNew = zeros(fix(length(TAxis(:,1))),1);
    OUT.OUTJ = complex(zeros(size(InitialField,1),length(TAxis(:,1))));
    OUT.OUTB = complex(zeros(size(InitialField,1),length(TAxis(:,1))));
end

Field = complex(zeros(size(InitialField,1),1));
OUT.ZAxis = zeros(length(IZ),1);
OUT.TAxis = zeros(length(TAxisNew),1);

% step = 1;
Field(:,1) = InitialField;
TAxisNew(1,1) = TAxis(1,1);

kpar2 = zeros(length(ZAxis),1);

DeltaZ = ZAxis(2) - ZAxis(1);
DeltaT = TAxis(2) - TAxis(1);

% SQR2 = sqrt(2);
SQR2M2 = 2.828427124746190;
SQR2D2 = 0.707106781186548;
SQRDT = sqrt(DeltaT);
SQRDZ = DeltaZ*DeltaZ;

C0 = -1i;
% C0 = 1i;
CR = 0;

WNz = -((0.666666666666667*C0*DeltaZ/DeltaT + kpar2(end)*DeltaZ/3) - 1/DeltaZ);
WNzm1 = -((C0/3*DeltaZ/DeltaT + kpar2(end-1)*DeltaZ/6) + 1/DeltaZ);

C2 = 1/sqrt(1i*pi);
% C2 = 1/sqrt(-1i*pi);

N = length(ZAxis);
A = complex(zeros(N,1));
B = complex(zeros(N-1,1));
C = complex(zeros(N-1,1));
D = complex(zeros(N,1));

A(1) = 1;
A(2:end-1) = -2*(1 - DeltaZ/DeltaT*C0*DeltaZ - DeltaZ*kpar2(2:end-1)*DeltaZ/2);
A(end) = 1 + 1.333333333333333*C2*WNz*SQRDT;

B(1) = 0;
B(2:end) = 1;

C(1:end-1) = 1;
C(end) = 1.333333333333333*C2*WNzm1*SQRDT;

M = spdiags([[C; 0] A [0 ;B]], -1:1, N, N);

WR = complex(zeros(length(TAxis)+1,1));
FNz = complex(zeros(length(TAxis)+1,1));
FNzm1 = complex(zeros(length(TAxis)+1,1));
JNz = complex(zeros(length(TAxis)+1,1));
JNzm1 = complex(zeros(length(TAxis)+1,1));
SigmaNz = complex(zeros(length(TAxis)+1,1));
SigmaNzm1 = complex(zeros(length(TAxis)+1,1));
steps = length(TAxis) - 1;
fmax = zeros(length(TAxis)+1, 1);
field = complex(zeros(length(Field),1));
field_p = complex(zeros(length(Field),1));
rfield_p = complex(zeros(length(Field),1));
lfield_p = complex(zeros(length(Field),1));
NewfieldRaw = complex(zeros(length(Field),1));
rNewfieldRaw = complex(zeros(length(Field),1));
lNewfieldRaw = complex(zeros(length(Field),1));
IROldPart = complex(zeros(1,1));
OldOldJ = complex(zeros(length(Field),1));
Jp = zeros(length(ZAxis),1);
J = zeros(length(ZAxis),1);

jout = 1;
OUT.OUTB(:, jout) = Field(IZ,1);
th0 = 2*pi*(0:Ne-1)/Ne;
theta = zeros(length(ZAxis), Ne);

% Initial values
field(:,1) = Field(:,1);
% J = zeros(length(ZAxis),1);
theta(:,:) = pendulumODE(field(:,1), ZAxis(:,1), th0(1,:), Delta);
J(:,1) = I * 2 / Ne * sum(exp(-1i*theta), 2);

IDX = @(j) (j + 1);

fmax(IDX(0)) = max(abs(field(:,1)));
FNz(IDX(0)) = field(end);
FNzm1(IDX(0)) = field(end-1);
JNz(IDX(0)) = J(end);
JNzm1(IDX(0)) = J(end-1);
SigmaNz(IDX(0)) = 0;
SigmaNzm1(IDX(0)) = 0;

WR(IDX(0)) = DeltaZ * (0.166666666666667*(2 * JNz(IDX(0)) + 2 * JNzm1(IDX(0))));

SHOW = 1;
if SHOW == 1
    [lhr, lha, hFig] = makeFig(ZAxis, TAxis);
end

fprintf('\n');
timerVal = tic;
for step=1:steps 
    
    if SHOW == 1
        lHandleB.YData = abs(field(:,1));
        lHandleJ.YData = abs(J(:,1));
        lhr.YData(1:step) = fmax(1:step);
        lhr.XData(1:step) = TAxis(1:step);
        lha.YData = abs(field);
        drawnow
    end        
    
%     C2mSQRDTm4d3 = C2 * SQRDT * 4 / 3;
%     SQRDZd2 = SQRDZ / 2;
%     C0mSQRDZdDeltaT = C0 * SQRDZ / DeltaT;
%     SQRDTm4d3 = 4 / 3 * SQRDT;
%     SQRDTm2d3 = 2 / 3 * SQRDT;
%     DeltaZmC0m2d3dDeltaT = DeltaZ / DeltaT * C0 * 2 / 3;
%     DeltaZmC0d3dDeltaT = DeltaZ / DeltaT * C0 / 3;
%     DeltaZd2 = DeltaZ / 2;
%     DeltaZd3 = DeltaZ / 3;
%     DeltaZd6 = DeltaZ / 6;
    
    if (step ~= 1) && (mod(step-1,INTERVALT) == 0)
        OUT.OUTJ(:,jout) = J(IZ,1);
    end
    
    WR(IDX(step)) = DeltaZ * ((C0 * 0.666666666666667 / DeltaT - kpar2(end) / 3) * FNz(IDX(step-1))...
        + (C0 / 3 / DeltaT - kpar2(end - 1) / 6) * FNzm1(IDX(step-1))...
        + 0.166666666666667*(4 * JNz(IDX(step-1)) + 2 * JNzm1(IDX(step-1))) - (2 * SigmaNz(IDX(step-1)) + SigmaNzm1(IDX(step-1))));

    u = @(j) (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j))).' .* exp(CR * DeltaT * (step - j));
    
    if step == 1
        IR = 0;
    elseif step == 2
        IR = 1.333333333333333 * SQRDT * (u(0)*(1 - SQR2D2) + u(1)*(SQR2M2 - 2.5));
    else
        j = 1:step-2;
        IR = 1.333333333333333 * SQRDT * (u(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
            + sum(u(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
            + u(step - 1)*(SQR2M2 - 2.5));
    end    
    
    D(1) = 0;
    %         D(1) = IN.TimeAmp * exp(1i * IN.TimeFreq * AxisTau(step));
    D(2:end - 1) = DeltaZ.^2 * (2*J(2:end-1)) ...
        + 2 * (1 + C0 * DeltaZ.^2/DeltaT - DeltaZ.^2 * kpar2(2:end-1)/2) .* field(2:end - 1)...
        - (field(1:end - 2) + field(3:end));
    D(end) = -C2*(IR + 1.333333333333333 * WR(IDX(step)) * SQRDT + 0.666666666666667 * SQRDT * (WNzm1 * field(end - 1)...
        + WNz * field(end) + WR(IDX(step-1))) * exp(CR * DeltaT));
    
    % nesamosoglasovannoe pole
    field_p = M \ D;
    %     rfield_p = rtridag(C,A,B,D);
    %     lfield_p = ltridag(C,A,B,D);
    %     field_p = (rfield_p + lfield_p)/2;
    
    while 1        
%         Jp = zeros(length(ZAxis),1);
        theta(:,:) = pendulumODE(field_p(:,1), ZAxis(:,1), th0(1,:), Delta);
        Jp(:,1) = I * 2 / Ne * sum(exp(-1i*theta), 2);               
        
        WR(IDX(step)) = DeltaZ * ((C0 * 0.666666666666667 / DeltaT - kpar2(end) / 3) * field(end)...
            + (C0 / 3 / DeltaT - kpar2(end - 1) / 6)*field(end - 1)...
            + 0.166666666666667 * (2 * Jp(end) + 2 * J(end) + Jp(end - 1) + J(end - 1)) - (2 * SigmaNz(IDX(step-1)) + SigmaNzm1(IDX(step-1))));        
        
        D(2:end - 1) = DeltaZ.^2 * (Jp(2:end - 1) + J(2:end - 1)) ...
            + 2 * (1 + C0 * DeltaZ.^2 / DeltaT - DeltaZ.^2 * kpar2(2:end - 1) / 2).*field(2:end - 1)...
            - (field(1:end - 2) + field(3:end));
        D(end) = -C2*(IR + 1.333333333333333 * WR(IDX(step)) * SQRDT + 0.666666666666667 * SQRDT * (WNzm1 * field(end - 1)...
            + WNz * field(end) + WR(IDX(step-1))) * exp(CR * DeltaT));
               
        
        % samosoglasovannoe pole
        field_p(:,1) = M \ D;
        %     rfield_p(:,1) = rtridag(C,A,B,D);
        %     lfield_p(:,1) = rtridag(C,A,B,D);
        %     field_p = (rNewfieldRaw + lNewfieldRaw)/2;
        
        
        fmax(IDX(step)) = max(abs(field(:,1)));
        
        if ((fmax(IDX(step)) - fmax(IDX(step-1)))/fmax(IDX(step))) < tol
            break
        end
    end
    
    field(:,1) = field_p(:,1);
%     J = zeros(length(ZAxis),1);
    theta(:,:) = pendulumODE(field(:,1), ZAxis(:,1), th0(1,:), Delta);
    J(:,1) = I * 2 / Ne * sum(exp(-1i*theta), 2);
    
    k = step + 1;
    
    if mod(step,INTERVALT) == 0
        jout = jout + 1;
        OUT.OUTB(:, jout) = field(IZ,1);
        TAxisNew(jout,1) = TAxis(k,1);
    end
    
    FNz(IDX(step)) =  field(end);
    FNzm1(IDX(step)) = field(end - 1);
    JNz(IDX(step)) = J(end);
    JNzm1(IDX(step)) = J(end - 1); 
    
    SigmaNz(IDX(step)) = -(kpar2(end)/6 + C0/3/DeltaT) * FNz(IDX(step)) ...
        + (C0/3/DeltaT - kpar2(end)/6) * FNz(IDX(step - 1)) ...
        + 0.166666666666667*(JNz(IDX(step)) + JNz(IDX(step - 1))) - SigmaNz(IDX(step - 1));
    SigmaNzm1(IDX(step)) = -(kpar2(end - 1)/6 + C0/3/DeltaT) * FNzm1(IDX(step)) ...
        + (C0/3/DeltaT - kpar2(end - 1)/6) * FNzm1(IDX(step - 1)) ...
        + 0.166666666666667*(JNzm1(IDX(step)) + JNzm1(IDX(step - 1))) - SigmaNzm1(IDX(step - 1));
    
    fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        'Step = %8i   Time = %8.4f   Bmax = %15.10f   Jmax = %15.10f'],...
        step, TAxis(k), fmax(k), max(abs(J(:,1))));
end

OUT.OUTJ(:,jout) = J(IZ,1); % для последнего jout см. стр.212

fprintf("\n\n\n");

ExecutionTime = toc(timerVal);

hours = fix(ExecutionTime/3600);
minutes = fix((ExecutionTime - fix(hours*3600))/60);
seconds = ExecutionTime - hours*3600 - minutes*60;

fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

OUT.ZAxis(:,1) = ZAxis(IZ,1);
OUT.TAxis(:,1) = TAxisNew(:,1);

fprintf(" \n\n");

% if SHOW == 1
%     close(hFig);
% end

end


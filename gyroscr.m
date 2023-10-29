function [OUTB, OUTJ] = gyroscr(Nz, Nt, Ne, ZAxis, TAxis, Delta, Ic, dt, dz, tol, INTT, INTZ, OUTNz, OUTNt, InitialField) %#codegen
% function OUT = gyroscr(IN) %#codegen

% fprintf('\nTime')

Field = complex(zeros(size(InitialField,1),1));
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
cu_p = zeros(length(ZAxis),1);
cu = zeros(length(ZAxis),1);
OUTB = zeros(OUTNz, OUTNt);
OUTJ = zeros(OUTNz, OUTNt);

if INTZ > 1
    IZ = 0:INTZ:length(ZAxis);
    IZ(1) = 1;
    SIZEZ = length(IZ);
else
    IZ = 1:INTZ:length(ZAxis);
end

Field(:,1) = InitialField;
kpar2 = zeros(length(ZAxis),1);
% N = length(ZAxis);
A = complex(zeros(Nz,1));
B = complex(zeros(Nz-1,1));
C = complex(zeros(Nz-1,1));
D = complex(zeros(Nz,1));

% SQR2 = sqrt(2.0D0);
SQR2M2 = 2.828427124746190;
SQR2D2 = 0.707106781186548;
SQRDT = sqrt(dt);
SQRDZ = dz*dz;

C0 = -1i;
% C0 = 1i;
CR = 0;
C2 = 1/sqrt(1i*pi);
% C2 = 1/sqrt(-1i*pi);
WNz = -((2.0D0/3.0D0*C0*dz/dt + kpar2(end)*dz/3.0D0) - 1.0D0/dz);
WNzm1 = -((C0/3.0D0*dz/dt + kpar2(end-1)*dz/6.0D0) + 1.0D0/dz);


A(1) = 1.0D0;
A(2:end-1) = -2.0D0*(1.0D0 - dz/dt*C0*dz - dz*kpar2(2:end-1)*dz/2.0D0);
A(end) = 1.0D0 + 4.0D0/3.0D0*C2*WNz*SQRDT;
B(1) = 0;
B(2:end) = 1.0D0;
C(1:end-1) = 1.0D0;
C(end) = 4.0D0/3.0D0*C2*WNzm1*SQRDT;

% M = spdiags([[C; 0] A [0 ;B]], -1:1, N, N);

jout = 1;
OUTB(:, jout) = Field(IZ,1);
th0 = 2.0D0*pi*(0:Ne-1)/Ne;
theta = zeros(length(ZAxis), Ne);

% Initial values
field(:,1) = Field(:,1);
theta(:,:) = pendulumODE(field(:,1), ZAxis(:,1), th0(1,:), Delta);
cu(:,1) = Ic * 2.0D0 / Ne * sum(exp(-1i*theta), 2);
OUTJ(:,jout) = cu(IZ,1);

IDX = @(j) (j + 1);

fmax(IDX(0)) = max(abs(field(:,1)));
FNz(IDX(0)) = field(end);
FNzm1(IDX(0)) = field(end-1);
JNz(IDX(0)) = cu(end);
JNzm1(IDX(0)) = cu(end-1);
SigmaNz(IDX(0)) = 0;
SigmaNzm1(IDX(0)) = 0;

WR(IDX(0)) = dz * (2.0D0/3.0D0*(2.0D0 * JNz(IDX(0)) + JNzm1(IDX(0))));

SHOW = 1;
if SHOW == 1
    [lhr, lha, hFig] = makeFig(ZAxis, TAxis);
end

fprintf('\n');
timerVal = tic;
for step=1:steps 
    
    if SHOW == 1
        lHandleB.YData = abs(field(:,1));
        lHandleJ.YData = abs(cu(:,1));
        lhr.YData(1:step) = fmax(1:step);
        lhr.XData(1:step) = TAxis(1:step);
        lha.YData = abs(field);
        drawnow
    end        
    
%     C2mSQRDTm4d3 = C2 * SQRDT * 4.0D0 / 3.0D0;
%     SQRDZd2 = SQRDZ / 2.0D0;
%     C0mSQRDZdDeltaT = C0 * SQRDZ / dt;
%     SQRDTm4d3 = 4.0D0 / 3.0D0 * SQRDT;
%     SQRDTm2d3 = 2.0D0 / 3.0D0 * SQRDT;
%     DeltaZmC0m2d3dDeltaT = dz / dt * C0 * 2.0D0 / 3.0D0;
%     DeltaZmC0d3dDeltaT = dz / dt * C0 / 3.0D0;
%     DeltaZd2 = dz / 2.0D0;
%     DeltaZd3 = dz / 3.0D0;
%     DeltaZd6 = dz / 6.0D0;
 
    WR(IDX(step)) = dz * ((C0 * 2.0D0/3.0D0/dt - kpar2(end)/3.0D0) * FNz(IDX(step-1))...
        + (C0/3.0D0/dt - kpar2(end - 1)/6.0D0) * FNzm1(IDX(step-1))...
        + 1.0D0/6.0D0*(4.0D0 * JNz(IDX(step-1)) + 2.0D0 * JNzm1(IDX(step-1))) - (2.0D0 * SigmaNz(IDX(step-1)) + SigmaNzm1(IDX(step-1))));

    u = @(j) (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j))).' .* exp(CR * dt * (step - j));
    
    if step == 1
        IR = 0;
    elseif step == 2
        IR = 4.0D0/3.0D0 * SQRDT * (u(0)*(1 - SQR2D2) + u(1)*(SQR2M2 - 2.5D0));
    else
%         j = 1:step-2;
%         IR = 4.0D0/3.0D0 * SQRDT * (u(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
%             + sum(u(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
%             + u(step - 1)*(SQR2M2 - 2.5));
        IR = 4.0D0/3.0D0 * SQRDT * (u(0)*((step - 1.0D0).^(1.5) - (step - 1.5D0)*sqrt(step)) + u(step - 1)*(SQR2M2 - 2.5D0));
        for j = 1:step-2
            IR = IR + 4.0D0/3.0D0 * SQRDT * (u(j).*((step - j - 1.0D0).^(1.5) - 2.0D0*(step - j).^(1.5) + (step - j + 1.0D0).^(1.5)));
        end
    end    
    
    D(1) = 0;
    %         D(1) = IN.TimeAmp * exp(1i * IN.TimeFreq * AxisTau(step));
    D(2:end - 1) = SQRDZ*(2.0D0*cu(2:end-1)) ...
        + 2.0D0 * (1.0D0 + C0 * SQRDZ/dt - SQRDZ * kpar2(2:end-1)/2.0D0) .* field(2:end - 1)...
        - (field(1:end - 2) + field(3:end));
    D(end) = -C2 *(IR + 4.0D0/3.0D0 * WR(IDX(step)) * SQRDT + 2.0D0/3.0D0 * SQRDT * (WNzm1 * field(end - 1)...
        + WNz * field(end) + WR(IDX(step-1))) * exp(CR * dt));
    
    % nesamosoglasovannoe pole
    %     field_p = M \ D;
    rfield_p = rtridag(C,A,B,D);
    lfield_p = ltridag(C,A,B,D);
    field_p = (rfield_p + lfield_p)/2.0D0;        
    
    maxfield = max(abs(field(:,1)));    
    while 1        
%         cu_p = zeros(length(ZAxis),1);
        theta(:,:) = pendulumODE(field_p(:,1), ZAxis(:,1), th0(1,:), Delta);
        cu_p(:,1) = Ic * 2.0D0 / Ne * sum(exp(-1i*theta), 2);               
        
        WR(IDX(step)) = dz * ((C0 * 2.0D0/3.0D0 / dt - kpar2(end) / 3.0D0) * field(end)...
            + (C0 / 3.0D0 / dt - kpar2(end - 1) / 6.0D0)*field(end - 1)...
            + 1.0D0/6.0D0 * (2.0D0 * cu_p(end) + 2.0D0 * cu(end) + cu_p(end - 1) + cu(end - 1)) - (2.0D0 * SigmaNz(IDX(step-1)) + SigmaNzm1(IDX(step-1))));        
        
        D(2:end - 1) = SQRDZ * (cu_p(2:end - 1) + cu(2:end - 1)) ...
            + 2.0D0 * (1 + C0 * SQRDZ / dt - SQRDZ * kpar2(2:end - 1) / 2.0D0).*field(2:end - 1)...
            - (field(1:end - 2) + field(3:end));
        D(end) = -C2*(IR + 4.0D0/3.0D0 * WR(IDX(step)) * SQRDT + 2.0D0/3.0D0 * SQRDT * (WNzm1 * field(end - 1)...
            + WNz * field(end) + WR(IDX(step-1))) * exp(CR * dt));
               
        
        % samosoglasovannoe pole
        %field_p(:,1) = M \ D;
        rfield_p(:,1) = rtridag(C,A,B,D);
        lfield_p(:,1) = rtridag(C,A,B,D);
        field_p = (rfield_p + lfield_p)/2.0D0;
        
        
        maxfield_p = max(abs(field_p(:,1)));
        maxdiff = abs(maxfield - maxfield_p)/maxfield;
        if maxdiff < tol
            break
        end
        maxfield = maxfield_p;
    end
    
    field(:,1) = field_p(:,1);
    theta(:,:) = pendulumODE(field(:,1), ZAxis(:,1), th0(1,:), Delta);
    cu(:,1) = Ic * 2.0D0 / Ne * sum(exp(-1i*theta), 2);
    fmax(IDX(step)) = max(abs(field(:,1)));
    
    k = step + 1;
    
    if mod(step,INTT) == 0
        jout = jout + 1;
        OUTB(:, jout) = field(IZ,1);
        OUTJ(:, jout) = cu(IZ,1);
    end
    
    FNz(IDX(step)) =  field(end);
    FNzm1(IDX(step)) = field(end - 1);
    JNz(IDX(step)) = cu(end);
    JNzm1(IDX(step)) = cu(end - 1); 
    
    SigmaNz(IDX(step)) = -(kpar2(end)/6.0D0 + C0/3.0D0/dt) * FNz(IDX(step)) ...
        + (C0/3.0D0/dt - kpar2(end)/6.0D0) * FNz(IDX(step - 1)) ...
        + 1.0D0/6.0D0*(JNz(IDX(step)) + JNz(IDX(step - 1))) - SigmaNz(IDX(step - 1));
    SigmaNzm1(IDX(step)) = -(kpar2(end - 1)/6.0D0 + C0/3.0D0/dt) * FNzm1(IDX(step)) ...
        + (C0/3.0D0/dt - kpar2(end - 1)/6.0D0) * FNzm1(IDX(step - 1)) ...
        + 1.0D0/6.0D0*(JNzm1(IDX(step)) + JNzm1(IDX(step - 1))) - SigmaNzm1(IDX(step - 1));
    
%     fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
%         '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
%         '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
%         '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
%         'Step = %8i   Time = %8.4f   Bmax = %15.10f   Jmax = %15.10f'],...
%         step, TAxis(k), fmax(k), max(abs(cu(:,1))));
    
    str = sprintf('Step = %8i   Time = %8.4f   Bmax = %15.10f   Jmax = %15.10f', step, TAxis(k), fmax(k), max(abs(cu(:,1))));
    
    disp(str);
    
end

OUTJ(:,jout) = cu(IZ,1);

fprintf("\n\n\n");

ExecutionTime = toc(timerVal);

hours = fix(ExecutionTime/3600);
minutes = fix((ExecutionTime - fix(hours*3600))/60);
seconds = ExecutionTime - hours*3600 - minutes*60;

fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

fprintf(" \n\n");
end



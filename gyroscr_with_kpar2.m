function [OUTF, OUTJ, Eff, Omega, jout] = gyroscr(Nz, Nt, Ne, ZAxis, TAxis, Delta, Ic, dt, dz, tol, kpar2, INTT, INTZ, OUTNz, OUTNt, InitialField) %#codegen

WR = complex(zeros(Nt,1));
FNz = complex(zeros(Nt,1));
FNzm1 = complex(zeros(Nt,1));
JNz = complex(zeros(Nt,1));
JNzm1 = complex(zeros(Nt,1));
SigmaNz = complex(zeros(Nt,1));
SigmaNzm1 = complex(zeros(Nt,1));
fmax = zeros(Nt, 1);
jmax = zeros(Nt, 1);
field = complex(zeros(Nz,1));
field_p = complex(zeros(Nz,1));
rfield_p = complex(zeros(Nz,1));
lfield_p = complex(zeros(Nz,1));
cu_p = complex(zeros(Nz,1));
cu = complex(zeros(Nz,1));
J_p = complex(zeros(Nz,1));
J = complex(zeros(Nz,1));
OUTF = complex(zeros(OUTNz, OUTNt));
OUTJ = complex(zeros(OUTNz, OUTNt));
Eff = complex(zeros(Nt,1));
Omega = complex(zeros(Nt,1));
% theta = zeros(Nz, Ne);
% p = zeros(Nz, Ne);
% pv = zeros(Nz, 2*Ne);
% p0 = zeros(Ne,1);
% p0v = zeros(2*Ne,1);
% reidx = zeros(1,Ne);
% imidx = zeros(1,Ne);

if INTZ > 1
    IZ = 0:INTZ:length(ZAxis);
    IZ(1) = 1;
    SIZEZ = length(IZ);    
else
    IZ = 1:INTZ:length(ZAxis);    
end

% kpar2 = zeros(length(ZAxis),1);
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

Nzm1 = Nz - 1;
C0 = 1.0D0;
CR = 0;
C2 = 1.0D0/sqrt(1i*pi);
WNz = -((-1i*2.0D0/3.0D0*C0*dz/dt + kpar2(Nz)*dz/3.0D0) - 1.0D0/dz);
WNzm1 = -((-1i*C0/3.0D0*dz/dt + kpar2(Nzm1)*dz/6.0D0) + 1.0D0/dz);

A(1) = 1.0D0;
% A(2:Nzm1) = -2.0D0*(1.0D0 + 1i * dz/dt*C0*dz - dz*kpar2(2:Nzm1)*dz/2.0D0);
A(2:Nzm1) = -2.0D0*(1.0D0 + 1i * SQRDZ/dt*C0 - dz*kpar2(2:Nzm1)*dz/2.0D0);
A(Nz) = 1.0D0 + 4.0D0/3.0D0*C2*WNz*SQRDT;
B(1) = 0;
B(2:Nzm1) = 1.0D0;
C(1:Nz-2) = 1.0D0;
C(Nzm1) = 4.0D0/3.0D0*C2*WNzm1*SQRDT;

M = spdiags([[C; 0] A [0 ;B]], -1:1, Nz, Nz);

% Initial values
jout = 1;
field(:,1) = InitialField;
OUTF(:, jout) = field(IZ,1);
th0 = 2.0D0*pi*(0:Ne-1)/Ne;
p0 = exp(1i*th0)';
p0v = [real(p0); imag(p0)];
reidx = 1:Ne;
imidx = Ne+1:2*Ne;
p = oscill_reim(field, Nz, ZAxis, Delta, p0v, reidx, imidx);
% p = oscill_cmplx(field, ZAxis, Delta, p0);
J(:,1) = Ic * trapz(th0, p, 2)  / (2*pi);
% J(:,1) = Ic * trpz(dz, p, Ne)  / (2*pi);
cu(:,1) = J(:) - 1i*kpar2(:).*field(:);
OUTJ(:,jout) = J(IZ,1);


IDX = @(j) (j + 1);

fmax(IDX(0)) = max(abs(field(:,1)));
jmax(IDX(0)) = max(abs(cu(:,1)));
FNz(IDX(0)) = field(Nz);
FNzm1(IDX(0)) = field(Nzm1);
JNz(IDX(0)) = cu(Nz);
JNzm1(IDX(0)) = cu(Nzm1);
SigmaNz(IDX(0)) = 0;
SigmaNzm1(IDX(0)) = 0;
Eff(IDX(0)) = 1 - trapz(th0, abs(p(Nz,:).^2))/(2*pi);
Omega(IDX(0)) = 0;

WR(IDX(0)) = dz * (2.0D0/3.0D0*(2.0D0 * JNz(IDX(0)) + JNzm1(IDX(0))));

SHOW = 1;
if SHOW == 1
    [lhfmax, lhfabs, lhjmax, lhjabs, hFig] = makeFig(ZAxis, TAxis);
end

%Coefficients
coeff_1i_m_C0_m_2_d_3_d_dt = 1i*C0*2.0D0/3.0D0/dt;
coeff_1i_m_C0_d_3_d_dt = 1i*C0/3.0D0/dt;
coeff_1i_d_6 = 1i/6.0D0;
coeff_4_d_3_m_SQRDT = 4.0D0/3.0D0*SQRDT;
coeff_2_d_3_m_SQRDT = 2.0D0/3.0D0*SQRDT;
coeff_1i_m_SQRDZ = 1i*SQRDZ;
coeff_1i_m_C0_m_SQRDZ_d_dt = 1i*C0*SQRDZ/dt;
coeff_C2_m_coeff_4_d_3_m_SQRDT = C2*coeff_4_d_3_m_SQRDT;
coeff_exp_CR_m_dt = exp(CR*dt);
coeff_CR_m_dt = CR*dt;
coeff_dz_m_coeff_1i_d_6 = dz*coeff_1i_d_6;

num_st_test_iter = 0;
fmax_glob_old = max(abs(field(:,1)));

fprintf('\n');
timerVal = tic;
for step=1:Nt-1
    
        if SHOW == 1           
            lhfmax.YData(1:step) = Omega(1:step);
            lhfmax.XData(1:step) = TAxis(1:step);
            lhfabs.YData = abs(field);
            
            lhjmax.YData(1:step) = Eff(1:step);
            lhjmax.XData(1:step) = TAxis(1:step);
            lhjabs.YData = abs(J);
            
            drawnow
        end            
    
    WR_PART = dz * ((-coeff_1i_m_C0_m_2_d_3_d_dt - kpar2(Nz) / 3.0D0) * field(Nz)...
            + (-coeff_1i_m_C0_d_3_d_dt - kpar2(Nzm1) / 6.0D0)*field(Nzm1)...
            - (2.0D0 * SigmaNz(IDX(step-1)) + SigmaNzm1(IDX(step-1))));
        
    WR(IDX(step)) = -coeff_dz_m_coeff_1i_d_6*(4.0D0 * cu(Nz) + 2.0D0 * cu(Nzm1)) + WR_PART;
                      
    
%     u = @(j) (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j))).' .* exp(CR*dt * (step - j));
    
    if step == 1
        IR = 0;
    elseif step == 2
        IR = coeff_4_d_3_m_SQRDT * (u(0)*(1 - SQR2D2) + u(1)*(SQR2M2 - 2.5D0));
    else
        j = 1:step-2;
        IR = coeff_4_d_3_m_SQRDT * (u(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
            + sum(u(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
            + u(step - 1)*(SQR2M2 - 2.5));
%         IR = coeff_4_d_3_m_SQRDT * (u(0)*((step - 1.0D0).^(1.5) - (step - 1.5D0)*sqrt(step)) + u(step - 1)*(SQR2M2 - 2.5D0));
%         for j = 1:step-2
%             IR = IR + coeff_4_d_3_m_SQRDT * (u(j).*((step - j - 1.0D0).^(1.5) - 2.0D0*(step - j).^(1.5) + (step - j + 1.0D0).^(1.5)));
%         end
    end
    
    D(1) = 0;
    %         D(1) = IN.TimeAmp * exp(1i * IN.TimeFreq * AxisTau(step));
    
    D_MIDDLE_PART = 2.0D0 * (1.0D0 - coeff_1i_m_C0_m_SQRDZ_d_dt - SQRDZ/2.0D0 * kpar2(2:Nzm1)) .* field(2:Nzm1)...
        - (field(1:Nz - 2) + field(3:Nz));
    
    D(2:Nzm1) = -coeff_1i_m_SQRDZ * (2.0D0*cu(2:Nzm1)) + D_MIDDLE_PART;
    
    D_END_PART = -C2 * (IR + coeff_2_d_3_m_SQRDT * (WNzm1 * field(Nzm1)...
        + WNz * field(Nz) + WR(IDX(step-1))) * coeff_exp_CR_m_dt);
                
    D(Nz) = -coeff_C2_m_coeff_4_d_3_m_SQRDT * WR(IDX(step)) + D_END_PART;
    
    % nesamosoglasovannoe pole
    field_p = M \ D;
%     rfield_p = rtridag(C,A,B,D);
%     lfield_p = ltridag(C,A,B,D);
%     field_p = (rfield_p + lfield_p)/2.0D0;
    
    num_insteps = 0;
    maxfield = max(abs(field_p(:,1)));
    while 1
        num_insteps = num_insteps + 1;
        p = oscill_reim(field_p, Nz, ZAxis, Delta, p0v, reidx, imidx);
%         p = oscill_cmplx(field_p, ZAxis, Delta, p0);
        J_p(:,1) = Ic * trapz(th0, p, 2)  / (2.0D0*pi);
        %         J_p(:,1) = Ic * trpz(dz, p, Ne)  / (2*pi);
        cu_p(:,1) = J_p(:) - 1i*kpar2(:).*field_p(:);
        
        WR(IDX(step)) = -coeff_dz_m_coeff_1i_d_6 * (2.0D0 * cu_p(Nz) + 2.0D0 * cu(Nz) + cu_p(Nzm1) + cu(Nzm1)) + WR_PART;
        
        D(2:Nzm1) = -coeff_1i_m_SQRDZ * (cu_p(2:Nzm1) + cu(2:Nzm1)) + D_MIDDLE_PART;            
        
        D(Nz) = -coeff_C2_m_coeff_4_d_3_m_SQRDT * WR(IDX(step)) + D_END_PART;
                
        % samosoglasovannoe pole
        field_p(:,1) = M \ D;
%         rfield_p(:,1) = rtridag(C,A,B,D);
%         lfield_p(:,1) = ltridag(C,A,B,D);
%         field_p = (rfield_p + lfield_p)/2.0D0;
        
        
        maxfield_p = max(abs(field_p(:,1)));
        maxdiff = abs(maxfield - maxfield_p)/maxfield;
        if maxdiff < tol
            break
        end
        maxfield = maxfield_p; 
        if num_insteps > 3000
            error('Too many inner steps!');
        end
    end
    
    field(:,1) = field_p(:,1);
    p = oscill_reim(field, Nz, ZAxis, Delta, p0v, reidx, imidx);
%     p = oscill_cmplx(field, ZAxis, Delta, p0);
    J(:,1) = Ic * trapz(th0, p, 2)  / (2.0D0*pi);
    %     J(:,1) = Ic * trpz(dz, p, Ne)  / (2*pi);
    cu(:,1) = J(:) - 1i*kpar2(:).*field(:);
    fmax(IDX(step)) = max(abs(field(:,1)));
    jmax(IDX(step)) = max(abs(cu(:,1)));    
    
    FNz(IDX(step)) =  field(Nz);
    FNzm1(IDX(step)) = field(Nzm1);
    JNz(IDX(step)) = cu(Nz);
    JNzm1(IDX(step)) = cu(Nzm1);
    
    Omega(IDX(step)) = (angle(field(Nz)) - angle(FNz(IDX(step-1))))/dt;
    Eff(IDX(step)) = 1 - trapz(th0, abs(p(Nz,:).^2))/(2*pi);
            
    if (mod(num_st_test_iter,1000))
        fmax_glob_new = max(abs(field(:,1)));
        if abs(fmax_glob_new - fmax_glob_old)/fmax_glob_old < tol
            jout = jout + 1;
            OUTF(:, jout) = field(IZ,1);
            OUTJ(:, jout) = J(IZ,1);
            return;
        end
        num_st_test_iter = num_st_test_iter + 1;
        fmax_glob_old = fmax_glob_new;
    end
    
    if mod(step,INTT) == 0
        jout = jout + 1;
        OUTF(:, jout) = field(IZ,1);
        OUTJ(:, jout) = J(IZ,1);
    end        
    
    SigmaNz(IDX(step)) = -(kpar2(Nz)/6.0D0 - coeff_1i_m_C0_d_3_d_dt) * field(Nz) ...
        + (-coeff_1i_m_C0_d_3_d_dt - kpar2(Nz)/6.0D0) * FNz(IDX(step - 1)) ...
        - coeff_1i_d_6*(cu(Nz) + JNz(IDX(step - 1))) - SigmaNz(IDX(step - 1));
    SigmaNzm1(IDX(step)) = -(kpar2(Nzm1)/6.0D0 -coeff_1i_m_C0_d_3_d_dt) * field(Nzm1) ...
        + (-coeff_1i_m_C0_d_3_d_dt - kpar2(Nzm1)/6.0D0) * FNzm1(IDX(step - 1)) ...
        - coeff_1i_d_6*(cu(Nzm1) + JNzm1(IDX(step - 1))) - SigmaNzm1(IDX(step - 1));    
    
    k = step + 1;
    
    fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...        
        'Step = %8i   Time = %10.4f   Bmax = %+15.10e   Jmax = %+15.10e   W = %+15.10e   E = %+15.10e'],...
        int64(step), TAxis(k), fmax(k), max(abs(cu(:,1))), Omega(IDX(step)), Eff(IDX(step)));
    
end

OUTJ(:,jout) = J(IZ,1);

fprintf("\n\n\n");

ExecutionTime = toc(timerVal);

hours = fix(ExecutionTime/3600);
minutes = fix((ExecutionTime - fix(hours*3600))/60);
seconds = ExecutionTime - hours*3600 - minutes*60;

fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

fprintf(" \n\n");

    function  f = u(j)
        f = (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j))).' .* exp(coeff_CR_m_dt * (step - j));
    end
end




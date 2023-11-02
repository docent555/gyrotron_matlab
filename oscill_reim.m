function p = oscill_reim(field, Nz, ZAxis, Delta, p0v, reidx, imidx)
fre = real(field);
fim = imag(field);
[reb,rec,red] = spline(Nz,ZAxis,fre);
[imb,imc,imd] = spline(Nz,ZAxis,fim);
S1 = @(z) seval_cmplx(z, Nz, ZAxis, fre, fim, reb, rec, red, imb, imc, imd);
[~, pv] = ode45(@(z, p) rhsv(z, p, Delta, S1, reidx, imidx) , ZAxis , p0v);
p = pv(:,reidx) + 1i*pv(:,imidx);
end


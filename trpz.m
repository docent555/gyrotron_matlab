function  u = trpz(dz, p, n) %%codegen
u = (p(:,1) + p(:,end))/2;
for i=2:n-1
    u = u + p(:,i);
end
u = u*dz;
end


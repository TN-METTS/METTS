
T = 0.5; 
beta = 1/T;
tmax = beta/2;
Nkeep = 30;
Nsweep = 4; 
dt =40; 
N = 100 ;

% Set H as MPO
Hloc = cell(5,5);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = S(:,:,1);
Hloc{3,1} = S(:,:,2);
Hloc{4,1} = S(:,:,3);
Hloc{5,2} = J*S(:,:,1)';
Hloc{5,3} = J*S(:,:,2)';
Hloc{5,4} = J*S(:,:,3)';
Hloc{end,end} = I;
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));

H_MPO = cell(1,N);
H_MPO(:) = {Hloc};
H_MPO{1} = H_MPO{1}(:,:,end,:); % choose the last components of the left leg
H_MPO{end} = H_MPO{end}(:,:,:,1); % choose the first components of the right leg
Hs = H_MPO;
%[H_MPO, lognorm_H_MPO] = MPO_canonForm(H_MPO, 0, [], 1e-8);
H2_MPO = mtimes_MPO (H_MPO, H_MPO, Nkeep, Nsweep);

[taus,lnZs,rho] = XTRG (Hs,dt,tmax,Nkeep,Nsweep);

lnZ = lnZs(end);

Sz = S(:,:,2);
Szloc=cell(2, 2);
Szloc{1, 1}=I;
Szloc{1, 2}=zeros(size(I));
Szloc{2, 1}=Sz;
Szloc{2, 2} = I; 
Szloc = cell2mat(reshape(Szloc, [1, 1, size(Szloc, 1), size(Szloc, 2)]));
Sz_MPO = cell(1, N);
Sz_MPO(:) = {Szloc};
Sz_MPO{1} = Sz_MPO{1}(:,:,end,:); % choose the last components of the left leg
Sz_MPO{end} = Sz_MPO{end}(:,:,:,1); % choose the first components of the right leg
Sz2_MPO = mtimes_MPO(Sz_MPO, Sz_MPO, Nkeep,Nsweep);






E=expVal(Hs, rho, lnZ)
E2=expVal(H2_MPO, rho, lnZ)
m = expVal(Sz_MPO, rho, lnZ)
m2 =expVal(Sz2_MPO, rho, lnZ)


Cv = (E2-E^2)/N/T^2 
chi = (m2-m^2)/N/T


function val = expVal(MPO, rho, lnZ)
    Tr = 1;
    for itN = (1:N)
        Tr = contract(Tr,2,2,rho{itN},4,3);
        Tr = contract(Tr,4,[2 3],MPO{itN},2,[2 1]);
    end
    val = Tr/exp(lnZ)
end 


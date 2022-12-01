function M = CPS_collapseonly(M,dim,s)
% M = CPS_collapse(M,dim,S,n) 
%
% Obtain CPS by collapsing METTS; collapsed state is the eigenstate of spin operator along the axis.
% the direction of axis is determined by input parameter 's'.
%
% 
%
%
% < Input >
% M : [1 x N cell array] Matrix product state to collapse 
%
%    1      2   1      2         1        2
%   ---M{1}---*---M{2}---* ... *---M{end}---
%       |          |                 |
%       ^3         ^3                ^3
% 
% dim : [integer] dimension of local Hilbert space (i.e. the dimension of physical
% leg)
% s : [integer] s = 1 -> only collapse to S_x eigenstates; s = 2 -> S_y; s
% = 3 -> S_z
% 
% Output :
% M : [1 x N cell array] collapsed matrix product space (to classical
% product space)
%

tobj = tic2;


N = numel(M); % the number of sites

C = cell(1,dim); % set for pure state projector

[S,~] = getLocalSpace('Spin',1);

if s == 1
    S = (S(:,:,1) + S(:,:,3))/2;
    axis = 'x';
elseif s == 2
    S = (S(:,:,3) - S(:,:,1))*1i/2;
    axis = 'y';
elseif s == 3
    S = S(:,:,2);
    axis = 'z';
end

[V,~] = eig(S);

for it = (1:dim)
    null = zeros(dim);
    null(:,it) = V(:,it);
    C{it} = null;
end


[M,~,~] = canonForm(M,N,[],0);
[M,~,~] = canonForm(M,1,[],0); % site-canonical form with orthogonality center site 1



for it = (1:N-1)
    nrm = sqrt(contract(M{it},3,[1 2 3],conj(M{it}),3,[1 2 3]));
    M{it} = M{it}/nrm; % normalize MPS
    
    % choose 1 eigenstate randomly
    idx = randi(3);
    M{it} = contract(M{it},3,3,C{idx},2,2);
    
    % move on to the next site (it+1) via SVD
    
    [U,S,Vd] = svdTr(M{it},3,[1 3],[],dim*eps);
    
    M{it} = permute(U,[1 3 2]);
    S = contract(diag(S),2,2,Vd,2,1);
    M{it+1} = contract(S,2,2,M{it+1},3,1);
end

disptime(['Collapsed along the',axis,'axis']);


toc2(tobj,'-v');

end


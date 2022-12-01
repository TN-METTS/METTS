function [ts,M,Ovals,EE,dw] = TS_1D(M,Hs,O,Nkeep,dt,tmax, print_log)
% < Description >
%
% [ts,M,Ovals,EE,dw] = TS_1D(M,Hs,O,Nkeep,dt,tmax, print_log)
%
% < Input >
% M : [cell] The initial state as the MPS. The length of M, i.e., numel(M),
%       defines the chain length. The leg convention of M{n} is as follows:
%
%    1      2   1      2         1        2
%   ---M{1}---*---M{2}---* ... *---M{end}---
%       |          |                 |
%       ^3         ^3                ^3
%
% Hs : [cell] Hamiltonian. Each cell element Hs{n} describes the two-site
%       interaction between site n and n+1. Thus, Hs(1:2:end) acts on odd
%       bonds; Hs(2:2:end) on even bonds. It should satisfy numel(M) ==
%       numel(Hs) + 1.
%       The leg convention of Hs{n} are as follows:
%
%       2      4      [legs 1 and 2 are for site n;
%       |      |       legs 3 and 4 are for site n+1]
%      [ Hs{n}  ]
%       |      |
%       1      3
%
% O : [matrix] Rank-2 tensor as a local operator acting on a site.
% Nkeep : [integer] Maximum bond dimension.
% dt : [numeric] Time step size. Each real-time evolution by step dt
%       consists of three exponential terms, exp(-1i*dt/2*Hodd) *
%       exp(-1i*dt*Heven) * exp(-1i*dt/2*Hodd).
% tmax : [numeric] Maximum time range.
%
% < Output >
% ts : [numeric] Row vector of discrete time values.
% M : [cell] The final MPS after real-time evolution.
%            In right(left)-canonical form if Nstep is even(odd) 
% Ovals : [matrix] Ovals(m,n) indicates the expectation value of local
%       operator O (input) at the site n and time ts(m).
% EE : [matrix] EE(m,n) indicates the entanglement entropy (with base 2) of
%       the MPS with respect to the bipartition cutting the bond between
%       the sites n and n+1, after applying the m-th row of time evolution
%       gates. For m's that are multiples of 3, EE(m/3,:) indicates the
%       entanglement entropy at time ts(m/3). Since the base 2 is chosen,
%       the value 1 of the entanglement entropy means one "ebit".
% dw : [matrix] Discarded weights (i.e., the sum of the squares of the
%       discarded singular values) after appying a row of time evolution
%       gates. dw(m,n) corresponds to the same bond and Trotter step
%       associated with EE(m,n).
%
% Written by S.Lee (Jun.19,2017); updated by S.Lee (Jun.22,2017)
% Updated by S.Lee (Jun.07,2019): Revised for SoSe 2019.
% Updated by S.Lee (Oct.08,2022): Revised for the course at SNU.

if print_log

    tobj = tic2;
end
% % % check the integrity of input
if numel(M) ~= (numel(Hs)+1)
    error('ERR: it should be: numel(M) == (numel(Hs)+1)');
elseif ~ismatrix(O)
    error('ERR: local operator O should be rank 2.');
end

for itN = (1:numel(Hs))
    if ~all(size(Hs{itN},1) == [size(Hs{itN},2) size(M{itN},3)])
        error(['ERR: The first and second legs of Hs{', ...
            sprintf('%i',itN),'} and the third leg of M{',sprintf('%i',itN), ...
            ' should have the same dimensions.']);
    elseif ~all(size(Hs{itN},3) == [size(Hs{itN},4) size(M{itN+1},3)])
        error(['ERR: The third and fourth legs of Hs{', ...
            sprintf('%i',itN),'} and the third leg of M{',sprintf('%i',itN+1), ...
            ' should have the same dimensions.']);
    end
end
% % % 
% save initial orthonormal state 
%initial_state = M(:); 

Nstep = ceil(tmax/dt);

% results
ts = dt*(1:Nstep);
if ~isempty(O) 
    if numel(O)==numel(M)
        Ovals = zeros(Nstep,1);
    elseif numel(O)==1
        Ovals = zeros(Nstep, numel(M));
    end
else 
    Ovals = [];
end
EE = zeros(3*Nstep,numel(M)-1);
dw = zeros(size(EE));

% show information
if print_log

    fprintf('TS-1D : Imaginary-time evolution with local measurements\n');
    fprintf(['N = ',sprintf('%i',numel(M)),', Nkeep = ',sprintf('%i',Nkeep), ...
        ', dt = ',sprintf('%.4g',dt),', tmax = ',sprintf('%g',ts(end)), ...
        ' (',sprintf('%.4g',Nstep),' steps)\n']);
end
% generate the unitray operator exp(-it*H) for each two-site pairs
expH = cell(1,numel(Hs));
for it1 = (1:numel(Hs))
    if ~isempty(Hs{it1})
        sdim = [size(Hs{it1},1) size(Hs{it1},3)];
        Htmp = permute(Hs{it1},[1 3 2 4]);
        Htmp = reshape(Htmp,(sdim(1)*sdim(2))*[1 1]);
        if mod(it1,2) == 1
            ttmp = dt/2; % half time step for odd bonds, as the time evolution steps for odd bonds will happen twice
        else
            ttmp = dt;
        end
        [VH,DH] = eig(Htmp);
        eH = VH*diag(exp((-ttmp)*diag(DH)))*VH';
        expH{it1} = reshape(eH,[sdim sdim]);
    end
end
if print_log
    disptime('Transform the MPS into right-canonical form.');
end
% since the first sweep is left-to-right, bring the input into
% right-canonical form, *without* truncation.
M = canonForm(M,0,[],[]);

if print_log
    disptime('Trotter steps: start');
end
for it1 = (1:3*Nstep)
% Here we use the 2nd order Trotter step exp(-dt/2*Hodd) * exp(-dt*Heven) *
% exp(-dt/2*Hodd). That is, for the case mod(it1,3) == 2, we act the
% unitary on even bonds. Otherwise, on odd bonds.
    expHtmp = cell(1,numel(Hs));
    if mod(it1,3) == 2 % even bonds
        expHtmp(2:2:end) = expH(2:2:end);
    else % odd bonds
        expHtmp(1:2:end) = expH(1:2:end);
    end
    
    % call local function tDMRG_1sweep which is written below in this file.
    [M,EE(it1,:),dw(it1,:)] = tDMRG_1sweep(M,expHtmp,Nkeep,mod(it1,2));
    if mod(it1,3) == 0
        if ~isempty(O)
            % evaluate expectation values
            Ovals(it1/3,:) = exp_val(M,O,mod(it1,2));
        end

    end

    M = normalize(M);
    if print_log 
        if (mod(it1,round(3*Nstep/10)) == 0) || (it1 == (3*Nstep))
            disptime(['#',sprintf('%i/%i',[it1/3,Nstep]), ...
                ' : t = ',sprintf('%g/%g',[ts(it1/3),ts(end)])]);
        end
    end
end

if print_log
    toc2(tobj,'-v');
    chkmem;
end
end

function [M,EE,dw] = tDMRG_1sweep (M,expH,Nkeep,isright)
% Apply exp(-it*Heven/odd), which is a row of two-site gates acting on
% either even or odd bonds, and then truncate bonds by using SVD. After
% applying this function, left-canonical state becomes right-canonical, and
% vice versa.
%
% < Input >
% M : [cell] Input MPS.
% expH : [cell] exp(-i*H*T) unitary operators for each bond. The length
%       should satisfy numel(expH) == numel(M)-1. And the every first (or
%       second) elements should be empty, since we act either even or odd
%       bonds at once.
% Nkeep : [numeric] Maximum bond dimension.
% isright : [logical] If true, we take left-to-right sweep. Otherwise, take
%       right-to-left sweep.
% 
% < Output >
% M : [cell] MPS after applying exp(-it*H) and truncating bonds.
% EE : [numeric vector] Entanglement entropy at each bond.
% dw : [numeric vector] Discarded weights when truncating the bond
%       dimensions.
    
N = numel(M);
EE = zeros(1,N-1);
dw = zeros(1,N-1);
Skeep = 1e-8;

if isright % left -> right
    for it = (1:N-1)
        % contract M{it} and M{it+1} with expH{it}
        T = contract(M{it},3,2,M{it+1},3,1,[1 2 4 3]);
        if ~isempty(expH{it})
            T = contract(expH{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
        end
        % SVD via svdTr
        [M{it},S,V,dw(it)] = svdTr(T,4,[1 2],Nkeep,Skeep);
        M{it} = permute(M{it},[1 3 2]);
        % update M{it+1}
        M{it+1} = contract(diag(S),2,2,V,3,1,[1 3 2]);

        % normalize the singular values, to normalize the norm of MPS
        S = S/norm(S);
        % compute entanglement entropy of base 2. Be aware of zero
        % singular values!
        Spart = -(S.^2).*log(S.^2)/log(2);
        EE(it) = sum(Spart(~isnan(Spart)));
    end
    M{end} = M{end}/norm(M{end}(:)); % to normalize the norm of MPS
else % right -> left
    for it = (N-1:-1:1)
        % contract M{it} and M{it+1} with expH{it}
        T = contract(M{it},3,2,M{it+1},3,1,[1 2 4 3]);
        if ~isempty(expH{it})
            T = contract(expH{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
        end
        % SVD via svdTr
        [U,S,M{it+1},dw(it)] = svdTr(T,4,[1 2],Nkeep,Skeep);
        M{it+1} = permute(M{it+1},[1 3 2]);
        % update M{it}
        M{it} = contract(U,3,3,diag(S),2,1,[1 3 2]);
        % normalize the singular values, to normalize the norm of MPS
        S = S/norm(S);
        % compute entanglement entropy of base 2. Be aware of zero
        % singular values!
        Spart = -(S.^2).*log(S.^2)/log(2);
        EE(it) = sum(Spart(~isnan(Spart)));
    end
end

end


function M = normalize(M)
    N = numel(M); 
    MM = 1;
    for itN=(1:N)
        T1 = contract(MM,3,3,M{itN},3,1); 
        MM = contract(conj(M{itN}),3,[1,3],T1,4,[1,4]);
    end 
    for itN=(1:N)

        M{itN} = M{itN}./sqrt(MM);
    end
end 



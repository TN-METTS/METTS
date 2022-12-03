function [M, isright] = TS_1D(M, Hs, Nkeep, dt, beta, print_log)
% < Description >
%
% [M, isright] = TS_1D(M,Hs, Nkeep, dt, beta, print_log)
%
% Imaginary time evolution using second order trotter decomposition. 
% To get the partial density matrix, evole (ibeta) time with Hamiltonian. 
% The final state is |phi> such that |phi> = e^(-beta*H)|phi_0>
% Using basis, approximate the density matrix (e^(- 2*beta*H)) using many
% number of |phi_0>s
%
% It is suitable for Hamiltonian which contains only nearest neighbor interaction 
% It is revised from tDMRG function by S.Lee(Oct.08,2022)
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
% Nkeep : [integer] Maximum bond dimension.
% dt : [numeric] Time step size. Each imaginary-time evolution by step dt
%       consists of three exponential terms, exp(-dt/2*Hodd) *
%       exp(-dt*Heven) * exp(-dt/2*Hodd).
% beta : [numeric] max imaginary time. Inverse of temperature(setting boltzman constant=1) 
%
% < Output >
% M : [cell] The final MPS after real-time evolution.
%            In right(left)-canonical form if Nstep is even(odd) 
% isright: [boolean] whether right(left)-canonical form

% Rewritten by M.Kim (Nov.30,2022)

if print_log
    tobj = tic2;
end

% % % check the integrity of input
if numel(M) ~= (numel(Hs)+1)
    error('ERR: it should be: numel(M) == (numel(Hs)+1)');
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

Nstep = ceil(beta/dt);
isright = mod(Nstep, 2); %In right(left)-canonical form if Nstep is even(odd) 
ts = dt*(1:Nstep);

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
    [M] = sweep(M,expHtmp,Nkeep,mod(it1,2));
    

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

function M = sweep (M,expH,Nkeep,isright)
% Apply exp(-it*Heven/odd), which is a row of two-site gates acting on
% either even or odd bonds, and then truncate bonds by using SVD. After
% applying this function, left-canonical state becomes right-canonical, and
% vice versa.
%
% < Input >
% M : [cell] Input MPS.
% expH : [cell] exp(-H*T) unitary operators for each bond. The length
%       should satisfy numel(expH) == numel(M)-1. And the every first (or
%       second) elements should be empty, since we act either even or odd
%       bonds at once.
% Nkeep : [numeric] Maximum bond dimension.
% isright : [logical] If true, we take left-to-right sweep. Otherwise, take
%       right-to-left sweep.
% 
% < Output >
% M : [cell] MPS after applying exp(-t*H) and truncating bonds.
%       dimensions.
    
N = numel(M);
Skeep = 1e-8;

if isright % left -> right
    for it = (1:N-1)
        % contract M{it} and M{it+1} with expH{it}
        T = contract(M{it},3,2,M{it+1},3,1,[1 2 4 3]);
        if ~isempty(expH{it})
            T = contract(expH{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
        end
        % SVD via svdTr
        [M{it},S,V,~] = svdTr(T,4,[1 2],Nkeep,Skeep);
        M{it} = permute(M{it},[1 3 2]);
        % update M{it+1}
        M{it+1} = contract(diag(S),2,2,V,3,1,[1 3 2]);


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
        [U,S,M{it+1},~] = svdTr(T,4,[1 2],Nkeep,Skeep);
        M{it+1} = permute(M{it+1},[1 3 2]);
        % update M{it}
        M{it} = contract(U,3,3,diag(S),2,1,[1 3 2]);

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
    M{N} = M{N}./sqrt(MM);
    
end 



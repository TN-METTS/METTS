function equality = Test_initial_guess()

% This function tests zipup algorithm "zipup_algo"
% [B, C] = zipup_algo(W, A, Nkeep)
% generate random MPS, MPO and compare whether the function works well. 
% This test might fail but it's okay! Since the function approximates the 
% multiplication and provides initial guess of the multiplication result. 
% Returns 
%   Equality(boolean) : whether the output of the function equals expected output 
% 
% Written by M.Kim (Nov.30,2022)

% generate sample input MPS A, MPO W 
N= 4;
A = cell(1, N);
% leg order : left, right, physical
A{1} = rand(1, 2, 3);
A{2} = rand(2, 2, 3);
A{3} = rand(2, 2, 3);
A{4} = rand(2, 1, 3);

W = cell(1, N); 
% leg order : bottom, up, left, right
W{1} = rand(3, 3, 1, 5);
W{2} = rand(3, 3, 5, 5);
W{3} = rand(3, 3, 5, 5);
W{4} = rand(3, 3, 5, 1);

Nkeep = 10; 
Skeep = 1e-8;
% B ~ WA Get the result from the function  
[B, C] = zipup_algo(W, A, Nkeep, Skeep);

% generate expected output 
B_expected = 1;
for itN = (1:4)
    tmp = contract(A{itN}, 3, 3, W{itN}, 4, 2 ); % AL, AR, WD, WL, WR
    tmp = contract(tmp, 5, [1, 4],  getIdentity(A{itN},1, W{itN}, 3),3, [1, 2], [4 2 1 3] ); % AWL, WD, AR, WR, 
    tmp = contract(tmp, 4, [3, 4], getIdentity(A{itN}, 2, W{itN}, 4), 3, [1, 2]); % AWL, WD, AWR
    B_expected = contract(B_expected, 1+itN, 1+itN, tmp, 3, 1); % BL, WD, AWR
end 
B_expected= permute(B_expected, [2 3 4 5 1]);
% extract norm
norm_ = norm(B_expected(:));
% Normalize
B_expected = B_expected./norm_;

% transform B to tensor 
B_ = permute(B{1}, [1 3 2]);
for itN=(2:N)
    B_ = contract(B_, 1+itN,1+itN, permute(B{itN}, [1 3 2]), 3, 1);
end 
B_ = permute(B_, [2 3 4 5 1]);


equality = false; 
tol = 2^(-13);
if size(B_) == size(B_expected)
    if max(abs(B_ - B_expected), [],'all')< tol % compare the normalized tensor
        if abs(C-norm_)<tol % compare the norm
            equality = true;
        else 
            disptime(['Norm difference: ', sprintf('%f\n',abs(C-norm_))])
        end  

    else    
        disptime(['Max difference : ', sprintf('%f\n', max(abs(B_ - B_expected), [],'all'))])
    end 
end 

end


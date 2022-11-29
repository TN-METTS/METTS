function equality = Test_sweeps()
%TEST_SWEEPS 이 함수의 요약 설명 위치
%   자세한 설명 위치

% generate sample input MPS A, MPO W 
A = cell(1, 4);
% leg order : left, right, physical
A{1} = rand(1, 2, 3);
A{2} = rand(2, 2, 3);
A{3} = rand(2, 2, 3);
A{4} = rand(2, 1, 3);

W = cell(1, 4); 
% leg order : bottom, up, left, right
W{1} = rand(3, 3, 1, 5);
W{2} = rand(3, 3, 5, 5);
W{3} = rand(3, 3, 5, 5);
W{4} = rand(3, 3, 5, 1);

% B = WA Get the result from the function  
B = MPO_multiplication_sweep(W, A, Nkeep, Nsweep);

% generate expected output 
B_expected = 1;
for itN = (1:4)
    tmp = contract(A{itN}, 3, 3, W{itN}, 4, 1 ); % AL, AR, WU WL, WR
    tmp = contract(tmp, 5, [1, 4],  getIdentity(A{itN},1, W{itN}, 3),3, [1, 2], [4 2 1 3] ); % AWL, WU, AR, WR, 
    tmp = contract(tmp, 4, [3, 4], getIdentity(A{itN}, 2, W{itN}, 4), 3, [1, 2]); % AWL, WU, AWR
    B_expected = contract(B_expected, 1+itN, 1+itN, tmp, 3, 1); % BL, WU, AWR
end 
B_expected= permute(B_expected, [2 3 4 5 1]);

% transform B to tensor 
B_ = permute(B{1}, [1 3 2]);
for itN=(2:4)
    B_ = contract(B_, 1+itN,1+itN, permute(B{itN}, [1 3 2]), 3, 1);
end 
B_ = permute(B_, [2 3 4 5 1]);

equality = false; 
tol = 2^(-13);
if size(B_) == size(B_expected)
    if max(abs(B_ - B_expected))< tol
        equality = true;
    end 
end 

end


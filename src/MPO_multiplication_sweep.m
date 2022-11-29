function B = MPO_multiplication_sweep(W, A, Nkeep, Nsweep)

%UNTITLED 이 함수의 요약 설명 위치
%   자세한 설명 위치

%
% To guarantee proper truncation in this naive approach, one would first sweep left to right,
% performing SVDs to make the basis to the left orthogonal, but leaving the bond dimension as
% mk since the basis to the right is not orthogonal. Once at the right edge, the SVDs in the reverse
% sweep could truncate the bond dimension to m or to some specified truncation error. However,
% this procedure scales as Nm3k3
% and would be highly inefficient if k ∼ 10–100.

% Form a random
% MPS |ψBi of bond dimension m and orthogonalize it to have any arbitrary orthogonality center.
% The optimal two-site wavefunction ψ(βi−1sisi+1βi+1) at this center is then found by forming
% tensors L and R representing the product in the basis for the left and right halves of the
% system and combining them with the local W and A matrices as in figure 4. This ψB minimizes
% k|ψBi − W|ψAik2
% . One may then split ψB into matrices B
% si and B
% si+1 and sweep back and forth
% through the system, repeating this process until it converges

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end


% File: gop.m (generalized operator)
% Date: 2005/07/14
% Author: Xinhua Peng
%
% gop(s,U)
% nqubits=7
% In: single-qubit unitary operator U, qubit name s 
% Out: n-spin unitary operator which acts on qubit s with U and
% trivially on the remaining qubits
% Order: 1,2,3,4 from high to low bit
% global nqubits;
function gop=gop(s,U)
global nqubit;

for position=1:(s-1)%if for j=1:0,then ª·÷±Ω”pass
    U=kron(eye(2),U);
end

for position=s+1:nqubit
    U=kron(U,eye(2));
end

gop=U;

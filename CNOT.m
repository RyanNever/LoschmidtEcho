function out=CNOT(Nqubits,cqubit,tqubit)
% control qubit index: cqubit / target qubit index: tqubit
a=[1 0;0 0]; b=[0 0;0 1]; E=eye(2); X=[0 1;1 0];
Part1=eye(2^(min(cqubit,tqubit)-1));
if cqubit<tqubit
    Part2=mykron(a,eye(2^(abs(cqubit-tqubit)-1)),E)+mykron(b,eye(2^(abs(cqubit-tqubit)-1)),X);
else
    Part2=mykron(E,eye(2^(abs(cqubit-tqubit)-1)),a)+mykron(X,eye(2^(abs(cqubit-tqubit)-1)),b);
end
Part3=eye(2^(Nqubits-max(cqubit,tqubit)));
out=mykron(Part1,Part2,Part3);
end

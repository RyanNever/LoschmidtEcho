function QFI = mixedQFI(Uc)
global nqubit probe Ge
D=diag(probe);
TrGe=sparse(Uc'*Ge*Uc);
mD=D.*ones(2^nqubit);
Pos=find((mD.'+mD)>0.00001);
A1=(2.*(mD.'-mD).^2./(mD.'+mD)).*(abs(TrGe).^2);
QFI=sum(A1(Pos));

% QFI=sum(sum((2.*(mD.'-mD).^2./(mD.'+mD)).*(abs(TrGe).^2)));
end


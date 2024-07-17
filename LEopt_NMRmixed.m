clear all
global nqubit probe Ge
gamP=1.61;
gamH=3.977;
% gamP=1;
% gamH=1;
nqubit=10;
I=[1,0;0,1]; Ix=[0,1;1,0]/2; Iy=[0,-1i;1i,0]/2; Iz=[1,0;0,-1]/2;
eps=1e-5;
Ge=0;
Jy=0;
Jz=0;
for i=2:nqubit
    Jy=Jy+gop(i,Iy);
    Jz=Jz+gop(i,Iz);
end
rho0=(gop(1,I)+eps.*(gamP.*gop(1,Iz)+gamH.*Jz))/(2^nqubit);
for i=1:nqubit
    Ge=Ge+gop(i,Iz);
end
probe=rho0;
Ub1=expm(-1.i.*(pi/2).*(gop(1,Iy)+Jy));
Ub2=1;
for i=2:nqubit
    Ub2=CNOT(nqubit,1,i)*Ub2;
end
Ub3=expm(-1.i.*(pi/2).*gop(1,Iy));
U_ben{1}=sparse(Ub1);
U_ben{2}=sparse(Ub1*Ub2);
U_ben{3}=sparse(Ub2*Ub3);
U_ben{4}=sparse(Ub2*Ub3*Ub2);
Fq0=mixedQFI(U_ben{1});
Numop=4;
rhot=U_ben{Numop}*rho0*U_ben{Numop}';
obsv=rhot;
% obsv=U_ben{Numop}*Ge*U_ben{Numop}';
for i=1:99
    phi(i)=0.01*i*pi;
    Uen=expm(-1.i*phi(i).*Ge);
    sqrM=trace(Uen*rhot*Uen'*obsv*obsv);
    M(i)=trace(Uen*rhot*Uen'*obsv);
    DetM2(i)=real(sqrM-M(i)^2);
    Devia(i)=real(-1.i.*trace(Ge*Uen*rhot*Uen'*obsv)+1.i.*trace(Uen*rhot*Ge*Uen'*obsv));
    Devia2(i)=(real(-1.i.*trace(Ge*Uen*rhot*Uen'*obsv)+1.i.*trace(Uen*rhot*Ge*Uen'*obsv)))^2;
    Fq1=mixedQFI(U_ben{Numop});
    bound_theory(i)=1/Fq1;
    bound_cal(i)=DetM2(i)/Devia2(i);
end
xas=[0.01:0.01:0.99];

subplot(3,2,[2 4 6])
plot(xas,bound_cal/(1/Fq0),'*--')
hold on
plot(xas,bound_theory/(1/Fq0))
hold off
legend('CFI bound','QFI bound')
xlabel('\phi(\pi)')
ylim([0 1.8])
% ylim([bound_theory(1)/(1/Fq0)*0.95,bound_theory(1)/(1/Fq0)*3])
title(['\epsilon=',num2str(eps)])
subplot(3,2,1)
plot(xas,DetM2)
xlabel('\phi(\pi)')
ylabel('\Delta^2M')

subplot(3,2,3)
plot(xas,real(M))
xlabel('\phi(\pi)')
ylabel('\langle M\rangle')
subplot(3,2,5)
plot(xas,Devia)
xlabel('\phi(\pi)')
ylabel('d\langle M\rangle/d\phi')


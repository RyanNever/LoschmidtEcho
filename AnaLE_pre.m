clear all
global nqubit probe Ge
nqubit=10;
I=[1,0;0,1]; Ix=[0,1;1,0]/2; Iy=[0,-1i;1i,0]/2; Iz=[1,0;0,-1]/2;
eps=1e-5;
rho0_single=0.5*(1+eps).*[1 0;0 0]+0.5*(1-eps).*[0 0;0 1];
rho0_total=1;
Ge=0;
Jy=0;
for i=1:nqubit
    Jy=Jy+gop(i,Iy);
    Ge=Ge+gop(i,Iz);
    rho0_total=kron(rho0_total,rho0_single);
end
probe=rho0_total;
Ub1=expm(-1.i.*(pi/2).*(Jy));
Ub2=1;
for i=2:nqubit
    Ub2=CNOT(nqubit,1,i)*Ub2;
end
Ub3=expm(-1.i.*(pi/2).*gop(1,Iy));
U_ben{1}=sparse(Ub1);
U_ben{2}=sparse(Ub1*Ub2);
U_ben{3}=sparse(Ub2*Ub3);
U_ben{4}=sparse(Ub2*Ub3*Ub2);
Numop=4;
for p=0:nqubit
    lmd_s(p+1)=(0.5+0.5*eps)^(nqubit-p)*(0.5-0.5*eps)^p;
end
Fq1=mixedQFI(U_ben{Numop});
for i=1:999
    phi(i)=0.001*i*pi;
    bound_theory(i)=1/Fq1;
    LE_a(i)=0;
    sqrM_a(i)=0;
    Devia_a(i)=0;
    for j=0:fix(nqubit/2)
        if j==fix(nqubit/2) && mod(nqubit,2)==0
            LE_a(i)=LE_a(i)+0.5*nchoosek(nqubit,j)*(0.5*(lmd_s(j+1)+lmd_s(nqubit-j+1))^2+0.5*cos((nqubit-2*j)*phi(i))*(lmd_s(j+1)-lmd_s(nqubit-j+1))^2);
            sqrM_a(i)=sqrM_a(i)+0.5*nchoosek(nqubit,j)*(0.5*(lmd_s(j+1)^2+lmd_s(nqubit-j+1)^2)*(lmd_s(j+1)+lmd_s(nqubit-j+1))+0.5*cos((nqubit-2*j)*phi(i))*(lmd_s(j+1)^2-lmd_s(nqubit-j+1)^2)*(lmd_s(j+1)-lmd_s(nqubit-j+1)));
            Devia_a(i)=Devia_a(i)+0.5*nchoosek(nqubit,j)*(-0.5*(nqubit-2*j)*sin((nqubit-2*j)*phi(i))*(lmd_s(j+1)-lmd_s(nqubit-j+1))^2);
        else
            LE_a(i)=LE_a(i)+nchoosek(nqubit,j)*(0.5*(lmd_s(j+1)+lmd_s(nqubit-j+1))^2+0.5*cos((nqubit-2*j)*phi(i))*(lmd_s(j+1)-lmd_s(nqubit-j+1))^2);
            sqrM_a(i)=sqrM_a(i)+nchoosek(nqubit,j)*(0.5*(lmd_s(j+1)^2+lmd_s(nqubit-j+1)^2)*(lmd_s(j+1)+lmd_s(nqubit-j+1))+0.5*cos((nqubit-2*j)*phi(i))*(lmd_s(j+1)^2-lmd_s(nqubit-j+1)^2)*(lmd_s(j+1)-lmd_s(nqubit-j+1)));
            Devia_a(i)=Devia_a(i)+nchoosek(nqubit,j)*(-0.5*(nqubit-2*j)*sin((nqubit-2*j)*phi(i))*(lmd_s(j+1)-lmd_s(nqubit-j+1))^2);
        end
    end
    bound_cal(i)=(sqrM_a(i)-LE_a(i)^2)/Devia_a(i)^2;
end
xas=[0.001:0.001:0.999];
subplot(3,2,[2 4 6])
plot(xas,bound_cal,'*--')
hold on
plot(xas,bound_theory)
hold off
legend('CFI bound','QFI bound')
xlabel('\phi(\pi)')
ylim([bound_theory(1)*0.95,bound_theory(1)*3])
title(['\epsilon=',num2str(eps)])
subplot(3,2,1)
plot(xas,sqrM_a-LE_a.^2)
xlabel('\phi(\pi)')
ylabel('\Delta^2M')

subplot(3,2,3)
plot(xas,LE_a)
xlabel('\phi(\pi)')
ylabel('\langle M\rangle')
subplot(3,2,5)
plot(xas,Devia_a)
xlabel('\phi(\pi)')
ylabel('d\langle M\rangle/d\phi')

clear all
global nqubit probe Ge
nqubit=8;
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
rhot=U_ben{Numop}*rho0_total*U_ben{Numop}';
for p=0:nqubit
    lmd_s(p+1)=(0.5+0.5*eps)^(nqubit-p)*(0.5-0.5*eps)^p;
end
for i=1:99
    phi(i)=0.01*i*pi;
    Uen=expm(-1.i*phi(i).*Ge);
%     obsv=U_ben{Numop}*Ge*U_ben{Numop}';
    obsv=rhot;
    sqrM(i)=trace(Uen*rhot*Uen'*obsv*obsv);
    M(i)=trace(Uen*rhot*Uen'*obsv);
    DetM2(i)=real(sqrM(i)-M(i)^2);
    Devia(i)=real(-1.i.*trace(Ge*Uen*rhot*Uen'*obsv)+1.i.*trace(Uen*rhot*Ge*Uen'*obsv));
    Devia2(i)=(real(-1.i.*trace(Ge*Uen*rhot*Uen'*obsv)+1.i.*trace(Uen*rhot*Ge*Uen'*obsv)))^2;
    Fq1=mixedQFI(U_ben{Numop});
    bound_theory(i)=1/Fq1;
    bound_cal(i)=DetM2(i)/Devia2(i);
    lmd=sort(diag(probe),'descend');
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
    
end
figure(1)
plot(phi,LE_a)
hold on
plot(phi,real(M),'*')
figure(2)
plot(phi,sqrM_a)
hold on
plot(phi,real(sqrM),'*')
figure(3)
plot(phi,Devia_a)
hold on
plot(phi,Devia,'*')
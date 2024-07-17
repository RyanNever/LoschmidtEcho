clear all
warning off
global nqubit probe Ge

I=[1,0;0,1]; Ix=[0,1;1,0]/2; Iy=[0,-1i;1i,0]/2; Iz=[1,0;0,-1]/2;
eps=1e-5;
rho0_single=0.5*(1+eps).*[1 0;0 0]+0.5*(1-eps).*[0 0;0 1];
lam0=0.5*(1+eps);
lam1=0.5*(1-eps);
for ss=1:30
    nqubit=ss;
    Fq1(ss)=0;
    for j=0:nqubit-1
        Fq1(ss)=Fq1(ss)+nchoosek(nqubit-1,j)*(nqubit-2*j)^2*(lam0^(nqubit-j)*lam1^j-lam0^j*lam1^(nqubit-j))^2/(lam0^(nqubit-j)*lam1^j+lam0^j*lam1^(nqubit-j));
    end
    for p=0:nqubit
        lmd_s(p+1)=(0.5+0.5*eps)^(nqubit-p)*(0.5-0.5*eps)^p;
    end
    for i=1:999
        phi(i)=0.001*i*pi;
        bound_theory(i)=1/Fq1(ss);
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
    opt_CFI(ss)=max(max(1./bound_cal));
end
plot([1:ss],Fq1)
hold on
plot([1:ss],opt_CFI)
legend('Optimal bound','Echo-based readout protocol')
xlim([1,100])
xlabel('N')
ylabel('1/\Delta\alpha')

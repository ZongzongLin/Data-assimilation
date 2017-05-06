clc;clear;
dbstop if error
Iteration=1000;
particle=1000;
n=8;
m=n*(n-1);
anum=n;
bnum=n;
obsnum=1;
pnum=m+n+anum+bnum;
x0=rand(n,1);
SG='Gaussian';
GA='Gaussian';
sigma=0.1;
gamma=0.1;
Sigma=sigma^2*eye(pnum);
Gamma=gamma^2*eye(obsnum);
H=[eye(obsnum),zeros(obsnum,pnum-obsnum)];
adjm=gegraph(n);
beta=rand(n,1);
alpha=rand(n,1);
matt=adjm;
mattv=mtov(matt,n);
for num=1:particle
  matev(:,num)=mattv;
  alphae(:,num)=rand(anum,1);
  betae(:,num)=rand(bnum,1);
end
X=zeros(pnum,Iteration);
Y=zeros(obsnum,Iteration);
X(:,1)=[x0;mattv;alpha;beta];
Y(:,1)=H*X(:,1);
HV=zeros(pnum,particle,Iteration);
V=zeros(pnum,particle,Iteration);
V(:,:,1)=[rand(n,particle);matev;alphae;betae];
HV(:,:,1)=V(:,:,1);
for j=1:Iteration-1
    x=virusdynamic(X(:,j),mattv,alpha,beta,n,SG,sigma);
    y=virusobserve(x,H,obsnum,GA,gamma);
    X(:,j+1)=cutoff(x);
    Y(:,j+1)=cutoff(y);
    Sig=inv(inv(Sigma^(2))+H'*inv(Gamma^(2))*H);
    d=zeros(pnum,particle);
    what=zeros(1,num);
    for num=1:particle
        HV(:,num,j+1)=virusdynamic_sirsop(V(1:n,num,j),V(n+1:n+m,j),V(n+m+1:m+2*n,num,j),V(m+2*n+1:m+3*n,num,j),n,pnum,Sigma,Gamma,Sig,H,y);
        d(:,num)=Y(:,j+1)-H*HV(:,num,j+1);
        what(num)=exp(-0.5*d(:,num)'*(Gamma^(-1)+H*Sigma^(-1)*H')*d(:,num));
    end
    w=what/sum(what);
    ws=cumsum(w);
    for num=1:particle
        ix=find(ws>rand,1,'first');
        V(:,num,j+1)=HV(:,ix,j+1);
    end
end
time=[1:Iteration];
figure(1)
subplot(2,2,1);
plot(aError);
title('alpha error')
subplot(2,2,2)
plot(bError)
title('beta error')
subplot(2,2,3)
plot(xError)
title('state error')
subplot(2,2,4)
plot(time,X(1,:),time,v(1,:))
title('state and estimation')
legend('node one','estimation')
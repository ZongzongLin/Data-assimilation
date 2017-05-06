clc;clear;
Iteration=1000;
particle=1000;
n=8;
m=n*(n-1);
anum=1;
bnum=1;
obsnum=1;
pnum=m+n+anum+bnum;
x0=rand(n,1);
GA='Gaussian';
gamma=0.1;
Gamma=gamma^2*eye(obsnum);
sigma=0.01;
H=[eye(obsnum),zeros(obsnum,m+n+anum+bnum-obsnum)];
adjm=gegraph(n);
alpha=rand(anum,1);
beta=rand(bnum,1);
matt=adjm;
mattv=mtov(matt,n);
for num=1:particle
  matev(:,num)=mtov(matt,n);
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
place=ones(Iteration,1);
for j=1:Iteration-1
    x=virusdynamic(X(1:n,j),X(n+1:n+m,j),X(n+m+1:n+m+anum,j),X(n+m+anum+1:n+m+anum+bnum,j),n);
    y=virusobserve(x,H,obsnum,GA,gamma);
    X(:,j+1)=cutoff(x);
    Y(:,j+1)=cutoff(y);
    d=zeros(obsnum,particle);
    what=zeros(particle,1);
    for num=1:particle
        V(n+m+1:pnum,num,j)=cutoff(V(n+m+1:pnum,num,j)+sigma*randn(anum+bnum,1));
        HV(:,num,j+1)=cutoff(virusdynamic(V(1:n,num,j),V(n+1:n+m,num,j),V(n+m+1:n+m+anum,num,j),V(n+m+anum+1:n+m+anum+bnum,num,j),n));
        d(:,num)=Y(:,j+1)-cutoff(virusobserve(HV(:,num,j+1),H,obsnum,0,0));
        what(num)=exp(-0.5*d(:,num)'*Gamma^(-1)*d(:,num));
    end
    w=what/sum(what);
    [value,place(j+1)]=max(w);
    ws=cumsum(w);
    for num=1:particle
        ix=find(ws>rand,1,'first');
        V(:,num,j+1)=HV(:,ix,j+1);
    end
end
for j=1:Iteration
    v(:,j)=mean(V(:,:,j),2);
    aError(j)=norm(X(n+m+1:n+m+anum,1)-v(n+m+1:n+m+anum,j));
    bError(j)=norm(X(n+m+anum+1:n+m+anum+bnum,1)-v(n+m+anum+1:n+m+anum+bnum,j));
    xError(j)=norm(X(1:n,j)-v(1:n,j));
    amaxError(j)=norm(X(n+m+1:n+m+anum,1)-V(n+m+1:n+m+anum,place(j),j));
    bmaxError(j)=norm(X(n+m+anum+1:n+m+anum+bnum,1)-V(n+m+anum+1:n+m+anum+bnum,place(j),j));
    xmaxError(j)=norm(X(1:n,j)-V(1:n,j));
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
figure(2)
subplot(2,2,1);
plot(amaxError);
title('ml alpha error')
subplot(2,2,2)
plot(bmaxError)
title('ml beta error')
subplot(2,2,3)
plot(xmaxError)
title('ml state error')
subplot(2,2,4)
plot(time,X(2,:),time,v(2,:))
title('state and estimation')
legend('node two','estimation')
%testdynamic(x0,mattv,alpha,v(m+n+1:m+n+anum,Iteration),beta,v(n+m+anum+1:pnum,Iteration),n,pnum,1000)


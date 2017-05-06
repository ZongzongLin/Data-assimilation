clc;clear;
Iteration=2000;
enum=200;
n=8;
m=n*(n-1);
anum=1;
bnum=1;
obsnum=1;
pnum=m+n+anum+bnum;
x0=rand(n,1);
GA='Gaussian';
gamma=0.1;
if strcmp(GA,'Gaussian')
    Gamma=gamma^2*eye(obsnum);
elseif strcmp(GA,'Normal')
    Gamma=1/3*sigma^2*eye(obsnum);
else
    Gamma=zeros(obsnum);
end
H=[eye(obsnum),zeros(obsnum,pnum-obsnum)];
adjm=gegraph(n);
alpha=rand(anum,1);
beta=rand(bnum,1);
matt=adjm;
mattv=mtov(matt,n);
for num=1:enum
    matev(:,num)=mtov(matt,n);
    alphae(:,num)=rand(anum,1);
    betae(:,num)=rand(bnum,1);
end
X=zeros(pnum,Iteration);
Y=zeros(obsnum,Iteration);
X(:,1)=[x0;mattv;alpha;beta];
Y(:,1)=H*X(:,1);
V=zeros(pnum,enum,Iteration);
HV=zeros(pnum,enum,Iteration);
V(:,:,1)=[rand(n,enum);matev;alphae;betae];
HV(:,:,1)=V(:,:,1);
hm=zeros(pnum,Iteration);
hm(:,1)=mean(HV(:,:,1),2);
for j=1:Iteration-1
    x=virusdynamic(X(1:n,j),X(n+1:m+n,j),X(n+m+1:n+m+anum,j),X(n+m+anum+1:n+m+anum+bnum,j),n);
    y=virusobserve(x,H,obsnum,GA,gamma);
    X(:,j+1)=cutoff(x);
    Y(:,j+1)=cutoff(y);
    for num=1:enum
        HV(:,num,j+1)=virusdynamic(V(1:n,num,j),V(n+1:n+m,num,j),V(n+m+1:n+m+anum,num,j),V(n+m+anum+1:n+m+anum+bnum,num,j),n);
    end
    HV(:,:,j+1)=cutoff(HV(:,:,j+1));
    hm(:,j+1)=mean(HV(:,:,j+1),2);
    Sum=zeros(pnum);
    for num=1:enum
        Sum=Sum+(HV(:,num,j+1)-hm(:,j+1))*(HV(:,num,j+1)-hm(:,j+1))';
    end
    HC=1/(enum-1)*Sum;
    S=H*HC*H'+Gamma;
    K=HC*H'*inv(S);
    for num=1:enum
        if strcmp(GA,'Gaussian')
            V(:,num,j+1)=(eye(pnum)-K*H)*HV(:,num,j+1)+K*(cutoff(Y(:,j+1)+gamma*randn(obsnum,1)));
        elseif strcmp(GA,'Normal')
            V(:,num,j+1)=(eye(pnum)-K*H)*HV(:,num,j+1)+K*(Y(:,j+1)+gamma*(2*rand(obsnum,1)-1));
        else
            V(:,num,j+1)=(eye(pnum)-K*H)*HV(:,num,j+1)+K*Y(:,j+1);
        end
    end
    V(:,:,j+1)=cutoff(V(:,:,j+1));
end
for j=1:Iteration
    v(:,j)=mean(V(:,:,j),2);
    aError(j)=norm(X(n+m+1:n+m+anum,1)-v(n+m+1:n+m+anum,j));
    bError(j)=norm(X(n+m+anum+1:n+m+anum+bnum,1)-v(n+m+anum+1:n+m+anum+bnum,j));
    xError(j)=norm(X(1:n,j)-v(1:n,j));
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
%testdynamic(x0,mattv,alpha,v(m+n+1:m+n+anum,Iteration),beta,v(n+m+anum+1:pnum,Iteration),n,pnum,1000)





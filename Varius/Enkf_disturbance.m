clc;clear;
Iteration=1000;
enum=200;
n=8;
m=n*(n-1);
pnum=m+3*n;
obsnum=8;
x0=[1,1,0,0,0,1,1,1]';
emat=repmat(x0,1,enum);
SG='';
GA='Gaussian';
sigma=0.01;
gamma=0.01;
if strcmp(SG,'Gaussian')
    Sigma=sigma^2*eye(pnum);
elseif strcmp(SG,'Normal')
    Sigma=1/3*sigma^2*eye(pnum);%[-0.01,0.01]
end
if strcmp(GA,'Gaussian')
    Gamma=gamma^2*eye(obsnum);
elseif strcmp(GA,'Normal')
    Gamma=1/3*sigma^2*eye(obsnum);%[-0.01,0.01];
else
    Gamma=zeros(obsnum);
end
%H=obsmatrix(pnum,obsnum);
H=[eye(obsnum),zeros(obsnum,m+3*n-obsnum)];
adjacent_matrix=ones(n)-eye(n);
beta=rand(n,1);
alpha=rand(n,1);
matt=zeros(n);
for i=1:n
    for j=1:n
        if j~=i
            matt(i,j)=rand;
        end
    end
end
mattv=mtov(matt,n);
mate=zeros(n);
matev=zeros(m,enum);
for num=1:enum
    for i=1:n
        for j=1:n
            if j~=i
                mate(i,j)=rand;
            end
        end
    end
    matev(:,num)=mtov(mate,n);
    alphae(:,num)=rand(n,1);
    betae(:,num)=rand(n,1);
end
X=zeros(pnum,Iteration);
Y=zeros(obsnum,Iteration);
X(:,1)=[x0;mattv;alpha;beta];
Y(:,1)=[H*X(:,1)];
V=zeros(pnum,enum,Iteration);
HV=zeros(pnum,enum,Iteration);
V(:,:,1)=[rand(n,enum);matev;alphae;betae];
HV(:,:,1)=V(:,:,1);
hm=zeros(pnum,Iteration);
hm(:,1)=mean(HV(:,:,1),2);
for j=1:Iteration-1
    x=virusdynamic(X(1:n,j),matt,alpha,beta,n,SG,sigma);
    x=cutoff([x;mattv;alpha;beta]);
    y=cutoff(virusobserve(x,H,obsnum,GA,gamma));
    X(:,j+1)=x;
    Y(:,j+1)=y;
    for num=1:enum
        HV(1:n,num,j+1)=virusdynamic(V(1:n,num,j),V(n+1:n+m,num,j),V(n+m+1:2*n+m,num,j),V(2*n+m+1:3*n+m,num,j),n,pnum,SG,sigma);
        HV(:,num,j+1)=[HV(1:n,num,j+1);V(n+1:n+m,num,j);V(n+m+1:2*n+m,num,j);V(2*n+m+1:3*n+m,num,j)];
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
            V(:,num,j+1)=(eye(pnum)-K*H)*HV(:,num,j+1)+K*(Y(:,j+1)+gamma*rand(obsnum,1));
        else
            V(:,num,j+1)=(eye(pnum)-K*H)*HV(:,num,j+1)+K*Y(:,j+1);
        end
    end
    V(:,:,j+1)=cutoff(V(:,:,j+1));
    for j=1:Iteration-1
            Error(j)=norm(X(n+1:end,1)-mean(V(n+1:end,:,j),2));
    end
end
plot(Error)
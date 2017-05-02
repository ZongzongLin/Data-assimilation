clc;clear;
dbstop if error
Iteration=1000;
particle=1000;
n=8;
m=n*(n-1);
pnum=m+3*n;
obsnum=1;
x0=[1,1,0,0,0,1,1,1]';
SG='Gaussian';
GA='Gaussian';
sigma=0.1;
gamma=0.1;
Sigma=sigma^2*eye(pnum);
Gamma=gamma^2*eye(obsnum);
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
matev=mtov(mate,n);
for num=1:particle
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
        d(:,num)=Y(:,j+1)-virusobserve(HV(:,num,j+1),H,obsnum,'',gamma);
        what(num)=exp(-0.5*d(:,num)'*(Gamma^(-1)+H*Sigma^(-1)*H')*d(:,num));
    end
    w=what/sum(what);
    ws=cumsum(w);
    for num=1:particle
        ix=find(ws>rand,1,'first');
        V(:,num,j+1)=HV(:,ix,j+1);
    end
end
for j=1:Iteration
    Error(j)=norm(X(n+1:end,j)-mean(V(n+1:end,:,j),2));
end
plot(Error)
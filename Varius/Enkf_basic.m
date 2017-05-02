clc;clear;
dbstop if error
Iteration=5;
enum=10000;
n=8;
m=n*(n-1);
pnum=m+3*n;
x0=[1,1,0,0,0,1,1,1]';
emat=repmat(x0,1,enum);
SG='';
GA='';
sigma=0;
gamma=0;
%H=[1,zeros(1,node*node+node-1)];
Hobs=eye(n);
H=[Hobs,zeros(n,m+2*n)];
obsnum=n;
%adjm=[0,1,1,1,0,0,0,0;1,0,1,0,1,0,0,0;1,1,0,1,1,1,1,0;
%    1,0,1,0,0,0,1,0;0,1,1,0,0,1,0,0;0,0,1,0,1,0,1,1;0,0,1,1,0,1,0,0;
%    0,0,0,0,0,1,0,0;];
adjacent_matrix=ones(n)-eye(n);
beta=rand(n,1);
alpha=rand(n,1);
%gamÒì¹¹µÄÍøÂç
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
HV=zeros(pnum,enum,Iteration);
V=zeros(pnum,enum,Iteration);
V(:,:,1)=[rand(n,enum);matev;alphae;betae];
hm=zeros(pnum,Iteration);
hm(:,1)=mean(HV(:,:,1),2);
for j=1:Iteration-1
    x=virusdynamic(X(:,j),mattv,alpha,beta,n,SG,sigma);
    y=Hobs*x;
    X(:,j+1)=[x;mattv;alpha;beta];
    Y(:,j+1)=y;
    for num=1:enum
        HV(1:n,num,j+1)=virusdynamic(V(1:n,num,j),V(n+1:n+m,num,j),V(n+m+1:2*n+m,num,j),V(2*n+m+1:3*n+m,num,j),n,SG,sigma);
        HV(:,num,j+1)=[HV(1:n,n,j+1);V(n+1:n+m,num,j);V(n+m+1:2*n+m,num,j);V(2*n+m+1:3*n+m,num,j)];
    end
    hm(:,j+1)=mean(HV(:,:,j+1),2);
    Sum=zeros(pnum);
    for num=1:enum
        Sum=Sum+(HV(:,num,j+1)-hm(:,j+1))*(HV(:,num,j+1)-hm(:,j+1))';
    end
    HC=1/(enum-1)*Sum;
    S=H*HC*H';
    K=HC*H'*inv(S);
    for num=1:enum
        V(:,num,j+1)=(eye(pnum)-K*H)*HV(:,num,j+1)+K*Y(:,j+1);
    end
end
for j=1:Iteration
    Error(j)=norm(X(:,j)-mean(V(:,:,j),2));
end
plot(Error)
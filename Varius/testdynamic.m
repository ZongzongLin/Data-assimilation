function testdynamic(x,lv,alpha,alpha1,beta,beta1,n,pnum,length)
x1=x;
y=zeros(pnum,length);
yc=zeros(pnum,length);
y(1:n,1)=x;
yc(1:n,1)=x1;
for i=1:length
    y(:,i)=virusdynamic(x,lv,alpha,beta,n);
    yc(:,i)=virusdynamic(x1,lv,alpha1,beta1,n);
    x=y(1:n,i);
    x1=yc(1:n,i);
end
time=[1:length];
for i=1:n
    subplot(2,4,i)
    plot(time,y(i,:),time,yc(i,:));
end



function y=cutoff(x)
s=size(x);
y=x;
for i=1:s(1)
    for j=1:s(2)
        if x(i,j)>1
            y(i,j)=1;
        elseif x(i,j)<0
            y(i,j)=0;
        end
    end
end
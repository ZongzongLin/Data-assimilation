function H=obsmatrix(all,obs)
H=zeros(obs,all);
S=randperm(all);
S=S(1:obs);
for i=1:obs
    for j=1:all
        H(i,S(i))=1;
    end
end
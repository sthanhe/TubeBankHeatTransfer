function lcm=lcm(x)
    %Normalize
    x=reshape(x,1,numel(x));


    %Prime factors of each number
    f=arrayfun(@(x) factor(x),x,'UniformOutput',false);
    

    %Occuring prime factors
    fmax=max(cellfun(@(x) max(x),f));
    p=primes(fmax);


    %Exponent matrix
    ex=zeros(numel(f),numel(p));
    for i=1:numel(f)
        ex(i,:)=arrayfun(@(p) nnz(f{i}==p),p);
    end


    %Least common multiple
    lcm=prod(p.^max(ex,[],1));
end
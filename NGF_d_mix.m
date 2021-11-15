function [a,k] = NGF_d_mix(N,s,beta,d,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If you use this code, please cite:
%   1. G. Bianconi and C. Rahmede
%  "Network geometry with flavour: from complexity to quantum geometry"
%  Physical Review E 93, 032315 (2016).
%  2. A. P. Millan, R. Ghorbanchian, N. Defenu, F. Battiston and G. Bianconi
% "Local topological moves determine global diffusion properties of hyperbolic higher-order networks"
%  Physical Review E 104, 054302 (2021).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code that generates NGF in dimension d and flavour s=-1,0,1 costructed
% with a mixture of d-dimensional ORTHOPLEXES and d-dimensional simplices

% a adjacency matrix
% k vector of degrees   of the nodes
% This code uses
% N maximal number of nodes in the NGF
% Flavour of the NGF  s=-1,0,1
% Inverse temperature: beta>0 or beta=0
% Dimension d with d>1
% Probability of hypercubes p
% (c)  Ginestra Bianconi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization with a simgle d-hypercube

     epsilon=floor(10*rand(N+2*d,1));
    
for n1=1:d,
    for s1=0:1,
        i1=(n1-1)*2+s1+1;
        for n2=1:(d),
            for s2=0:1,
                if(abs(n2-n1)>0),
                i2=(n2-1)*2+s2+1;
                a(i1,i2)=1;
                a(i2,i1)=1;
                end
            end
         end
        end
    end


for nt=1:(2^d),
    at(nt)=1;
    a_occ(nt)=1;
    a_occ3(nt)=1;
    v=de2bi(nt-1,d);
    for i=1:d,
        j=(i-1)*2+v(i)+1;
        node(nt,i)=j;
        at(nt)=at(nt)*exp(-beta*epsilon(j));        
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At each time t=in-8 we attach a new d-hypercube with probability p otherwise we add a d-simplex


it=2*d;
while it <=(N-d),
    % Choose (d-1) face to which to attach the new d-orthoplex
    
    [I,J,V]=find(at.*a_occ);
    
    norm=sum(V);
    x=rand(1)*norm;
    for nj1=1:numel(V),
        x=x-V(nj1);
        if x<0,
            nj=J(nj1);
            break;
        end
    end
    
    
    
    a_occ(nj)=a_occ(nj)+s; %update of 1+sn
    a_occ3(nj)=a_occ3(nj)+1;%update of the generalized degree
    
    xrand=rand(1);
    
    if(xrand<p)
    %Add a new orthoplex

        for n1=1:d,
            I2b((n1-1)*2+1)=node(nj,n1);
            I2b((n1-1)*2+2)=it+n1;
        end
        it=it+d;
    
        a(I2b,I2b)=a([1:(2*d)],[1:(2*d)]);
        
        for nt2=2:2^(d),
            nt=nt+1;
            at(nt)=1;
            a_occ(nt)=1;
            a_occ3(nt)=1;
            v=de2bi(nt2-1,d);
            for i=1:d,
                j=I2b((i-1)*2+v(i)+1);
                node(nt,i)=j;
                at(nt)=at(nt)*exp(-beta*epsilon(j));
            end
        end
    end
    if(xrand>=p)
        %add a simplex
        it=it+1;
        for n1=1:d,
            j=node(nj,n1);
            a(it,j)=1;
            a(j,it)=1;
        end
        for n1=1:d,
            nt=nt+1;
            at(nt)=1;
            a_occ(nt)=1;
            a_occ3(nt)=1;
            node(nt,1)=it;
            at(nt)=1;
            j=1;
            for n2=1:d,
                if(abs(n2-n1)>0)
                    j=j+1;
                    node(nt,j)=node(nj,n2);
                    at(nt)=at(nt)*exp(-beta*epsilon(node(nj,n2)));
                    a(it,node(nj,n2))=1;
                    a(node(nj,n2),it)=1;
                end
            end
        end

    end
    k=sum(a>0);
end

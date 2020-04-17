function [dpo] = stehfest(rw,t,qinj,kj,hj,etaj,lohat)
% Function that computes the deltapo term using Stehfest algorithm
    
    % setting the number N for Stehfest algorithm
    N=12;
    % initializing the deltapo term as zero
    dpo=0.0;
    % entering Stehfest loop
    for j=1:N
        % setting the value of y
        u=log(2)*j/t;
        % setting the value of x
        x=Vj(N,j)*f(u,rw,qinj,kj,hj,etaj,lohat);
        % incrementing deltapo
        dpo=dpo+x;
    end
    
    dpo=dpo.*log(2.0)./t;
end

% function that computes the Vj coefficient, required by Stehfest algorithm
function [vj]=Vj(N,j)
    % initializing vj as zero
    vj=0.0;
    % computing the value of k1
    k1=floor((j+1)/2.0);
    % computing the value of kn
    kn=min(j,N/2);
    % implementing the sum
    for k=k1:kn
        % computing the 1st factorial expression
        a1=power(k,N*0.5)*factorial(2*k);
        % computing the 2nd factorial expression
        a2=factorial(N*0.5-k);
        % computing the 3rd factorial expression
        a3=factorial(k)*factorial(k-1);
        % computing the 4th factorial expression
        a4=factorial(j-k)*factorial(2*k-j);
        % incrementing the sum
        vj=vj+double(a1/a2/a3/a4);
    end
    % fixing the sign of vj
    vj=vj*power(-1,round(j+N/2));
end

% function in Laplace domain that must be inverted to the real field
function[y]=f(u,rw,qinj,kj,hj,etaj,lohat)
    % initializing the function's value as 0
    y=0.0;
    % getting the number of layers
    nlayers=length(kj);
    % setting the unit conversion constants
    global alphap    alphat
    % computing the contribuition of each layer to the pressure drop according to equation 14 @pdf 5 from file "02B-142746"
    for j=1:nlayers
        % sum numerator
        aux=lohat(j)*kj(j)*hj(j)*rw*sqrt(u/etaj(j)/alphat)*besselk(1,sqrt(u/etaj(j)/alphat)*rw);
        % sum denominator
        aux=aux/besselk(0,sqrt(u/etaj(j)/alphat)*rw);
        % incremeting the sum
        y=y+aux;
    end
    % inverting the sum
    y=1/y;
    % multiplying the sum by the term alphap*q*Bw/u
    y=y*alphap*qinj/u;
%     y=1./u/u;
end
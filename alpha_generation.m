function [alpha, N, gama] = alpha_generation(a, Na, Ne, Nd, d_initial, dobs, LCd)

trun = 12;
dD = zeros(Nd - trun+1, Ne);

for i = 1:Ne
    dD(:, i) = d_initial(trun:end, i) - mean(d_initial(trun:end, :), 2);
end

dD = (1/sqrt(Ne-1))*dD;
Gd = LCd*dD;
[U, S, V] = svd(Gd, 'econ');

y = LCd*(dobs(trun:end) - mean(d_initial(trun:end, :), 2));

switch a
    
    case 0
        alpha = ones(Na, 1)*Na;
        N = Na;
        gama = 1;
    case 1
        % METHOD PROPOSED BY THIAGO 1
        v = diag(S(1:Na, 1:Na));
        vsq = v.^2;
        alpha = vsq*sum(1./vsq);
        N = Na;
        gama = 1;
        
    case 3
        % METHOD PROPOSED BY RAFIEE 2017
        v = diag(S);
        sum_v = sum(v);
        a = 0.99;
        
        for i = 1:length(v)
            
            sum_aux = sum(v(1:i));
            if (sum_aux > a*sum_v)
                p = i;
                break;
            end
            
        end
        
        vm = mean(v(1:p));
        vmsq = vm^2;
        f = @(x)(fsum(x,vmsq, Na));
        
        auxc1 = 1;
        auxc2 = 10e90;
        
        c1 = f(auxc1);
        c2 = f(auxc2);
        den = 2;
        
        
        if(c1*c2 > 0)
            while (c1*c2 > 0)
                auxc2 = auxc1;
                auxc1 = auxc1/den;
                c1 = f(auxc1);
            end
        end
        
        tic
        gama = bisection(f, auxc1, auxc2);
        toc
        
        alpha = compute_alpha(gama, vmsq, Na);
        N = Na;
        
        
    case 4
        % METHOD PROPOSED BY EMERICK 2019
        anam = 1.5; amax = 100000; Namax = 100;
        
        g = @(x)( V*((S'*S + x*eye(size(S'*S)))\(S'*U'*y)) );
        f = @(x)(norm( S*V'*g(x) - U'*y ));
        emi = @(x)(f(x) - Nd);
        alphaprime = bisection(emi, 1, amax);
        if(or(alphaprime > amax, alphaprime == 0.15021992))
            alphaprime = amax;
        end
        a1 = alphaprime - 1;
        
        while (a1 < alphaprime)
            gam = @(x)(fsum_emi(x, anam, Na));
            gama = bisection(gam, 0.01, 1);
            a1 = anam*gama^(1 - Na);
            if (a1 < alphaprime)
                Na = Na + 1;
            end
            if (Na > Namax)
                Na = Namax;
                break;
            end
        end
        alpha = compute_alpha(gama, a1, Na);
        N = Na;
end
end

function [aprime] = newton_emerick(Na, Namax, amax, U, S, y)

a = Na;
ha = h_emerick(a, U, S, y);

if (ha > 0)
    aprime = a;
else
    
    a = (Na + amax)/2;
    
    for i = 1:Namax
        ha = h_emerick(a, U, S, y);
        hda = hder_emerick(a, U, S, y);
        
        aux = a;
        a = a - ha/hda;
        
        if( a > amax )
            aprime = amax;
        else if ( abs(a - aux) < 10e-3 )
                aprime = a;
            end
        end
    end
end
end

function [h] = h_emerick(a, U, S, V, y)

Nd = length(y);
Id = eye(size(S));
x = V*( S*S + a*Id )*S*U*y;

A = norm( S*V'*x - U'*y )^2 ;
h = A - Nd;

end

function [h] = hder_emerick(a, U, S, y)

v = diag(S);
v = v(v>1);
Nr = length(v);
ha = zeros(Nr, 1);
for i = 1:Nr
    ha(i) = ( ( 2*a*v(i)^2 )/( (a + v(i)^2)^3 )) * ( dot(U(:, i), y) )^2;
end

h = sum(ha);

end

function p = bisection(f,a,b)
if f(a)*f(b)>0
%     disp('Wrong choice bro')
    p = 0.15021992;
else
    p = (a + b)/2;
    err = abs(f(p));
    while err > 1e-5
        if f(a)*f(p)<0
            b = p;
        else
            a = p;
        end
        p = (a + b)/2;
        err = abs(f(p));
    end
end
end

function [alp] = compute_alpha(a, a1, Na)

alp = zeros(Na, 1);
for i = 1:Na
    alp(i) = fteste(a, a1, i);
end

end

function [fsum] = fsum(a, a1, Na)

alpha = compute_alpha(a, a1, Na);

aux = sum(1./alpha);

fsum = aux - 1;

end

function [fsum_m] = fsum_emi(gama, ana, Na)

fsum_m = gama.^Na - ana.*gama + ana - 1;

end

function [teste] = fteste(a, a1, k)

if(k == 1)
    teste = a1;
else
    teste = a.*(fteste(a, a1, k-1));
    
end
end

function [Na] = descobrir_Na(a1, ana)

a = (a1 - a1.*ana)./(ana - a1.*ana); b = ana./a1;

Na = log(b)./log(a) + 1;

end

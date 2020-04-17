function [m, d_k] = KF_EnKF_kj(rskinj, kskinj, kj, rw, nlayers, t, deltat, qinj, hj, phij, etaj, ct, mio, miw, lohat, lwhat, sw, dfw, lambdat, alphal, d_obs)

[Ne, numrsj] = size(rskinj);
[~, numksj] = size(kskinj);
[~, numkj] = size(kj);
% numkj = 0;
tam = length(d_obs);

d_k = zeros(tam, Ne);
for i = 1:Ne

    Sj=(kj(i,:)./kskinj(i, :)-1).*log(rskinj(i, :)./rw);
    [~,~,d_k(:, i),~,~,~] = compute_pwf(nlayers,length(t),t(1),deltat,rw,qinj,kj(i,:),hj,phij,etaj,ct,kskinj(i, :),rskinj(i, :),Sj,mio,miw,lohat,lwhat,sw,dfw,lambdat);
    
    %  Sj=(kj./kskinj(i, :)-1).*log(rskinj(i, :)./rw);
    % [~,~,d_k(:, i),~,~,~] = compute_pwf(nlayers,length(t),t(1),deltat,rw,qinj,kj,hj,phij,etaj,ct,kskinj(i, :),rskinj(i, :),Sj,mio,miw,lohat,lwhat,sw,dfw,lambdat);

end

% figure
% loglog(t, d_k)


% numero de pontos no 1o ciclo logaritmico (pontos desprezados pq apresentam instabilidade no modelo analitico)
a=11;

Cd = (d_obs(a+1:tam)*d_obs(a+1:tam)' + cov(d_obs(a+1:tam))*eye(size(d_obs(a+1:tam)*d_obs(a+1:tam)')));

d_obsj = zeros(size(d_k));
for i = 1:Ne
    d_obsj(:,i) = d_obs + sqrt(alphal)*sqrt(cov(d_obs))*randn(tam, 1);
end

m = [rskinj'; kskinj'; kj']; 
% Y = [rskinj'; kskinj'; kj'; d_k(a+1:tam, :)];
% Y = [rskinj'; kskinj'; d_k(a+1:tam, :)];

% H = [zeros(tam, 3) eye(tam)];
% Y_mean = mean(Y, 2);
% Cy = (1/Ne)*(Y - Y_mean)*(Y-Y_mean)';

Cdd = zeros(size(d_k(a+1:tam, 1)*d_k(a+1:tam, 1)'));
for i = 1:Ne
    Cdd = Cdd + (d_k(a+1:tam, i) - mean(d_k(a+1:tam, i)))*(d_k(a+1:tam, i)-mean(d_k(a+1:tam, i)))';
end

if (mod(tam,2)==0)
    Cmd = zeros(numrsj + numksj + numkj,tam/2-a);
else
    Cmd = zeros(numrsj + numksj + numkj,tam-a);
end

for i = 1:Ne
    Cmd = Cmd + (1/Ne)*(m(:, i) - mean(m, 2))*(d_k(a+1:tam, i) - mean(d_k(a+1:tam, :), 2))';
end

K = Cmd/(Cdd + alphal*Cd);

for i = 1:Ne
    m(:,i) = m(:,i) + K*(d_obsj(a+1:tam,i) - d_k(a+1:tam, i));
end

end
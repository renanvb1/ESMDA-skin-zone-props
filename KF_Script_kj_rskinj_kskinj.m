clc
close all;
fclose all;

Renan_multilayer_VW_1;
Sj_true = Sj;
seq=dot(qj(end,:),Sj)/sum(qj(end,:));
rskin_true=rskinj;
kskin_true=kskinj;
deltat=10.0^(0.1);
pwf = pwf + 0.015*var(pwf)*randn(length(pwf),1);
Ne = 100;
Na = 4;
marker = 3;
alphal = Na:-1:1;
alphal = alphal*(sum(1./alphal));
nparameters = 2*nlayers; %QUANTOS PARAMETROS ESTA ESTIMANDO?

% rskinj=randn(nlayers,Ne);
% rskinj=abs(rskinj)./max(abs(rskinj));
% rskinj = rskinj*0.6 + rw;
rskinj=abs(randn(nlayers, Ne)*(max(rskin_true)*1.2-min(rskin_true)*0.1)) + rw;
rskinj = rskinj';
kskinj = abs(randn(nlayers, Ne)*(max(kj)*1.2-min(kj)*0.1) + min(kj)*0.05);
kskinj= kskinj';

% kj = rand(nlayers, Ne)*(max(kj)*1.2-min(kj)*0.1) + min(kj)*0.1;
% kj=kj';

dobs = pwf;
Nd = length(dobs);
d = zeros(Nd, Ne);
tr = 11;
Cd = var(dobs(tr:end))*eye(Nd-tr);
LCd = chol(Cd, 'lower');
Cdd = zeros(Nd-tr);
Cmd = zeros(nparameters, Nd-tr);
dobsj = zeros(Nd-tr, Ne);
% Nd = length(dobs(a+1:Nd));

tic
fprintf('Iteration: 1 \n')
for i = 1:Ne
    %         Sj=(kj(i,:)./kskinj(i, :)-1).*log(rskinj(i, :)./rw);
    %         [~,~,d(:, i),~,~,~] = compute_pwf(nlayers,length(t),t(1),deltat,rw,qinj,kj(i,:),hj,phij,etaj,ct,kskinj(i, :),rskinj(i, :),Sj,mio,miw,lohat,lwhat,sw,dfw,lambdat);
    Sj=(kj./kskinj(i, :)-1).*log(rskinj(i, :)./rw);
    [~,~,d(:, i),~,~,~] = compute_pwf(nlayers,length(t),t(1),deltat,rw,qinj,kj,hj,phij,etaj,ct,kskinj(i, :),rskinj(i, :),Sj,mio,miw,lohat,lwhat,sw,dfw,lambdat);
end
%     [m, d] = KF_EnKF_kj(rskinj, kskinj, kj, rw, nlayers, t, deltat, qinj, hj, phij, etaj, ct, mio, miw, lohat, lwhat, sw, dfw, lambdat, alphal, pwf);

d_initial = d;

[alphal, Na, gama] = alpha_generation(marker, Na, Ne, Nd, d_initial, dobs, LCd);

for ii = 1:Na
    
   
    if(ii>1)
        fprintf('Iteration: %d \n', ii)
        for i = 1:Ne
            %         Sj=(kj(i,:)./kskinj(i, :)-1).*log(rskinj(i, :)./rw);
            %         [~,~,d(:, i),~,~,~] = compute_pwf(nlayers,length(t),t(1),deltat,rw,qinj,kj(i,:),hj,phij,etaj,ct,kskinj(i, :),rskinj(i, :),Sj,mio,miw,lohat,lwhat,sw,dfw,lambdat);
            Sj=(kj./kskinj(i, :)-1).*log(rskinj(i, :)./rw);
            [~,~,d(:, i),~,~,~] = compute_pwf(nlayers,length(t),t(1),deltat,rw,qinj,kj,hj,phij,etaj,ct,kskinj(i, :),rskinj(i, :),Sj,mio,miw,lohat,lwhat,sw,dfw,lambdat);
        end
    end
    %     [m, d] = KF_EnKF_kj(rskinj, kskinj, kj, rw, nlayers, t, deltat, qinj, hj, phij, etaj, ct, mio, miw, lohat, lwhat, sw, dfw, lambdat, alphal, pwf);
    
    %     m = [rskinj'; kskinj'; kj'];
    m = [rskinj'; log(kskinj')];
    
    Cdd = zeros(Nd-tr);
    Cmd = zeros(nparameters, Nd-tr);
    for i = 1:Ne
        Cdd = Cdd + (1/Ne)*(d(tr+1:Nd, i) - mean(d(tr+1:Nd, :), 2))*(d(tr+1:Nd, i) - mean(d(tr+1:Nd, :), 2))';
        Cmd = Cmd + (1/Ne)*(m(:, i) - mean(m, 2))*(d(tr+1:Nd, i) - mean(d(tr+1:Nd, :), 2))';
        %        dobsj(:, i) = dobs + var(dobs)*alphal(ii)*normrnd(0, 1, [tam 1]);
    end
    
    K = Cmd/(Cdd + alphal(ii)*Cd);
    %+ sqrt(alphal(ii))*LCd*randn(Nd, 1)
    
    m = real(m + K*(dobs(tr+1:Nd) - d(tr+1:Nd, :)));
    
    rskinj = (m(1:nlayers, :)');
    kskinj = exp(m(nlayers+1:2*nlayers, :)');
    
    kskinj(kskinj>1000) = 1000;
    %     kj = (m(2*nlayers+1:3*nlayers, :)');
end

atime = toc;
% clock

Untitled5;
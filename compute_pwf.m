function [t,tp,pwf,deltapo,deltapl,qj] = compute_pwf(nlayers,dim,t0,deltat,rw,qinj,kj,hj,phij,etaj,ct,kskinj,rskinj,Sj,mio,miw,lohat,lwhat,sw,dfw,lambdat)
    % function that computes the wellbore pressure data
    
    % setting the unit conversion constant
    global alphap
    
    % initializing the time, pressure and pressure derivative vectors
    if (mod(dim,2)==0)
        t=zeros(dim/2,1);
        pwf=t;
        % dpwf=zeros(dim,1);
        deltapo=t;
        deltapl=t;
        qj=zeros(dim/2,nlayers);
        Aj=ones(dim/2,nlayers);
        % filling the time vector
        t(1)=t0;
        for ii=2:dim/2
            t(ii)=t(ii-1)*deltat;
        end
        tp=t(dim/2);
        for ii=dim/2+1:dim
            t(ii)=tp+t(ii-dim/2);
        end
        
        % computing the reservoir equivalent properties
        [h,keq,phict,flowcap,lohatm]=compute_equi_props(kj,hj,phij,ct,lohat);
        
        %initializing deltaPskin as zero
        dpskin=0.0;
        
        % initializing the waterfront profile in all layers
        rf=zeros(length(sw),nlayers);
        
        %entering the time loop during the injection period
        for ii=1:dim/2
            % before computing any pressure data, updating flow-rates in each layer
            qj(ii,:)=update_flowrate(ii,dpskin,Aj,pwf,deltapo,kj,hj,flowcap,qinj);
            % during injection, update the reservoir mechanical skin an waterfront radii
            S=dot(Sj,qj(ii,:))/qinj;
            % computing deltaPskin as described by Hawkins (file 732...)
            dpskin=alphap*qinj*S/keq/h/lohatm;
            % computing the deltapo term
            deltapo(ii)=stehfest(rw,t(ii),qinj,kj,hj,etaj,lohat);
            % computing the aj coefficients in all layers
            for j=1:nlayers
                % for each layer, computing the Aj coefficient
                [Aj(ii,j),v]=Aj_Rj(ii,j,rw,t,tp,qinj,qj,kj,hj,phij,keq,phict,kskinj,rskinj,lohat,lwhat,sw,dfw,lambdat);
                rf(:,j)=v;
                % incrementing the sum
                deltapl(ii)=deltapl(ii)+1/Aj(ii,j);
            end
            % computing the deltaplamda term
            deltapl(ii)=qinj/deltapl(ii);
            % computing the wellbottom hole pressure
            pwf(ii)=deltapo(ii)+deltapl(ii)+dpskin;
        end
        
        for ii=dim/2+1:dim
            % computing the deltapo term applying the superposition principle
            a=stehfest(rw,t(ii),qinj,kj,hj,etaj,lohat);
            b=stehfest(rw,t(ii)-tp,-qinj,kj,hj,etaj,lohat);
            deltapo(ii)=a+b;
            % computing the rj coefficients in all layers
            for j=1:nlayers
                % for each layer, computing the rj coefficient
                [Aj(ii,j),v]=Aj_Rj(ii,j,rw,t,tp,qinj,qj,kj,hj,phij,keq,phict,kskinj,rskinj,lohat,lwhat,sw,dfw,lambdat);
                rf(:,j)=v;
                %incrementing the sum
                deltapl(ii)=deltapl(ii)+1/Aj(ii,j);
            end
            % computing the deltaplamda term
            deltapl(ii)=1/deltapl(ii);
            % computing the wellbottom hole pressure
            pwf(ii)=deltapo(ii)+deltapl(ii);
        end
    else
        t=zeros(dim,1);
        pwf=t;
        % dpwf=zeros(dim,1);
        deltapo=t;
        deltapl=t;
        qj=zeros(dim,nlayers);
        Aj=ones(dim,nlayers);
        % filling the time vector
        t(1)=t0;
        for ii=2:dim
            t(ii)=t(ii-1)*deltat;
        end
        tp=t(dim);
        
        % computing the reservoir equivalent properties
        [h,keq,phict,flowcap,lohatm]=compute_equi_props(kj,hj,phij,ct,lohat);
        
        %initializing deltaPskin as zero
        dpskin=0.0;
        
        % initializing the waterfront profile in all layers
        rf=zeros(length(sw),nlayers);
        
        %entering the time loop during the injection period
        for ii=1:dim
            % before computing any pressure data, updating flow-rates in each layer
            qj(ii,:)=update_flowrate(ii,dpskin,Aj,pwf,deltapo,kj,hj,flowcap,qinj);
            % during injection, update the reservoir mechanical skin an waterfront radii
            S=dot(Sj,qj(ii,:))/qinj;
            % computing deltaPskin as described by Hawkins (file 732...)
            dpskin=alphap*qinj*S/keq/h/lohatm;
            % computing the deltapo term
            deltapo(ii)=stehfest(rw,t(ii),qinj,kj,hj,etaj,lohat);
            % computing the aj coefficients in all layers
            for j=1:nlayers
                % for each layer, computing the Aj coefficient
                [Aj(ii,j),v]=Aj_Rj(ii,j,rw,t,tp,qinj,qj,kj,hj,phij,keq,phict,kskinj,rskinj,lohat,lwhat,sw,dfw,lambdat);
                rf(:,j)=v;
                % incrementing the sum
                deltapl(ii)=deltapl(ii)+1/Aj(ii,j);
            end
            % computing the deltaplamda term
            deltapl(ii)=qinj/deltapl(ii);
            % computing the wellbottom hole pressure
            pwf(ii)=deltapo(ii)+deltapl(ii)+dpskin;
        end
    end
end

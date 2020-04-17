function [aj,rf] = Aj_Rj(ii,j,rw,t,tp,qinj,qj,kj,hj,phij,keq,phict,kskinj,rskinj,lohat,lwhat,sw,dfw,lambdat)
    % function that computes the Aj or Rj coefficient
    
    global alphap
    compsw=length(sw);
    % during injection, compute the aj coefficient=========================
    if t(ii)<=tp
        % updating the waterfront profile
        rf=computerf(ii,j,rw,t,qj,hj,phij,dfw);
        % initializing aj as zero
        aj=0.0;
        % computing the integral expression required to evaluate aj
        for u=compsw:-1:2
            % numerically integrating expression (lambdahat/lambdatotal - 1) d ln(r)/r using the trapeze rule
            aj=aj+log(rf(u-1)/rf(u))*(lohat(j)/lambdat(u,j)-1+lohat(j)/lambdat(u-1,j)-1)/2;
        end
        % if there is formation damage, applying the proper adjustment
        if rskinj(j)>rw
            % if the waterfront is within the damaged zone, use equation 44 @pdf 29
            if rf(1)<rskinj(j)
                % equation 44 @ pdf 29 from file "03 ... n-camadas.pdf"
                aj=aj*alphap/kskinj(j)/hj(j)/lohat(j);
            else
                % initializing the second integral in equation 45 @pdf 29
                intsj=0.0;
                % initializing u as the water saturation vector length
                u=compsw;
                % incrementing the auxiliary integral for all radii smaller than rskinj
                while(rf(u)<=rskinj(j))
                    % checking if the next radial step is out of the damaged zone
                    if (rf(u-1)>rskinj(j))
                        % if the next radius is out of the damaged zone, interpolate the total mobility
                        lambdatskin=lambdat(u)+(lambdat(u-1) - lambdat(u))*(rskinj(j) - rf(u)) / (rf(u-1) - rf(u));
                        % incrementing the auxiliary variable using the interpolated total mobility
                        intsj=intsj+log(rskinj(j)/rf(u))*(lohat(j)/lambdat(u)-1+lohat(j)/lambdatskin-1)/2;
                    % if the next radius is inside the damaged zone, no interpolation is required
                    else
                        % incrementing the auxiliary integral
                        intsj=intsj+log(rf(u-1)/rf(u))*(lohat(j)/lambdat(u)-1+lohat(j)/lambdat(u-1)-1)/2;
                    end
                    % decreasing the waterfront radius index until the skin radius is reached
					u = u - 1;
                end
                % multiplying the second integral by (kj/kskinj-1)
                intsj=intsj*(kj(j)/kskinj(j)-1);
                % adding the second integral to aj
                aj=aj+intsj;
                % multiplying aj by alphap/(k*h*lohat)
                aj=aj*alphap/kj(j)/hj(j)/lohat(j);
            end
        else
            % if there is no formation damage, multiplying aj by alphap/(k*h*lohat)
            aj=aj*alphap/kj(j)/hj(j)/lohat(j);
        end
%during falloff, compute the rj coefficient====================================================================================
    else
        % initializing the variables that will store the rj coefficients
        rj=0.0;
        rjskin=0.0;
        % updating the waterfront profile
        rf=computerf(length(t)/2,j,rw,t,qj,hj,phij,dfw);
        % selecting the desired total flow-rate approximation (AS APROXIMACOES SEGUEM A ORDEM DO ARQUIVO "04-FUNDAMENTOS...")
        approx=1;
        % defining the number of partitions between the wellbore and the first point inside the flooded region
        n=1000;
        % defining the radial step
        dr = (rf(compsw-1) - rf(compsw)) / n;
        % calculating the flow-rate integral between the wellbore and the first point inside the flooded region
        for u=0:n-1
            % setting the value of the current radial step
            r1 = rf(compsw) + u * dr;
            % setting the value of the next radial step
			r2 = r1 + dr;
            % computing the oil flow-rate at the current radial step
			qos1 = approxlog(ii, r1, t, tp, rw, qinj, keq, lohat(j), phict);
            % computing the oil flow-rate at the next radial step
            qos2 = approxlog(ii, r2, t, tp, rw, qinj, keq, lohat(j), phict);
            % computing the oil flow-rate at the midpoint
            qos3 = approxlog(ii, (r1+r2)/2, t, tp, rw, qinj, keq, lohat(j), phict);
            % computing total flow-rate at the current radial step
			qts1 = calcqts(ii, r1, approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
            % computing total flow-rate at the next radial step
            qts2 = calcqts(ii, r2, approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
            % computing total flow-rate at the midpoint
            qts3 = calcqts(ii, (r1+r2)/2, approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
            
            % calculating the total mobility at the intermediate points using a linear interpolation
            lt1=lambdat(compsw,j)+(lambdat(compsw-1,j)-lambdat(compsw))*(r1-rf(compsw))/(rf(compsw-1)-rf(compsw));
            lt2=lambdat(compsw,j)+(lambdat(compsw-1,j)-lambdat(compsw))*(r2-rf(compsw))/(rf(compsw-1)-rf(compsw));
            lt3=(lt1+lt2)*0.5;
            
            % numerically integrating the expression using the Simpson rule
            rj=rj+log(r2/r1)*(qts1/lt1-qos1/lohat(j)+qts2/lt2-qos2/lohat(j)+4*qts3/lt3-4*qos3/lohat(j))/6;
            
            % if there is formation damage, then the auxiliary rj must also be updated
            if (rskinj(j)>rw)
                if (r2<rskinj(j))
                    % if both points are inside the damaged zone:
                    rjskin = rjskin + log(r2 / r1)*(qts1 / lt1 - qos1 / lohat(j) + qts2 / lt2 - qos2 / lohat(j) + 4 * qts3 / lt3 - 4 * qos3 / lohat(j)) / 6;
                else
                    % if only the current radial step is inside the damaged zone:
                    if (r1<rskinj(j))
                        % estimating oil and total flow-rates at the skin radius
                        qosskin = approxlog(ii, rskinj(j), t, tp, rw, qinj, keq, lohat(j), phict);
                        qtsskin = calcqts(ii, rskinj(j), approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
                        % estimating  total mobility at the skin radius
                        ltskin = lambdat(compsw,j) + (lambdat(compsw-1,j) - lambdat(compsw))*(rskinj(j) - rf(compsw)) / (rf(compsw-1) - rf(compsw));
                        % if layer skin radius is between r1 and r2, using the trapeze rule to calculate this integral term
                        rjskin = rjskin + log(rskinj(j) / r1)*(qts1 / lt1 - qos1 / lohat(j) + qtsskin / ltskin - qosskin / lohat(j)) / 2;
                    end
                end
            end
        end
        % calculating the flow-rate integral for the other radii
        for u=compsw-1:-1:2
            % updating the radial step
            r1=rf(u);
            % updating the next radius
            r2=rf(u-1);
            % computing qos at the 1st radius
            qos1=approxlog(ii, r1, t, tp, rw, qinj, keq, lohat(j), phict);
            % computing qos an the 2nd radius
            qos2 = approxlog(ii, r2, t, tp, rw, qinj, keq, lohat(j), phict);
            % computing qos at the midpoint
            qos3=approxlog(ii, (r1+r2)*0.5, t, tp, rw, qinj, keq, lohat(j), phict);
            % computing total flow-rate at the current radial step
			qts1 = calcqts(ii, r1, approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
            % computing total flow-rate at the next radial step
            qts2 = calcqts(ii, r2, approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
            % computing total flow-rate at the midpoint
            qts3 = calcqts(ii, (r1+r2)/2, approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
            % updating total mobility at the current radial step
            lt1=lambdat(u,j);
            % updating total mobility at the next radial step
            lt2=lambdat(u-1,j);
            % interpolating total mobility at the midpoint
            lt3=(lt1+lt2)/2;%lt1+(lt2-lt1)/(r2-r1)*((r1+r2)/2-r1);%
            
            % numerically integrating the expression using the Simpson rule
            rj=rj+log(r2/r1)*(qts1/lt1-qos1/lohat(j)+qts2/lt2-qos2/lohat(j)+4*qts3/lt3-4*qos3/lohat(j))/6;
            
            % if there is formation damage, then the auxiliary rj must also be updated
            if (rskinj(j)>rw)
                if (r2<rskinj(j))
                    % if both points are inside the damaged zone:
                    rjskin = rjskin + log(r2 / r1)*(qts1 / lt1 - qos1 / lohat(j) + qts2 / lt2 - qos2 / lohat(j)) / 2;
                elseif (r1<rskinj(j))
                    % calculating flow-rates at the skin radius
                    qosskin = approxlog(ii, rskinj(j), t, tp, rw, qinj, keq, lohat(j), phict);
                    qtsskin = calcqts(ii, rskinj(j), approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
                    % calculating total mobility at the skin radius
                    ltskin = lambdat(compsw,j) + (lambdat(compsw-1,j) - lambdat(compsw))*(rskinj(j) - rf(compsw)) / (rf(compsw-1) - rf(compsw));
                    % if layer skin radius is between r1 and r2, using the trapeze rule to calculate this integral term
                    rjskin = rjskin + log(rskinj(j) / r1)*(qts1 / lt1 - qos1 / lohat(j) + qtsskin / ltskin - qosskin / lohat(j)) / 2;
                end
            end
        end
        % multiplying rjskin by (kj/kskinj-1)
        if rjskin~=0
            rjskin=rjskin*(kj(j)/kskinj(j)-1);
        end
        % multiplying rj by alphap/(kj*hj)
        aj=(rj+rjskin)*alphap/kj(j)/hj(j);
    end
end

% function that computes the waterfront profile at a given time
function [rf]=computerf(ii,j,rw,t,qj,hj,phij,dfw)
    % computing the integral at the first timestep
    integral=qj(1,j)*t(1);
    % numerically integrating (using the trapeze rule) the expression required to determine waterfront radius
    for u=1:ii-1
        integral=integral+(t(u+1)-t(u))*(qj(u+1,j)+qj(u,j))/2;
    end
    % calculating the waterfront radii using Buckley-Leverett theory for radial geometry
    rf=sqrt(integral.*dfw(:,j)./24./pi./phij(j)./hj(j)+rw.*rw);
end

% function that computes the flow-rate at a given point of the reservoir using the line-source logarithmic approximation
function [qd]=approxlog(index, ray, t, tp, rw, qinj, keq, lambdahat, phict)
    global alphat
    % % initializing the dimensionless flow-rate as zero
    % qd=0.0;
    % calculating the dimensionless radius
    rd=ray/rw;
    % % initializing the dimensionless time (and dimensionless deltaT)
    % td=1.0; deltatd=1.0;
    % computing the dimensionless times
    td=alphat*keq*lambdahat*t(index)/(phict*rw*rw);
    deltatd=alphat*keq*lambdahat*(t(index)-tp)/(phict*rw*rw);
    % calculating the dimensionless flow-rate with respect to td using the logarithmic approximation
    qd=exp(-rd*rd /(4*td));
    % during falloff, adjust the dimensionless flow-rate to account for deltatd's contribuition
    if t(index)>tp
        qd=qd-exp(-rd*rd/4/deltatd);
    end
    % dimensionalizing qd
    qd=qd*qinj;
end

% function that computes the total flow-rate at a given radius during falloff period
function [qts]=calcqts(index, ray, approx, t, tp, rw, qinj, keq, lohatm, phict, lwhatm)
    % initializing qts as zero
    qts=0.0;
    if (approx==1)
        % according to the 1st aproximation, qts = qos
        qts=approxlog(index, ray, t, tp, rw, qinj, keq, lohatm, phict);
    end
    if (approx==2)
        % according to the 2nd approximation, qts = qws
        qts=approxlog(index, ray, t, tp, rw, qinj, keq, lwhatm, phict);
    end
    if (approx==3)
        % according to the 3rd approximation, qts is an averaged mean using water and oil properties @ their endpoint saturations
        qos=approxlog(index, ray, t, tp, rw, qinj, keq, lohatm, phict);
        qws=approxlog(index, ray, t, tp, rw, qinj, keq, lwhatm, phict);
        qts=(qos*lohatm+qws*lwhatm)/(lohatm+lwhatm);
    end    
end







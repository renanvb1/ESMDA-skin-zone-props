function [dpwf]=compute_derivative(t,pwf,tp)
% Function that computes the pressure derivative
    
    % setting the vector length
    dim=length(pwf);
    % initializing the pressure derivative vector
    dpwf=zeros(size(pwf));
    
    if (mod(dim,2)==0)
        %computing the pressure derivative at the first point
        dpwf(1)=(pwf(2)-pwf(1))/log(t(2)/t(1));
        % computing the pressure derivative for the 1st half of the vector
        for ii=2:dim/2-1
            % for every point with at least 2 neighbours, compute the derivative using Bourdet
            a = (pwf(ii + 1) - pwf(ii)) / log(t(ii + 1) / t(ii))*log(t(ii) / t(ii - 1));
            b = (pwf(ii) - pwf(ii - 1)) / log(t(ii) / t(ii - 1))*log(t(ii + 1) / t(ii));
            % Bourdet derivative available @pdf 74 from file "notas cap 3"
            dpwf(ii)=(a+b)/log(t(ii + 1) / t(ii-1));
        end
        % computing the derivative at the last point of the injection period
        a=(pwf(dim/2)-pwf(dim/2-2))/log(t(dim/2)/t(dim/2-2));
        dpwf(dim/2)=a;
        
        % computing pressure derivative at the first falloff point
        teq0=tp*(t(dim/2+1)-tp)/t(dim/2+1);
        teq1=tp*(t(dim/2+2)-tp)/t(dim/2+2);
        dpwf(dim/2+1)=(pwf(dim/2+2)-pwf(dim/2+1))/log(teq0/teq1);
        % for every point with at least 2 neighbours, compute the derivative using Bourdet
        for ii=dim/2+2:dim-1
            % computing Agarwal's equivalent times
            teq0=tp*(t(ii-1)-tp)/t(ii-1);
            teq1=tp*(t(ii-0)-tp)/t(ii-0);
            teq2=tp*(t(ii+1)-tp)/t(ii+1);
            % computing the 2 terms in Bourdet derivative
            a = (pwf(ii+1) - pwf(ii))*log(teq0 / teq1) / (log(teq2 / teq1)*log(teq2 / teq0));
            b = (pwf(ii) - pwf(ii-1))*log(teq1 / teq2) / (log(teq1 / teq0)*log(teq2 / teq0));
            dpwf(ii)=a+b;
        end
        % computing the derivative at the last point
        teq1=tp*(t(dim-1)-tp)/t(dim-1);
        teq2=tp*(t(dim)-tp)/t(dim);
        dpwf(dim)=(pwf(dim)-pwf(dim-1))/log(teq1/teq2);
    else
        %computing the pressure derivative at the first point
        dpwf(1)=(pwf(2)-pwf(1))/log(t(2)/t(1));
        % computing the pressure derivative for the 1st half of the vector
        for ii=2:dim-1
            % for every point with at least 2 neighbours, compute the derivative using Bourdet
            a = (pwf(ii + 1) - pwf(ii)) / log(t(ii + 1) / t(ii))*log(t(ii) / t(ii - 1));
            b = (pwf(ii) - pwf(ii - 1)) / log(t(ii) / t(ii - 1))*log(t(ii + 1) / t(ii));
            % Bourdet derivative available @pdf 74 from file "notas cap 3"
            dpwf(ii)=(a+b)/log(t(ii + 1) / t(ii-1));
        end
        % computing the derivative at the last point of the injection period
        a=(pwf(dim)-pwf(dim-2))/log(t(dim)/t(dim-2));
        dpwf(dim)=a;
    end
end


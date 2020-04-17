function [sw,lambdat,dfw,ct]=fill_data(compsw,krw,kro,swi,sor,mio,miw,co,cw,cr)
    % function receives the relative permeability and viscosity data
    % function computes the fractional flow and total mobility data
    
    % intializing the output and auxiliary parameters
    sw=zeros(compsw,1);
    fw=zeros(compsw,1);
    dfw=zeros(compsw,1);
    lambdat=zeros(compsw,1);
    
    % computing total compressibility
    ct=cr+swi*cw+(1-sor)*co;
    
    deltasw=(1-sor-swi)/(compsw-1);
    
%     for j=1:nlayers
%         % filling the water saturation
%         sw(:,j)=swi:deltasw:1-sor;
%         % filling the total mobility data
%         lambdat(:,j)=kro./mio(j)+krw./miw;
%         % filling the fractional flow data
%         fw(:,j)=(krw(:)./miw)./lambdat(:,j);
%     end
    % filling the water saturation
    sw(:)=swi:deltasw:1-sor;
    % filling the total mobility data
    lambdat(:)=kro./mio+krw./miw;
    % filling the fractional flow data
    fw(:)=(krw(:)./miw)./lambdat(:);
    
    % numerically computing the fractional flow derivative at the 1st point
    dfw(1,:)=(fw(2,:)-fw(1,:))/(sw(2,:)-sw(1,:));
    % the last fractional flow derivative must be equal to zero
    dfw(compsw,:)=0;
    % for each point with 2 neighbours, computing the derivative using Bourdet
    for ii=2:compsw-1
        % computing the auxiliary variables required by the Bourdet derivative
        a=((fw(ii+1,:)-fw(ii,:))/(sw(ii+1,:)-sw(ii,:)))*(sw(ii,:)-sw(ii-1,:));
        b=((fw(ii,:)-fw(ii-1,:))/(sw(ii,:)-sw(ii-1,:)))*(sw(ii+1,:)-sw(ii,:));
        % computing the Bourdet derivative
        dfw(ii,:)=(a+b)/(sw(ii+1,:)-sw(ii-1,:));
    end
    
    % using the Welge method to adjust the fractional flow derivative
    
%     for j=1:nlayers
%         % starting at the end of the derivative curve
%         index=compsw;
%         % computing the auxiliary fractional flow derivative
%         dfwaux=(fw(index,j)-fw(1,j))/(sw(index,j)-sw(1,j));
%         % while the auxiliary derivative is higher than the calculated derivative, keep searching for the tangency point
%         while (dfw(index,j)<dfwaux && index >1)
%             % updating the auxiliary fractionary flow derivative
%             dfwaux=(fw(index,j)-fw(1,j))/(sw(index,j)-sw(1,j));
%             % updating the saturation index
%             index=index-1;
%         end
%         % only update the fractional flow derivative data if the tangency point is not at the first saturation
%         if (index>1)
%             % replacing all fractional flow derivatives below the tangency point
%             for ii=1:index
%                 dfw(ii,j)=dfwaux;
%             end
%         end
%     end
    % starting at the end of the derivative curve
    index=compsw;
    % computing the auxiliary fractional flow derivative
    dfwaux=(fw(index)-fw(1))/(sw(index)-sw(1));
    % while the auxiliary derivative is higher than the calculated derivative, keep searching for the tangency point
    while (dfw(index)<dfwaux && index >1)
        % updating the auxiliary fractionary flow derivative
        dfwaux=(fw(index)-fw(1))/(sw(index)-sw(1));
        % updating the saturation index
        index=index-1;
    end
    % only update the fractional flow derivative data if the tangency point is not at the first saturation
    if (index>1)
        % replacing all fractional flow derivatives below the tangency point
        for ii=1:index
            dfw(ii)=dfwaux;
        end
    end
    
end


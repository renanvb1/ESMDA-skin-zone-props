% function that updates the flow-rate at a given instant
function [qaux]=update_flowrate(ii,dpskin,Aj,pwf,deltapo,kj,hj,flowcap,qinj)
    % computing the initial flow-rates (initial flow-rate is computed in the same way both for horizontal and vertical wells)
    if (ii==1)
        qaux=qinj.*kj.*hj./flowcap;
%         n = length(hj);
%         qaux = (qinj/n)*ones(size(kj));
    % updatinf flow-rates for every other timestep
    else
        %for each layer, updating the flow-rates using the formulation proposed by Barreto et al.
        qaux=(pwf(ii-1)-deltapo(ii-1)-dpskin)./Aj(ii-1,:);
    end
end


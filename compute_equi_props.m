function [h,keq,phict,flowcap,lohatm] = compute_equi_props(kj,hj,phij,ct,lohat)
% Function that computes the reservoir equivalent properties

% computing the reservoir total thickness
h=sum(hj);
% computing the reservoir flow capacity
flowcap=dot(kj,hj);%kj*hj;
% computing the reservoir equivalent permeability
keq=flowcap/h;
% computing the averaged compressibility-thickness product
phict=sum(phij.*ct.*hj)/h;
% computing the averaged endpoint oil mobility
lohatm=dot(lohat,hj)/h;
end


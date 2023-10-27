function [niceRho,niceM] = PutStrategyInNiceBasis(rho,M)
%PutStrategyInNiceBasis Puts strategy in a nicer basis, with first measurement being in Z-basis
%   Assumes for now that everything is qubits.

    n = length(M);
    d = 2;

    UAll = cell(1,n); % Unitaries for each party taking measurement to computational basis
    for k = 1:n
        U = zeros(d,d);
        for a = 1:2
            % We play it safe and check which eigenvalue is biggest
            [V,eigvals] = eig(M{k}(:,:,a,1),'vector');
            if max(eigvals) - min(eigvals) > 1 - 1e-5 % can only do this when measurement rank-1!
                [~, ind] = sort(eigvals);
                U(a,:) = V(:,ind(end))'; % row of U is given by max eigenvector of measurement
            else
                id = eye(2);
                U(a,:) = id(a,:);
            end
        end
        UAll{k} = U;
    end
    
    niceRho = Tensor(UAll)*rho*Tensor(UAll)';
    
    niceM = cell(1,n);
    for k = 1:n
        niceM{k} = zeros(d,d,2,2);
        for t = 1:2
            for a = 1:2
                niceM{k}(:,:,a,t) = UAll{k}*M{k}(:,:,a,t)*UAll{k}';
            end
        end
    end
    
    % We also make sure 1st outcome of 2nd measurement has positive off-diagonals, to force same 
    % convention between outcomes
    sz = [1,0;0,-1];
    URho = 1;
    for k = 1:n
        if real(niceM{k}(1,2,1,2)) < 0
            % Still change the 1st measurement in case we couldn't make it a Z-basis
            niceM{k}(:,:,1,1) = sz*niceM{k}(:,:,1,1)*sz;
            niceM{k}(:,:,2,1) = sz*niceM{k}(:,:,2,1)*sz;
            niceM{k}(:,:,1,2) = sz*niceM{k}(:,:,1,2)*sz;
            niceM{k}(:,:,2,2) = sz*niceM{k}(:,:,2,2)*sz; 
            URho = Tensor(URho,sz);
        else
            URho = Tensor(URho,eye(2));
        end
    end
    niceRho = URho*niceRho*URho';
    
end


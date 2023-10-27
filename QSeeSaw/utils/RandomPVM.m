function M = RandomPVM(d,k,ranks)
%RandomPVM Generates a random PVM (projection-valued measurement)
%   M = RandomPVM(d,k,ranks) takes as input the dimension d and
%   the ranks of the k corresponding PVM elements
%   Must satisfiy sum(ranks) == d and length(ranks) == k
%   By default rank is shared as equally as possible
    arguments
       d (1,1) double {mustBePositive, mustBeInteger}
       k (1,1) double {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(k,d)}
       ranks (1,:) double {mustBePositive, mustBeInteger, ranksMustBeValid(d,k,ranks)} = [floor(d/k)*ones(1,k-1), d-(k-1)*floor(d/k)]
    end
    
    M = cell(1,k);
    U = RandomUnitary(d);
    for i = 1:k
        M{i} = zeros(d,d);
        for r = sum(ranks(1:i-1))+1:sum(ranks(1:i))
           M{i} = M{i} + U(:,r)*U(:,r)'; 
        end
    end
end

function ranksMustBeValid(d,k,ranks)
    if length(ranks) ~= k
       eid = 'Size:invalidRanks';
       msg = 'Must specify rank of each output.';
       throwAsCaller(MException(eid,msg))
    end
    if sum(ranks) ~= d
       eid = 'Size:invalidRanks';
       msg = 'Ranks must add to dimension.';
       throwAsCaller(MException(eid,msg))
    end
end
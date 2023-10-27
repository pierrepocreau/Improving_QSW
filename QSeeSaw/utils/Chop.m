function x = Chop(xin,tol)
%Chop Rounds numbers close to 0 to precisely 0
%   x = Chop(xin,tol) chops any real or imaginary parts smaller than tol
%   to be 0.
%   Default tolerance is 1e-9

    if nargin < 2
        tol = 1e-9;
    end
    
    realToChop = abs(real(xin)) < tol;
    imagToChop = abs(imag(xin)) < tol;

    x = xin - real(xin).*realToChop - 1i*imag(xin).*imagToChop;
end


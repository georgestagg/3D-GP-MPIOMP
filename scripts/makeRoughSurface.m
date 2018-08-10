function R = makeRoughSurface(gridx,gridy,l,sigma,number)
    function K = Kexpquad(posvec)
        K=sigma^2*exp(-((posvec(1)-posvec(3))^2+(posvec(2)-posvec(4))^2)/(2*l^2));
    end
    lx = length(gridx);
    ly = length(gridy);
    kargs = combvec(combvec(gridx,gridy),combvec(gridx,gridy));
    K = cellfun(@Kexpquad, num2cell(kargs, 1));
    K = reshape(K,lx*ly,lx*ly);
    mu = zeros(1,lx*ly);
    R = mvnrnd(mu,K,number);
    R = reshape(R,number,lx,ly);
end



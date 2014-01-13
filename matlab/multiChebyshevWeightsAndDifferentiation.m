function [x, w, D] = multiChebyshevWeightsAndDifferentiation(Nx, xMin, xMax, NIntervals)

D=zeros(Nx,Nx);
w=zeros(1,Nx);
x=zeros(1,Nx);

xMins=zeros(1,NIntervals);
xMaxs=zeros(1,NIntervals);

NxPerInterval=floor(Nx/NIntervals)+1;
firstNx = Nx - (NxPerInterval-1)*(NIntervals-1);
xMins(1) = xMin;
xMaxs(1) = xMin + (xMax-xMin)*(firstNx-1)/(Nx-1);
for interval=2:NIntervals
    xMins(interval) = xMaxs(interval-1);
    xMaxs(interval) = xMaxs(interval-1) + (xMax-xMin)*(NxPerInterval-1)/(Nx-1);
end


[x1, w1]=singleChebyshevWeights(firstNx, xMins(1), xMaxs(1));
[x1, D1]=singleChebyshevDifferentiation(firstNx, xMins(1), xMaxs(1));
w(1:firstNx) = w1;
x(1:firstNx) = x1;
D(1:(firstNx-1), 1:firstNx) = D1(1:(end-1),:);
D(firstNx, 1:firstNx) = D1(end,:)/2;

nextI=firstNx;
for interval=2:NIntervals
    [x1, w1]=singleChebyshevWeights(NxPerInterval, xMins(interval), xMaxs(interval));
    [x1, D1]=singleChebyshevDifferentiation(NxPerInterval, xMins(interval), xMaxs(interval));
    indices=nextI:(nextI+NxPerInterval-1);
    shortIndices=(nextI+1):(nextI+NxPerInterval-2);
    w(indices) = w(indices)+w1;
    x(indices) = x1;
    D(nextI, indices) = D(nextI, indices) + D1(1,:)/2;
    D(shortIndices, indices) = D1(2:(end-1), :);
    D(nextI+NxPerInterval-1, indices) = D1(end,:)/2;
    nextI=nextI+NxPerInterval-1;
end
a=size(D1, 2);
D(end, (end-a+1:end)) = D(end, (end-a+1:end)) + D1(end,:)/2;

%{
figure(1)
clf
plot(x,w,'o-')
%}

end

function [x, w]=singleChebyshevWeights(N1, a, b)
N=N1-1; bma=b-a;
c=zeros(N1,2);
c(1:2:N1,1)=(2./[1 1-(2:2:N).^2 ])'; c(2,2)=1;
f=real(ifft([c(1:N1,:);c(N:-1:2,:)]));
w=bma*([f(1,1); 2*f(2:N,1); f(N1,1)])/2;
x=0.5*((b+a)+N*bma*f(1:N1,2));
x=x(end:-1:1);
w=w(end:-1:1)';
end

function [x, D]=singleChebyshevDifferentiation(N, xMin, xMax)
N1=N-1;
if N1==0, D=0; x=1; return, end
x = cos(pi*(0:N1)/N1)';
c = [2; ones(N1-1,1); 2].*(-1).^(0:N1)';
X = repmat(x,1,N1+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N1+1)));  % off-diagonal entries
D = D - diag(sum(D'));             % diagonal entries
D = D * 2/(xMax-xMin);
x = (x+1) * (xMax-xMin)/2 + xMin;
D=fliplr(flipud(D));
x=fliplr(x');
end
function [R2Star,S0img] = createR2StarMap(recoInfo, FMImages, mask, varargin)
%createR2StarMap - Creates an R2* Map from multiple echo
%                   time data.
%
% Syntax:  [R2Star] = createFieldMap(rInfo, FMImages, mask);
%
% Inputs:
%    rInfo - recoInfo object that has been initialized with a Siemens
%    VB/VD/VE .dat file that can be read via mapVBVD
%    FMImages - SENSE Reconstructed multiple echo time images from Field
%    Map
%    mask - mask corresponding to one from SENSE map generation
%
% Outputs:
%    R2Star - R2Star map of size in rad/s [N,N,NSlices]
%    FMImages - SENSE reconstructed images [N,N,NSlices,NTEs]
%
% Example:
%    rInfo          = recoInfo(filename);
%    images         = gridCoilImages(rInfo);
%    [sen, mask]    = createSenMap(images,1);
%    [FM, FMImages] = createFieldMap(rInfo,sen,mask,1);
%    [R2Star]       = createR2StarMap(rInfo, FMImages, mask);
%
% Other m-files required: recoInfo.m, solve_pwls_pcg.m, Robject.m
% Subfunctions: none
% MAT-files required: none
%
% Author:
% University of Illinois at Urbana-Champaign
% email address:
% Website:
% September 2016; Last revision: 4-Sep-2016

%% Deal with input parsing using inputParser Class
p = inputParser;

p.addOptional('TEMask',@isscalar);

p.parse(varargin);

inputs = p.Results;

if isempty(p.Results.TEMask)
    TEMask = logical(ones(recoInfo.nEchoes,1));
else
    TEMask = logical(inputs.TEMask{2});
end


%% create field maps
SizeFMImages = size(FMImages);
S0img = zeros(SizeFMImages(1:3));
R2Star = zeros(SizeFMImages(1:3));
TEs = recoInfo.TE;
TEs = TEs - TEs(1); % Obtain relevant echo times for asymmetric spin echo
TEs = abs(TEs); %Negative TEs are equivalent to positive TEs for R2* purposes

sizes = size(FMImages);

f = @(x,xdata)(x(1).*exp(-x(2).*xdata));
options = optimoptions('lsqcurvefit','Display','off');
for ii = 1:sizes(1)
    for jj = 1:sizes(2)
        for kk = 1:sizes(3)
            if mask(ii,jj)
                sigvals = col(1E2*abs(FMImages(ii,jj,kk,TEMask)))./abs(max(col(FMImages)));
                %[result,resnorm] = lsqcurvefit(f,[0;50;],col(TEs([1,4])/1000)*1E-3,sigvals,[0,0],[1,1000],options);
                
                result = polyfit(col(TEs(TEMask)/1000),col(log(sigvals)),1);
                
                R2Star(ii,jj,kk) = -1*result(1);
                S0img(ii,jj,kk) = exp(result(2));
            end
        end
    end
    ii
end

end


